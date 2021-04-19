
// Uncomment this line for image exporting functionality
#define TC_IMAGE_IO

// Note: You DO NOT have to install taichi or taichi_mpm.
// You only need [taichi.h] - see below for instructions.
#include "../taichi/taichi.h"
#include <math.h>
#include <iostream>

using namespace taichi;
using namespace std;

using Vec = Vector2;
using Mat = Matrix2;

// Window
const int window_size = 800;

// Grid resolution (cells)
const int n = 80;

//number of particles per object
const int num_particles = 1000.f;

const real dt = 1e-4_f;
const real frame_dt = 1e-3_f;
const real dx = 1.0_f / n;
const real inv_dx = 1.0_f / dx;

const real alpha = 0.95;

// Snow material properties
const auto particle_mass = 1.0_f;
const auto vol = 1.0_f; // Particle Volume
const auto xi = 10.0_f; // Snow hardening factor
const auto E = 1e4_f;   // Young's Modulus
const auto nu = 0.2_f;  // Poisson ratio
const bool plastic = true;

// Initial Lamé parameters
const real mu_0 = E / (2 * (1 + nu));
const real lambda_0 = E * nu / ((1 + nu) * (1 - 2 * nu));

//neighbor grid
const int neighbor = 2;
int iteration = 0;

struct Particle
{
    // Position and velocity
    Vec x, v;
    //volume
    real vol;
    // Deformation gradient
    Mat F;
    // Affine momentum from APIC
    Mat C;
    // Determinant of the deformation gradient (i.e. volume)
    real Jp;
    // Color
    int c;
    Mat F_e; //elastic and plastic components of the deformation gradient
    Mat F_p;

    Particle(Vec x, int c, Vec v = Vec(0)) : x(x),
                                             v(v),
                                             vol(0),
                                             F(1),
                                             C(0),
                                             Jp(1),
                                             F_e(1),
                                             F_p(1),
                                             c(c)
    {
    }
};

real weight(real x)
{
    real abs_x = abs(x);
    // printf("abs x: %f\n", x);
    // if (abs_x < 1.f / n)
    if (abs_x < 1.f)

    {
        // printf("in return 1");
        return 0.5f * pow(abs_x, 3.0) - pow(abs_x, 2.0) + 2.0f / 3.0f;
    }
    if (abs_x < 2.f)
    {
        // printf("in return 2");
        return -1.0f / 6.0f * pow(abs_x, 3.0) + pow(abs_x, 2.0) - 2.0f * abs_x + 4.0f / 3.0f;
    }
    return 0;
}

real N_partial_derivative(real x)
{
    real abs_x = abs(x);
    if (abs_x < 1.f)
    {
        return 1.5f * pow(abs_x, 2) - 2 * abs_x;
    }
    if (abs_x < 2.f)
    {
        return -0.5f * pow(abs_x, 2) + 2 * abs_x - 2;
    }
    return 0;
}

//in: particle position, converted to grid bases (ie 1/h * (x_p - ih)...)
//out: delta w_i_p or delta w_i_p^T?
Vec weight_gradient(Vec pos)
{
    real N_x_real = N_partial_derivative(pos[0]) * weight(pos[1]);
    real N_y_real = weight(pos[0]) * N_partial_derivative(pos[1]);
    Vector2 vec(N_x_real, N_y_real);
    return vec;
}

//output: m * v_transpose = shape(1,2)
Vec multiply_vec_transpose(Mat m, Vec v)
{
    // printf("in vec: %f,%f\n", v[0], v[1]);
    // printf("in mat: %f,%f, %f, %f\n", m[0][0], m[0][1], m[1][0], m[1][1]);
    Vec vec = Vec(m[0][0] * v[0] + m[1][0] * v[0], m[0][1] * v[1] + m[1][1] * v[1]);
    return vec;
}

std::vector<Particle> particles;

// Vector3: [velocity_x, velocity_y, mass]
Vector3 grid[n + 1][n + 1];

void initialize()
{

    memset(grid, 0, sizeof(grid));

    //initialize particle weights and set mass of grid
    for (auto &p : particles)
    {

        // element-wise floor
        Vector2i base_coord = (p.x * inv_dx - Vec(0.5f)).cast<int>();

        //loop through neighboring grids [-2,2]
        //add weight * particle_mass to all neighboring grid cells
        for (int i = -2; i < 3; i++)
        {
            for (int j = -2; j < 3; j++)
            {
                Vector2i curr_grid = Vector2i(base_coord.x + i, base_coord.y + j);
                if (curr_grid.x >= 0 && curr_grid.x <= n && curr_grid.y >= 0 && curr_grid.y <= n)
                { //check bounds 0 to n in both dirs
                    //p.x = [0,1], base_coord = [0,n], fx is distance from particle position to nearest grid coordinate

                    Vec fx = p.x * inv_dx - curr_grid.cast<real>();
                    //compute weight

                    real N = weight(fx[0]) * weight(fx[1]);

                    // printf("N: %f\n", N);

                    //add mass to grid
                    grid[curr_grid[0]][curr_grid[1]].z += N * particle_mass;
                    // printf("mass: %f\n", grid[curr_grid[0]][curr_grid[1]].z);
                }
            }
        }
    }
    //initializing particle volumes based on density of grid
    for (auto &p : particles)
    {
        real density;
        Vector2i base_coord = (p.x * inv_dx - Vec(0.5f)).cast<int>();
        for (int i = -2; i < 3; i++)
        {
            for (int j = -2; j < 3; j++)
            {
                Vector2i curr_grid = Vector2i(base_coord.x + i, base_coord.y + j);
                if (curr_grid.x >= 0 && curr_grid.x <= n && curr_grid.y >= 0 && curr_grid.y <= n)
                { //check bounds 0 to n in both dirs
                    //p.x = [0,1], base_coord = [0,n], fx is distance from particle position to nearest grid coordinate
                    Vec fx = p.x * inv_dx - curr_grid.cast<real>();
                    //compute weight

                    real N = weight(fx[0]) * weight(fx[1]);

                    // h^3 = dx*dx*dx
                    // h^2 for 2d
                    // density is sum of grid masses multiplied by weight divided by vol/area of cell
                    density += grid[curr_grid.x][curr_grid.y].z * N / (dx * dx);
                    // printf("curr grid mass: %f\n", grid[curr_grid.x][curr_grid.y].z);
                }
            }
        }
        p.vol = particle_mass / density;
        // printf("mass, %f, density, %f\n", particle_mass, density);
        // printf("volume: %f\n", p.vol);
    }
}

void update(real dt)
{

    //store old grid velocities
    Vector2f oldVelocities[n + 1][n + 1];
    for (int i = 0; i <= n; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            auto &g = grid[i][j];
            oldVelocities[i][j] = Vector2f(g.x, g.y); //store old velocity
            // printf("old velocity: %f\n", oldVelocities[i][j].y);
            // if (g.z > 0.0f)
            // {                                                  //only if denominator is not 0
            //     g.x = g.x + dt * -1.0f * forces[i][j].x / g.z; //equation 10. update velocity (force is negative of sum in eq 6)
            //     g.y = g.y + dt * -1.0f * forces[i][j].y / g.z;
            // }
        }
    }

    // Reset grid
    std::memset(grid, 0, sizeof(grid));

    for (auto &p : particles)
    {
        //transfer mass
        // element-wise floor
        Vector2i base_coord = (p.x * inv_dx - Vec(0.5f)).cast<int>();

        //loop through neighboring grids [-2,2]
        //add weight * particle_mass to all neighboring grid cells
        for (int i = -2; i < 3; i++)
        {
            for (int j = -2; j < 3; j++)
            {
                Vector2i curr_grid = Vector2i(base_coord.x + i, base_coord.y + j);
                if (curr_grid.x >= 0 && curr_grid.x <= n && curr_grid.y >= 0 && curr_grid.y <= n)
                { //check bounds 0 to n in both dirs
                    //p.x = [0,1], base_coord = [0,n], fx is distance from particle position to nearest grid coordinate

                    Vec fx = p.x * inv_dx - curr_grid.cast<real>();
                    //compute weight

                    real N = weight(fx[0]) * weight(fx[1]);
                    //add mass to grid

                    grid[curr_grid[0]][curr_grid[1]].z += N * particle_mass;
                }
            }
        }
    }
    //step one of paper. transfering velocity
    for (auto &p : particles)
    {
        Vector2i base_coord = (p.x * inv_dx - Vec(0.5f)).cast<int>();

        //loop through neighboring grids [-2,2]
        //add velocity  to all neighboring grid cells
        for (int i = -2; i < 3; i++)
        {
            for (int j = -2; j < 3; j++)
            {
                Vector2i curr_grid = Vector2i(base_coord.x + i, base_coord.y + j);
                if (curr_grid.x >= 0 && curr_grid.x <= n && curr_grid.y >= 0 && curr_grid.y <= n)
                { //check bounds 0 to n in both dirs
                    Vec fx = p.x * inv_dx - curr_grid.cast<real>();
                    real N = weight(fx[0]) * weight(fx[1]);
                    //sum of particle's velocities
                    if (grid[curr_grid[0]][curr_grid[1]].z > 0.0f)
                    { //only if denominator is not 0
                        grid[curr_grid[0]][curr_grid[1]].x += N * p.v.x * particle_mass / grid[curr_grid[0]][curr_grid[1]].z;
                        grid[curr_grid[0]][curr_grid[1]].y += N * p.v.y * particle_mass / grid[curr_grid[0]][curr_grid[1]].z;
                    }
                }
            }
        }
    }

    //data structure to store grid forces
    Vector2f forces[n + 1][n + 1]; //same dimensions as grid
    for (int i = 0; i <= n; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            forces[i][j] = Vector2f(0.f, 0.f);
        }
    }
    //compute forces
    for (auto &p : particles)
    {
        //loop through neighbourhood [-2, 2]
        Vector2i base_coord = (p.x * inv_dx - Vec(0.5f)).cast<int>();

        for (int i = -2; i < 3; i++)
        {
            for (int j = -2; j < 3; j++)
            {
                Vector2i curr_grid = Vector2i(base_coord.x + i, base_coord.y + j);
                if (curr_grid.x >= 0 && curr_grid.x <= n && curr_grid.y >= 0 && curr_grid.y <= n)
                { //check bounds 0 to n in both dirs
                    Vec fx = p.x * inv_dx - curr_grid.cast<real>();
                    real V_p_n = determinant(p.F) * p.vol; //page 6
                    Mat F_hat_E_p = p.F_e;                 // equation 4 (x_hat = x_i)

                    Mat Re;
                    Mat Se;
                    polar_decomp(p.F_e, Re, Se);

                    real J_e = determinant(p.F_e);
                    real J = determinant(p.F);
                    real J_p = determinant(p.F_p);
                    Mat stress;
                    if (J != 0.0f)
                    {
                        // cout << "non zero J" << std::endl;
                        // printf("line 315");
                        real mu = mu_0 * exp(xi * (1 - J_p));
                        // if (mu != 0.0f)
                        // {
                        //     cout << "non zero mu" << std::endl;
                        // }
                        real lambda = lambda_0 * exp(xi * (1 - J_p));
                        // if (lambda != 0.0f)
                        // {
                        //     cout << "non zero lambda" << std::endl;
                        // }
                        // if (determinant(p.F_e) != 0.0f)
                        // {
                        //     cout << "non zero P fe" << std::endl;
                        // }
                        // cout << "mu: " << mu << ", J: " << J << ", p.F_e: " << p.F_e << ", Re: " << Re << ", lambda: " << lambda << ", J_e: " << J_e << endl;
                        // if (determinant(p.F_e - Re) != 0.0f)
                        // {
                        //     cout << "non zero difference!" << endl; //TODO!! this is always zero
                        // }
                        // if (J_e - 1 != 0.0f)
                        // {
                        //     cout << "non zero difference!" << endl; //TODO!! this is always zero too
                        // }
                        // cout << p.F_e - Re << endl;
                        stress = ((2.0f * mu) / J) * (p.F_e - Re) * transposed(p.F_e) + (lambda * 1.0f / J) * (J_e - 1) * (J_e * Mat(1)); // above quation 6
                        // stress = ((2.0f * mu) / J) * (p.F_e) * transposed(p.F_e) + (lambda * 1.0f / J) * (J_e - 1) * (J_e * Mat(1)); // above quation 6

                        // cout << "stress: " << stress << endl;
                    }
                    else
                    {
                        stress = Mat(0.0f);
                    }
                    if (determinant(stress) != 0)
                    {
                        // cout << "non zero stress: " << stress << std::endl;
                    }

                    Vec N = weight_gradient(fx);
                    // printf("curr coords, %d,%d", curr_grid.x, curr_grid.y);
                    Vector2f force_at_grid_by_particle = V_p_n * multiply_vec_transpose(stress, N); // equation 6
                    // force_at_grid_by_particle = Vector2f(-0.1f, 0.1f);

                    // printf("stress*n: %f, %f\n", multiply_vec_transpose(stress, N)[0], multiply_vec_transpose(stress, N)[1]);
                    // printf("Vpn: %f\n", V_p_n);
                    real epsilon = 1e-4f;
                    // if (abs(force_at_grid_by_particle[1]) > epsilon)
                    // {
                    //     printf("force: %f, %f\n", force_at_grid_by_particle[0], force_at_grid_by_particle[1]);
                    // }
                    forces[curr_grid.x][curr_grid.y] += force_at_grid_by_particle;
                }
            }
        }
    }

    //Step 4
    for (int i = 0; i <= n; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            if (grid[i][j][2] > 0.f)
            {
                // cout << "grid forces: " << forces[i][j][0] << ", " << forces[i][j][1] << std::endl;

                grid[i][j][0] += forces[i][j][0] * dt; //(1.0f / grid[i][j][2]) * dt;
                grid[i][j][1] += forces[i][j][1] * dt; //(1.0f / grid[i][j][2]) * dt;
                // if (forces[i][j][1] > 0.f)
                // {
                // std::cout << forces[i][j][1] << std::endl;
                // }
            }
        }
    }

    // For all grid nodes: GRAVITY
    for (int i = 0; i <= n; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            auto &g = grid[i][j];
            // No need for epsilon here
            if (g[2] > 0)
            {
                // Normalize by mass
                // g /= g[2];
                // Gravity
                g += dt * Vector3(0, -200, 0);

                //copied from taichi example
                // boundary thickness
                real boundary = 0.05;
                // Node coordinates
                real x = (real)i / n;
                real y = real(j) / n;

                //added a sphere collider (hemisphere)
                Vec circleCenter = Vec(0.5, 0 + boundary);
                real circleRadius = 0.1;
                real mu = 0.1;

                //if inside the sphere...
                if ((x - circleCenter.x) * (x - circleCenter.x) + (y - circleCenter.y) * (y - circleCenter.y) < circleRadius * circleRadius)
                {
                    Vec n = normalized(Vec(x, y) - circleCenter);
                    Vec v = Vec(g.x, g.y);
                    real v_dot_n = v.dot(n);
                    if (v_dot_n < 0)
                    { //section 8 body collision
                        Vec v_t = v - n * v_dot_n;
                        real v_t_norm = pow(v_t.dot(v_t), 0.5);
                        if (v_t_norm > 0)
                        {
                            Vec v_prime = v_t + mu * v_dot_n * v_t / v_t_norm;
                            g.x = v_prime.x;
                            g.y = v_prime.y;
                        }
                    }
                }
                // Sticky boundary
                if (x < boundary || x > 1 - boundary || y > 1 - boundary)
                {
                    g = Vector3(0);
                }
                // Separate boundary
                if (y < boundary)
                {
                    g[1] = std::max(0.0f, g[1]);
                }
            }
            //copy end
        }
    }

    for (auto &p : particles)
    {
        Vector2i base_coord = (p.x * inv_dx - Vec(0.5f)).cast<int>();

        Mat v_p_n_plus_1;
        for (int i = -2; i < 3; i++)
        {
            for (int j = -2; j < 3; j++)
            {
                Vector2i curr_grid = Vector2i(base_coord.x + i, base_coord.y + j);
                if (curr_grid.x >= 0 && curr_grid.x <= n && curr_grid.y >= 0 && curr_grid.y <= n)
                { //check bounds 0 to n in both dirs
                    Vec fx = p.x * inv_dx - curr_grid.cast<real>();
                    Vec grid_velocity(grid[curr_grid.x][curr_grid.y].x, grid[curr_grid.x][curr_grid.y].y);
                    v_p_n_plus_1 += Mat::outer_product(grid_velocity, weight_gradient(fx));
                }
            }
        }

        //update force

        // std::cout << v_p_n_plus_1 << std::endl;
        //update elastic compoenent - before we do plasticity the elastic component gets all of the F
        p.F_e = (Mat(1) + dt * v_p_n_plus_1) * p.F_e;
        // cout << p.F_e << endl;
        // p.F_e = (Mat(1)) * p.F_e;
        // printf("p.fe Before: %f\n", determinant(p.F_e));

        //plastic component - compiles but is does not change sim.
        Mat F_hat_P_p = p.F_p; // Section 7
        Mat U_p, Sig_hat_p, V_transpose_p;
        svd(p.F_e, U_p, Sig_hat_p, V_transpose_p); //compute singular value decomposition
        real ten = 10;
        real theta_c = 2.5f * pow(ten, -2.0); //move up later
        real theta_s = 7.5f * pow(ten, -3.0);
        // real theta_c = 0.f; //move up later
        // real theta_s = 0.f;

        Mat Sig_p;
        Sig_p = Sig_hat_p;
        // std::cout << "409: " << (Sig_p) << std::endl;

        Sig_p[0][0] = clamp(Sig_p[0][0], 1 - theta_c, 1 + theta_s); //clamp
        Sig_p[1][1] = clamp(Sig_p[1][1], 1 - theta_c, 1 + theta_s);

        // std::cout << "414: " << (Sig_p) << std::endl;
        p.F = p.F_e * p.F_p;

        p.F_e = U_p * Sig_p * V_transpose_p;
        p.F_p = transposed(V_transpose_p) * inversed(Sig_p) * transposed(U_p) * p.F; //

        // assert(p.F == p.F_e * p.F_p); //TODO

        // p.F = p.F_p;
        // p.F = Mat(500);
        // printf("p.fe After: %f\n", determinant(p.F_e));
        // printf("p.fp After: %f\n", determinant(p.F_p));

        // printf("p.fe: %f\n", determinant(p.F));
        // std::cout << p.F_p << std::endl;

        //update particle velocities
        Vec v_PIC(0, 0);
        Vec v_FLIP = p.v;
        // printf("initial flip: %f\n", p.v[1]);
        for (int i = -2; i < 3; i++)
        {
            for (int j = -2; j < 3; j++)
            {

                Vector2i curr_grid = Vector2i(base_coord.x + i, base_coord.y + j);
                if (curr_grid.x >= 0 && curr_grid.x <= n && curr_grid.y >= 0 && curr_grid.y <= n)
                { //check bounds 0 to n in both dirs
                    Vec fx = p.x * inv_dx - curr_grid.cast<real>();
                    real N = weight(fx[0]) * weight(fx[1]);
                    //update PIC velocity
                    v_PIC.x += grid[curr_grid.x][curr_grid.y].x * N;
                    v_PIC.y += grid[curr_grid.x][curr_grid.y].y * N;

                    //update FLIP velocity
                    v_FLIP.x += (grid[curr_grid.x][curr_grid.y].x - oldVelocities[curr_grid.x][curr_grid.y].x) * N;
                    v_FLIP.y += (grid[curr_grid.x][curr_grid.y].y - oldVelocities[curr_grid.x][curr_grid.y].y) * N;
                }
            }
        }
        real epsilon = 1e-4;
        // if (abs(v_PIC.y) > epsilon)
        // {
        //     printf("pic y: %f\n", v_PIC.y);
        // }
        // if (abs(v_FLIP.y) > epsilon)
        // {
        //     printf("flip y: %f\n", v_FLIP.y);
        // }

        //update particle velocities
        // p.v = (1 - alpha) * v_PIC + alpha * v_FLIP;
        p.v = v_PIC;
        // std::cout << p.v << std::endl;
        // printf("P v: %f, %f\n", p.v[0], p.v[1]);

        //update particle positions
        p.x += p.v * dt;
    }
}
// Seed particles with position and color
void add_object(Vec center, int c)
{
    // Randomly sample num_particles particles in the square
    for (int i = 0; i < num_particles; i++)
    {
        particles.push_back(Particle((Vec::rand() * 2.0f - Vec(1)) * 0.08f + center, c));
    }
}

// from: https://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c
std::vector<std::string>
getNextLineAndSplitIntoTokens(std::istream &str)
{
    std::vector<std::string> result;
    std::string line;
    std::getline(str, line);

    std::stringstream lineStream(line);
    std::string cell;

    while (std::getline(lineStream, cell, ','))
    {
        result.push_back(cell);
    }
    // This checks for a trailing comma with no data after it.
    if (!lineStream && cell.empty())
    {
        // If there was a trailing comma then add an empty element.
        result.push_back("");
    }
    return result;
}

void add_from_csv(char *infile_path, Vec center, int c)
{
    std::ifstream infile(infile_path);

    std::vector<std::string> split_line = getNextLineAndSplitIntoTokens(infile);
    split_line = getNextLineAndSplitIntoTokens(infile); //read two lines because first line is labels
    while (split_line.size() > 1)
    {
        real x_pos = std::stof(split_line[0].c_str()) / 50.f;
        real y_pos = std::stof(split_line[1].c_str()) / 50.f - 0.5f;                 //change to split_line[2] to get a bird's eye view cow :o
        particles.push_back(Particle(Vec(x_pos + center[0], y_pos + center[1]), c)); //using x and y coordinates for 2d
        split_line = getNextLineAndSplitIntoTokens(infile);
    }
}

int main(int argc, char *argv[])
{

    GUI gui("Real-time 2D MLS-MPM", window_size, window_size);
    auto &canvas = gui.get_canvas();

    if (argc == 1) //default
    {
        add_object(Vec(0.55, 0.45), 0xFFFAFA);
        add_object(Vec(0.45, 0.65), 0xFFFAFA);
        add_object(Vec(0.55, 0.85), 0xFFFAFA);
    }
    else
    {
        assert(argc == 2);
        add_from_csv(argv[1], Vec(0.55, 0.85), 0xF2B134);
    }

    int frame = 0;

    //initialize particle values
    initialize();

    // Main Loop
    for (int step = 0;; step++)
    {
        // Advance simulation
        update(dt);

        // Visualize frame
        if (step % int(frame_dt / dt) == 0)
        {
            // Clear background
            canvas.clear(0x112F41);
            // Box
            canvas.rect(Vec(0.04), Vec(0.96)).radius(2).color(0x4FB99F).close();
            // Particles
            for (auto p : particles)
            {
                canvas.circle(p.x).radius(2).color(p.c);
            }
            // Update image
            gui.update();

            // Write to disk (optional)
            // canvas.img.write_as_image(fmt::format("tmp/{:05d}.png", frame++));
        }
    }
}

/* -----------------------------------------------------------------------------
** Reference: A Moving Least Squares Material Point Method with Displacement
              Discontinuity and Two-Way Rigid Body Coupling (SIGGRAPH 2018)
  By Yuanming Hu (who also wrote this 88-line version), Yu Fang, Ziheng Ge,
           Ziyin Qu, Yixin Zhu, Andre Pradhana, Chenfanfu Jiang
** Build Instructions:
Step 1: Download and unzip mls-mpm88.zip (Link: http://bit.ly/mls-mpm88)
        Now you should have "mls-mpm88.cpp" and "taichi.h".
Step 2: Compile and run
* Linux:
    g++ mls-mpm88-explained.cpp -std=c++14 -g -lX11 -lpthread -O3 -o mls-mpm
    ./mls-mpm
* Windows (MinGW):
    g++ mls-mpm88-explained.cpp -std=c++14 -lgdi32 -lpthread -O3 -o mls-mpm
    .\mls-mpm.exe
* Windows (Visual Studio 2017+):
  - Create an "Empty Project"
  - Use taichi.h as the only header, and mls-mpm88-explained.cpp as the only source
  - Change configuration to "Release" and "x64"
  - Press F5 to compile and run
* OS X:
    g++ mls-mpm88-explained.cpp -std=c++14 -framework Cocoa -lpthread -O3 -o mls-mpm
    ./mls-mpm
** FAQ:
Q1: What does "1e-4_f" mean?
A1: The same as 1e-4f, a float precision real number.
Q2: What is "real"?
A2: real = float in this file.
Q3: What are the hex numbers like 0xED553B?
A3: They are RGB color values.
    The color scheme is borrowed from
    https://color.adobe.com/Copy-of-Copy-of-Core-color-theme-11449181/
Q4: How can I get higher-quality?
A4: Change n to 320; Change dt to 1e-5; Change E to 2e4;
    Change particle per cube from 500 to 8000 (Ln 72).
    After the change the whole animation takes ~3 minutes on my computer.
Q5: How to record the animation?
A5: Uncomment Ln 2 and 85 and create a folder named "tmp".
    The frames will be saved to "tmp/XXXXX.png".
    To get a video, you can use ffmpeg. If you already have taichi installed,
    you can simply go to the "tmp" folder and execute
      ti video 60
    where 60 stands for 60 FPS. A file named "video.mp4" is what you want.
Q6: How is taichi.h generated?
A6: Please check out my #include <taichi> talk:
    http://taichi.graphics/wp-content/uploads/2018/11/include_taichi.pdf
    and the generation script:
    https://github.com/yuanming-hu/taichi/blob/master/misc/amalgamate.py
    You can regenerate it using `ti amal`, if you have taichi installed.
Questions go to yuanming _at_ mit.edu
                            or https://github.com/yuanming-hu/taichi_mpm/issues.
                                                      Last Update: March 6, 2019
                                                      Version 1.5
----------------------------------------------------------------------------- */
