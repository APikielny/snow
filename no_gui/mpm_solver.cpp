#include "mpm_solver.h"
//#include "../taichi/taichi.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include "Eigen/SVD"


using namespace Eigen;
using namespace std;

bool sphere_collision = true;
double force_factor = 10.0f;
bool render = false;
bool write = false;



static double weight(double x)
{
    double abs_x = abs(x);
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

static double N_partial_derivative(double x)
{
    double abs_x = abs(x);
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
static Vec weight_gradient(Vec pos)
{
    double N_x_double = N_partial_derivative(pos[0]) * weight(pos[1]) * weight(pos[2]);
    double N_y_double = weight(pos[0]) * N_partial_derivative(pos[1]) * weight(pos[2]);
    double N_z_double = weight(pos[0]) * weight(pos[1]) * N_partial_derivative(pos[2]);
    Vec vec(N_x_double, N_y_double, N_z_double);
    return vec;
}

static void polar_decomp(const Mat m, Mat& R, Mat& S){
//    JacobiSVD<MatrixXd> svd(*m, ComputeThinU | ComputeThinV); //compute singular value decomposition
//    Mat U = -1.0*svd.matrixU();
//    Mat Sigma = svd.singularValues().asDiagonal();
//    Sigma.col(1)[1] = -1.0*Sigma.col(1)[1];
//    Mat V_transpose = svd.matrixV().transpose();
//    V_transpose.col(0) = -1.0*V_transpose.col(0);


//        double x = m(0, 0) + m(1, 1);
//        double y = m(1, 0) - m(0, 1);
//        double scale = 1.0 / sqrt(x * x + y * y);
//        double c = x * scale, s = y * scale;
//        R(0, 0) = c;
//        R(0, 1) = -s;
//        R(1, 0) = s;
//        R(1, 1) = c;
//        S = R.transpose() * m;

    //see relation of SVD to polar decomp: https://en.wikipedia.org/wiki/Polar_decomposition
    JacobiSVD<MatrixXd> svd(m, ComputeThinU | ComputeThinV); //compute singular value decomposition
    Matrix3d singular_value_matrix;
    singular_value_matrix.diagonal()<<svd.singularValues()[0], svd.singularValues()[1], svd.singularValues()[2];
    R = svd.matrixV() * singular_value_matrix * svd.matrixV().transpose();
    S = svd.matrixU() * svd.matrixV().transpose();

}

void collide(Vec coords, Vec &velocity){
    // Vec out_vel(0.0f, 0.0f);

        //copied from taichi example
        // boundary thickness
        double boundary = 0.25;

        //added a sphere collider (hemisphere)
        //initial z coordinate to be z coordinate of the value that was passed in
        Vec circleCenter = Vec(0.5, 0 + boundary, coords.z());
        double circleRadius = 0.1f;
        double mu = 0.1;

        //if inside the sphere...
        if (sphere_collision)
        {
            if ((coords.x() - circleCenter.x()) * (coords.x() - circleCenter.x()) + (coords.y() - circleCenter.y()) * (coords.y() - circleCenter.y()) < circleRadius * circleRadius)
            {
                Vec n = (coords - circleCenter).normalized();
                double v_dot_n = velocity.dot(n);
                if (v_dot_n < 0)
                { //section 8 body collision
                    Vec v_t = velocity - n * v_dot_n;
                    double v_t_norm = pow(v_t.dot(v_t), 0.5);
                    if (v_t_norm > 0)
                    {
                        Vec v_prime = v_t + mu * v_dot_n * v_t / v_t_norm;
                        velocity = v_prime;
                    }
                }
            }
        }

//        // Sticky boundary

//        if (coords.x < boundary || coords.x > 1 - boundary || coords.y > 1 - boundary)
//        {
//            velocity = Vec(0.f, 0.f);
//        }
        // Separate boundary
        if (coords.y() < boundary)
        {
            velocity[0] = 0.0f;
            velocity[1] = std::max(0.0f, (float)velocity.y());
            velocity[2] = 0.0f;
        }
}

mpm_solver::mpm_solver()
{
}





void mpm_solver::initialize()
{

//    std::cout << "initializing\n" << std::endl;

    //create shapes
    add_object(Vec(0.55, 0.45, 0.5f), 0xFFFAFA);
    add_object(Vec(0.45, 0.65, 0.5f), 0xFFFAFA);
    add_object(Vec(0.55, 0.85, 0.5f), 0xFFFAFA);

    std::cout << "Number of initialized particles: " << particles.size() << std::endl;

    memset(grid, 0, sizeof(grid));

    /**************************************/
    /******** INITIALIZE GRID MASS ********/
    /**************************************/
    for (auto &p : particles)
    {
        Vector3i base_coord = (p.x * inv_dx - Vec(0.5, 0.5, 0.5)).cast<int>();

        //loop through neighboring grids [-2,2]
        //add weight * particle_mass to all neighboring grid cells
        for (int i = -neighbor; i <= neighbor; i++)
        {
            for (int j = -neighbor; j <= neighbor; j++)
            {
                for (int k = -neighbor; k <= neighbor; k++)
                {
                    Vector3i curr_grid = Vector3i(base_coord.x() + i, base_coord.y() + j, base_coord.z() + k);
                    if (curr_grid.x() >= 0 && curr_grid.x() <= n && curr_grid.y() >= 0 && curr_grid.y() <= n  && curr_grid.z() >= 0 && curr_grid.z() <= n)
                    { //check bounds 0 to n in both dirs
                        //p.x = [0,1], base_coord = [0,n], fx is distance from particle position to nearest grid coordinate

                        Vec fx = p.x * inv_dx - curr_grid.cast<double>();

                        //compute weight
                        double N = weight(fx[0]) * weight(fx[1]) * weight(fx[2]);

                        //add mass to grid
                        grid[curr_grid[0]][curr_grid[1]][curr_grid[2]].w() += N * particle_mass;    //changed to w

                    }
                }
            }
        }
    }


    /**************************************/
    /***** INITIALIZE PARTICLE VOLUME *****/
    /**************************************/
    for (auto &p : particles)
    {
        double density = 0.0f;
        Vector3i base_coord = (p.x * inv_dx - Vec(0.5, 0.5, 0.5)).cast<int>();
        for (int i = -neighbor; i <= neighbor; i++)
        {
            for (int j = -neighbor; j <= neighbor; j++)
            {
                for (int k = -neighbor; k <= neighbor; k++)
                {
                    Vector3i curr_grid = Vector3i(base_coord.x() + i, base_coord.y() + j, base_coord.z() + k);
                    if (curr_grid.x() >= 0 && curr_grid.x() <= n && curr_grid.y() >= 0 && curr_grid.y() <= n && curr_grid.z() >= 0 && curr_grid.z() <= n)
                    { //check bounds 0 to n in both dirs
                        //p.x = [0,1], base_coord = [0,n], fx is distance from particle position to nearest grid coordinate
                        Vec fx = p.x * inv_dx - curr_grid.cast<double>();

                        //compute weight
                        double N = weight(fx[0]) * weight(fx[1]) * weight(fx[2]);

                        // h^3 = dx*dx*dx
                        // h^2 for 2d
                        // density is sum of grid masses multiplied by weight divided by vol/area of cell
                        density += grid[curr_grid.x()][curr_grid.y()][curr_grid.z()].w() * N; // / (dx * dx * dx);
                    }
                }
            }
        }
        p.vol = particle_mass / density;
        // printf("mass, %f, density, %f\n", particle_mass, density);
        // printf("volume: %f\n", p.vol);
    }


    /**************************************/
    /************ RENNDERING **************/
    /**************************************/
    if(render){
        // std::vector<Vector3f> vertices;
        // std::vector<Vector3i> faces;
        // int j = 0;
        // //vertices of tetrahedron centered at origin with side length ~=offset
        // float offset = 0.01;
        // Vector3f v_1 = Vector3f(offset, offset, offset);
        // Vector3f v_2 = Vector3f(-offset, -offset, offset);
        // Vector3f v_3 = Vector3f(-offset, offset, -offset);
        // Vector3f v_4 = Vector3f(offset, -offset, -offset);

        // for(int i = 0; i < particles.size(); i++){
        //     Vec p = particles[i].x;
        //     Vector3f pos = Vector3f(5*(float)p.x(), 5*(float)p.y(), 5*(float)p.z());

        //     vertices.push_back(v_1+pos);
        //     vertices.push_back(v_2+pos);
        //     vertices.push_back(v_3+pos);
        //     vertices.push_back(v_4+pos);

        //     Vector4i t = Vector4i(j, j+1, j+2, j+3);
        //     Vector3i f_1 = Vector3i(t[3], t[1], t[2]); //opposite t[0]
        //     Vector3i f_2 = Vector3i(t[2], t[0], t[3]); //opposite t[1]
        //     Vector3i f_3 = Vector3i(t[3], t[0], t[1]); //opposite t[2]
        //     Vector3i f_4 = Vector3i(t[1], t[0], t[2]); //opposite t[3]

        //     faces.emplace_back(f_1);
        //     faces.emplace_back(f_2);
        //     faces.emplace_back(f_3);
        //     faces.emplace_back(f_4);

        //     j+=4;
        // }
        // m_shape.init(vertices, faces);
    }
}

void mpm_solver::update(double dt)
{

//    std::cout<<"update"<<std::endl;

    // Reset grid
    std::memset(grid, 0, sizeof(grid));


    /**************************************/
    /******* TRANSFER MASS TO GRID ********/
    /**************************************/
    for (auto &p : particles)
    {
        //transfer mass
        // element-wise floor
        Vector3i base_coord = (p.x * inv_dx - Vec(0.5, 0.5, 0.5)).cast<int>();

        //loop through neighboring grids [-2,2]
        //add weight * particle_mass to all neighboring grid cells
        for (int i = -neighbor; i <= neighbor; i++)
        {
            for (int j = -neighbor; j <= neighbor; j++)
            {
                for (int k = -neighbor; k <= neighbor; k++)
                {
                    Vector3i curr_grid = Vector3i(base_coord.x() + i, base_coord.y() + j, base_coord.z() + k);
                    if (curr_grid.x() >= 0 && curr_grid.x() <= n && curr_grid.y() >= 0 && curr_grid.y() <= n && curr_grid.z() >= 0 && curr_grid.z() <= n)
                    { //check bounds 0 to n in both dirs
                        //p.x = [0,1], base_coord = [0,n], fx is distance from particle position to nearest grid coordinate

                        Vec fx = p.x * inv_dx - curr_grid.cast<double>();

                        //compute weight
                        double N = weight(fx[0]) * weight(fx[1]) * weight(fx[2]);

                        //add mass to grid
                        grid[curr_grid[0]][curr_grid[1]][curr_grid[2]].w() += N * particle_mass;    //changed to w
                    }
                }
            }
        }
    }
//    std::cout << "transferred mass\n" << std::endl;



    /**************************************/
    /******* TRANSFER VEL TO GRID ********/
    /**************************************/

    //step one of paper. transfering velocity
    for (auto &p : particles)
    {
        Vector3i base_coord = (p.x * inv_dx - Vec(0.5, 0.5, 0.5)).cast<int>();

//        add velocity  to all neighboring grid cells
        for (int i = -neighbor; i <= neighbor; i++)
        {
            for (int j = -neighbor; j <= neighbor; j++)
            {
                for (int k = -neighbor; k <= neighbor; k++)
                {
                    Vector3i curr_grid = Vector3i(base_coord.x() + i, base_coord.y() + j, base_coord.z() + k);
                    if (curr_grid.x() >= 0 && curr_grid.x() <= n && curr_grid.y() >= 0 && curr_grid.y() <= n && curr_grid.z() >= 0 && curr_grid.z() <= n)
                    { //check bounds 0 to n in both dirs
                        Vec fx = p.x * inv_dx - curr_grid.cast<double>();
                        double N = weight(fx[0]) * weight(fx[1]) * weight(fx[2]);
                        //sum of particle's velocities
                        if (grid[curr_grid[0]][curr_grid[1]][curr_grid[2]].w() > 0.0f)  //changed to w
                        { //only if denominator is not 0
                            grid[curr_grid[0]][curr_grid[1]][curr_grid[2]].x() += N * p.v.x() * particle_mass / grid[curr_grid[0]][curr_grid[1]][curr_grid[2]].w(); //changed to w
                            grid[curr_grid[0]][curr_grid[1]][curr_grid[2]].y() += N * p.v.y() * particle_mass / grid[curr_grid[0]][curr_grid[1]][curr_grid[2]].w();
                            grid[curr_grid[0]][curr_grid[1]][curr_grid[2]].z() += N * p.v.z() * particle_mass / grid[curr_grid[0]][curr_grid[1]][curr_grid[2]].w();
                        }
                    }
                }
            }
        }
    }
//    std::cout << "step 1\n" << std::endl;



    /**************************************/
    /******** STORE OLD VELOCITIES ********/
    /**************************************/
    Vec oldVelocities[n + 1][n + 1][n + 1];
    for (int i = 0; i <= n; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            for (int k = 0; k <= n; k++)
            {
                auto &g = grid[i][j][k];
                oldVelocities[i][j][k] = Vec(g.x(), g.y(), g.z()); //store old velocity
            }
        }
    }


    //data structure to store grid forces
    Vec forces[n + 1][n + 1][n + 1]; //same dimensions as grid
    std::memset(forces, 0, sizeof(forces)); //should be the same as looping

//    for (int i = 0; i <= n; i++)
//    {
//        for (int j = 0; j <= n; j++)
//        {
//            for (int k = 0; k <= n; k++)
//            {
//                forces[i][j][k] = Vec(0.f, 0.f, 0.f);
//            }
//        }
//    }


    /**************************************/
    /********* COMPUTE GRID FORCES ********/
    /**************************************/
    for (auto &p : particles)
    {

        Vector3i base_coord = (p.x * inv_dx - Vec(0.5, 0.5, 0.5)).cast<int>();

        double V_p_n = p.vol; //page 6

        Mat Re;
        Mat Se;
        polar_decomp(p.F_e, Re, Se);

        double J_e = p.F_e.determinant();
        double J = p.F.determinant();
        double J_p = p.F_p.determinant();
        Mat stress;

        double mu = mu_0 * exp(xi * (1.0f - J_p));

        double lambda = lambda_0 * exp(xi * (1.0f - J_p));

        stress = (2.0f * mu) * (p.F_e - Re) * p.F_e.transpose() + lambda * (J_e - 1) * (J_e * Mat(MatrixXd::Identity(3, 3))); // above equation 6

        //loop through neighbourhood [-2, 2]
        for (int i = -neighbor; i <= neighbor; i++)
        {
            for (int j = -neighbor; j <= neighbor; j++)
            {
                for (int k = -neighbor; k <= neighbor; k++)
                {
                    Vector3i curr_grid = Vector3i(base_coord.x() + i, base_coord.y() + j, base_coord.z() + k);
                    if (curr_grid.x() >= 0 && curr_grid.x() <= n && curr_grid.y() >= 0 && curr_grid.y() <= n && curr_grid.z() >= 0 && curr_grid.z() <= n)
                    { //check bounds 0 to n in both dirs
                        Vec fx = p.x * inv_dx - curr_grid.cast<double>();
                        Vec N = weight_gradient(fx);
                        // printf("curr coords, %d,%d", curr_grid.x, curr_grid.y);
                        Vec force_at_grid_by_particle = Vec(stress*N)*V_p_n; // equation 6

                        forces[curr_grid.x()][curr_grid.y()][curr_grid.z()] += force_at_grid_by_particle;
                    }
                }
            }
        }
    }
//    std::cout << "compute forces\n" << std::endl;




    /**************************************/
    /******* UPDATE GRID VELOCITIES *******/
    /**************************************/

    //Step 4
    for (int i = 0; i <= n; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            for (int k = 0; k <= n; k++)
            {
                if (grid[i][j][k].w() > 0.f)
                {
                    grid[i][j][k].x() += force_factor*forces[i][j][k].x() * (1.0f / grid[i][j][k].w()) * dt;
                    grid[i][j][k].y() += force_factor*forces[i][j][k].y() * (1.0f / grid[i][j][k].w()) * dt;
                    grid[i][j][k].z() += force_factor*forces[i][j][k].z() * (1.0f / grid[i][j][k].w()) * dt;
                }
            }
        }
    }

    /**************************************/
    /********* ADD GRAVITY TO GRID ********/
    /**************************************/
    for (int i = 0; i <= n; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            for (int k = 0; k <= n; k++)
            {
                auto &g = grid[i][j][k];
                if (g[3] > 0)   //changed to w
                {
                    // g /= g[2];
                    // Gravity
                    g += dt * Vector4d(0, -20000, 0, 0);   //not sure

                    // Node coordinates
                    double x = double(i) / n;
                    double y = double(j) / n;
                    double z = double(k) / n;

                    Vec out_vel = Vec(g.x(), g.y(), g.z());
                    collide(Vec(x, y, z), out_vel);

                    g.x() = out_vel.x();
                    g.y() = out_vel.y();
                    g.z() = out_vel.z();
                }
            }
        }
    }
//   std::cout << "added gravity\n" << std::endl;



    /**************************************/
    /**** UPDATE DEFORMATION GRADIENT *****/
    /**************************************/
    for (auto &p : particles)
    {
        Vector3i base_coord = (p.x * inv_dx - Vec(0.5, 0.5, 0.5)).cast<int>();

        Mat v_p_n_plus_1 = Matrix3d::Zero();
        for (int i = -neighbor; i <= neighbor; i++)
        {
            for (int j = -neighbor; j <= neighbor; j++)
            {
                for (int k = -neighbor; k <= neighbor; k++)
                {
                    Vector3i curr_grid = Vector3i(base_coord.x() + i, base_coord.y() + j, base_coord.z() + k);
                    if (curr_grid.x() >= 0 && curr_grid.x() <= n && curr_grid.y() >= 0 && curr_grid.y() <= n && curr_grid.z() >= 0 && curr_grid.z() <= n)
                    { //check bounds 0 to n in both dirs

                        Vec fx = p.x * inv_dx - curr_grid.cast<double>();
                        Vec grid_velocity(grid[curr_grid.x()][curr_grid.y()][curr_grid.z()].x(), grid[curr_grid.x()][curr_grid.y()][curr_grid.z()].y(), grid[curr_grid.x()][curr_grid.y()][curr_grid.z()].z());
                        v_p_n_plus_1 += grid_velocity*(weight_gradient(fx).transpose());
                    }
                }
            }
        }

//        update elastic compoenent - before we do plasticity the elastic component gets all of the F
        p.F_e = (Mat(MatrixXd::Identity(3, 3)) + dt * v_p_n_plus_1) * p.F_e;

        //plastic component - compiles but is does not change sim.
        JacobiSVD<MatrixXd> svd(p.F_e, ComputeThinU | ComputeThinV); //compute singular value decomposition
        Mat U_p = svd.matrixU();
        Mat Sig_hat_p = svd.singularValues().asDiagonal();
        Mat V_transpose_p = svd.matrixV().transpose();

//        //ERROR CHECKING SVD - https://stackoverflow.com/questions/34392115/eigen-jacobi-svd
//        MatrixXd Cp = U_p * Sig_hat_p * V_transpose_p;
//        MatrixXd diff = Cp - p.F_e;
//        cout << "ERROR CHECKING SVD: " << diff.array().abs().sum() << "\n";


        double theta_c = 2.5f * pow(10, -2.0); //move up later
        double theta_s = 7.5f * pow(10, -3.0);


        Mat Sig_p;
        Sig_p = Sig_hat_p;

        // Sig_p.col(0)[0] = clamp(Sig_p.col(0)[0], 1 - theta_c, 1 + theta_s); //clamp
        // Sig_p.col(1)[1] = clamp(Sig_p.col(1)[1], 1 - theta_c, 1 + theta_s);
        Sig_p.col(0)[0] = std::max(Sig_p.col(0)[0],1 - theta_c);
        Sig_p.col(0)[0] = std::min(Sig_p.col(0)[0],1 + theta_s);
        Sig_p.col(1)[1] = std::max(Sig_p.col(1)[1],1 - theta_c);
        Sig_p.col(1)[1] = std::min(Sig_p.col(1)[1],1 + theta_s);

        p.F = p.F_e * p.F_p;

        p.F_e = U_p * Sig_p * V_transpose_p;
        p.F_p = V_transpose_p.transpose() * Sig_p.inverse() * U_p.transpose() * p.F; //



        //update particle velocities
        Vec v_PIC(0, 0, 0);
        Vec vfLIP = p.v;
        // printf("initial flip: %f\n", p.v[1]);
        for (int i = -neighbor; i <= neighbor; i++)
        {
            for (int j = -neighbor; j <= neighbor; j++)
            {
                for (int k = -neighbor; k <= neighbor; k++)
                {

                Vector3i curr_grid = Vector3i(base_coord.x() + i, base_coord.y() + j, base_coord.z() + k);
                if (curr_grid.x() >= 0 && curr_grid.x() <= n && curr_grid.y() >= 0 && curr_grid.y() <= n && curr_grid.z() >= 0 && curr_grid.z() <= n)
                { //check bounds 0 to n in both dirs
                    Vec fx = p.x * inv_dx - curr_grid.cast<double>();
                    double N = weight(fx[0]) * weight(fx[1]) * weight(fx[2]);
                    //update PIC velocity
                    v_PIC.x() += grid[curr_grid.x()][curr_grid.y()][curr_grid.z()].x() * N;
                    v_PIC.y() += grid[curr_grid.x()][curr_grid.y()][curr_grid.z()].y() * N;
                    v_PIC.z() += grid[curr_grid.x()][curr_grid.y()][curr_grid.z()].z() * N;

                    //update FLIP velocity
                    vfLIP.x() += (grid[curr_grid.x()][curr_grid.y()][curr_grid.z()].x() - oldVelocities[curr_grid.x()][curr_grid.y()][curr_grid.z()].x()) * N;
                    vfLIP.y() += (grid[curr_grid.x()][curr_grid.y()][curr_grid.z()].y() - oldVelocities[curr_grid.x()][curr_grid.y()][curr_grid.z()].y()) * N;
                    vfLIP.z() += (grid[curr_grid.x()][curr_grid.y()][curr_grid.z()].z() - oldVelocities[curr_grid.x()][curr_grid.y()][curr_grid.z()].z()) * N;
                }
                }
            }
        }

        //update particle velocities
//        p.v = (1 - alpha) * v_PIC + alpha * vfLIP;
        p.v = v_PIC;
//        std::cout << p.v << std::endl;

        //update particle positions
        p.x += p.v * dt;
//        std::cout << p.x << std::endl;
    }
//    std::cout << "update forces\n" << std::endl;




    /**************************************/
    /************ RENNDERING **************/
    /**************************************/
    if(render){
    //     //for graphical rendering using shaders
    //     std::vector<Vector3f> vertices;
    //     //vertices of tetrahedron centered at origin with side length ~=offset
    //     float offset = 0.01;
    //     Vector3f v_1 = Vector3f(offset, offset, offset);
    //     Vector3f v_2 = Vector3f(-offset, -offset, offset);
    //     Vector3f v_3 = Vector3f(-offset, offset, -offset);
    //     Vector3f v_4 = Vector3f(offset, -offset, -offset);

    //     for(int i = 0; i < particles.size(); i++){
    //         Vector3d p = particles[i].x;
    // //        cout << p << endl;
    //         Vector3f pos = Vector3f(5*(float)p.x(), 5*(float)p.y(), 5*(float)p.z());

    //         vertices.push_back(v_1+pos);
    //         vertices.push_back(v_2+pos);
    //         vertices.push_back(v_3+pos);
    //         vertices.push_back(v_4+pos);
    //     }
    //     m_shape.setVertices(vertices);
    }

    if(write){
        write_to_CSV();
        iteration++;
    }
}


// Seed particles with position and color
void mpm_solver::add_object(Vec center, int c)
{
    // Randomly sample num_particles particles in the square
    for (int i = 0; i < num_particles; i++)
    {
        particles.push_back(Particle((Vec(rand()/(float)RAND_MAX, rand()/(float)RAND_MAX, rand()/(float)RAND_MAX) * 2.0f - Vec(1, 1, 1)) * 0.08 + center, c));
    }
}

// void mpm_solver::draw(Shader *shader)
// {
//     m_shape.draw(shader);
// }

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

void mpm_solver::add_from_csv(char *infile_path, Vec center, int c)
{
    std::ifstream infile(infile_path);

    std::vector<std::string> split_line = getNextLineAndSplitIntoTokens(infile);
    split_line = getNextLineAndSplitIntoTokens(infile); //read two lines because first line is labels
    while (split_line.size() > 1)
    {
        double x_pos = std::stof(split_line[0].c_str()) / 50.f - 0.5f;
        double y_pos = std::stof(split_line[1].c_str()) / 50.f - 0.5f;
        double z_pos = std::stof(split_line[2].c_str()) / 50.f - 0.5f;
        particles.push_back(Particle(Vec(x_pos + center[0], y_pos + center[1], z_pos + center[2]), c)); //using x and y coordinates for 2d
        split_line = getNextLineAndSplitIntoTokens(infile);
    }
}

void mpm_solver::write_to_CSV(){
    std::ofstream myFile;
    myFile.open("/Users/yuna.hiraide/Desktop/snow_maya/square" + to_string(iteration)+".csv");

    for (auto &p : particles)
    {
        cout << "writinng " << endl;
        myFile << p.x.x() << ", " << p.x.y() << ", " << p.x.z() << "\n";
    }
    myFile.close();
}

int mpm_solver::run(int argc, char* argv[])
{

//    GUI gui("double-time 2D MLS-MPM", window_size, window_size);
//    auto &canvas = gui.get_canvas();

//    if (argc == 1) //default
//    {
//        add_object(Vec(0.55, 0.45, 0.f), 0xFFFAFA);
////        add_object(Vec(0.45, 0.65, 0.f), 0xFFFAFA);
////        add_object(Vec(0.55, 0.85, 0.f), 0xFFFAFA);
//    }
//    else
//    {
//        assert(argc == 2);
//        add_from_csv(argv[1], Vec(0.55, 0.85, 0.f), 0xF2B134);
//    }

//    int frame = 0;

    //initialize particle values
    initialize();

    // Main Loop
//    for (int step = 0;; step++)
//    {
        // Advance simulation
//        update(dt);

        // Visualize frame
//        if (step % int(frame_dt / dt) == 0)
//        {
            // Clear background
//            canvas.clear(0x112F41);
            // Box
//            canvas.rect(Vec(0.04), Vec(0.96)).radius(2).color(0x4FB99F).close();
            // Particles
//            for (auto p : particles)
//            {
//                canvas.circle(p.x).radius(2).color(p.c);
//            }
            // Update image
//            gui.update();

            // Write to disk (optional)
            // canvas.img.write_as_image(fmt::format("tmp/{:05d}.png", frame++));
//        }
//    }
return 0;
}

int main(){
    std::cout<<"753"<<std::endl;
    mpm_solver ms = mpm_solver();
//         std::cout<<"755"<<std::endl;

//     ms.initialize();
//     for (int step = 0;; step++)
//    {
//        ms.update(ms.dt);
//    }
    return 0;

}