#ifndef MPM_SOLVER_H
#define MPM_SOLVER_H

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include "graphics/shape.h"

using namespace Eigen;
using namespace std;

using Vec = Vector3d;
using Mat = Matrix3d;



struct Particle
{
    // Position and velocity
    Vec x, v;
    //volume
    double vol;
    // Deformation gradient
    Mat F;
    // Affine momentum from APIC
    Mat C;
    // Determinant of the deformation gradient (i.e. volume)
    double Jp;
    Mat F_e; //elastic and plastic components of the deformation gradient
    Mat F_p;
    // Color
    int c;


    Particle(Vec x, int c, Vec v = Vec(0, 0, 0)) : x(x),
                                             v(v),
                                             vol(0),
                                             F(Matrix3d::Identity(3, 3)),
                                             C(Matrix3d::Zero(3, 3)),
                                             Jp(1),
                                             F_e(Matrix3d::Identity(3, 3)),
                                             F_p(Matrix3d::Identity(3, 3)),
                                             c(c) {}
};

class mpm_solver        
{
public:
    mpm_solver();
    void initialize();
    void update(double dt);
    int run(int argc, char* argv[]);
    void draw(Shader *shader);
    void add_from_csv(char *infile_path, Vec center, int c);


    std::vector<Particle> particles;

    void vdb_write(std::string file);

private:
    //for rendering
    Shape m_shape;

    // Window
    const int window_size = 800;

    // Grid resolution (cells)
    const static int n = 40;

    //number of particles per object
    const int num_particles = 50.f;

    const double dt = 1e-4f;
    const double frame_dt = 1e-3f;
    const double dx = 1.0f / n;
    const double inv_dx = 1.0f / dx;

    const double alpha = 0.95;

    // Snow material properties
    const double particle_mass = 1.0f;
    const double vol = 1.0f; // Particle Volume
    const double xi = 10.0f; // Snow hardening factor
    const double E = 1e4f;   // Young's Modulus
    const double nu = 0.2f;  // Poisson ratio
    const bool plastic = true;

    // Initial Lamé parameters
    const double mu_0 = E / (2.0f * (1.0f + nu));
    const double lambda_0 = E * nu / ((1.0f + nu) * (1.0f - 2.0f * nu));

    //neighbor grid
    const int neighbor = 2;
    int iteration = 0;



    // Vector3: [velocity_x, velocity_y, mass]
//    Vector3d grid[n + 1][n + 1][n + 1];
    Vector4d grid[n + 1][n + 1][n + 1];



    void add_object(Vec center, int c);
    void write_to_CSV();
};

#endif // MPM_SOLVER_H
