//#ifndef SIMULATION_H
//#define SIMULATION_H

//#include "graphics/shape.h"
//#include <unordered_set>
//#include <string>
//#include "mpm_solver.h"


//using namespace Eigen;

//class Shader;

//struct particle{
//    Vector3d materialSpacePos;
//    Vector3d pos;
//    Vector3d midpointPos;
//    Vector3d vel;
//    Vector3d midpointVel;
//    Vector3d force;
//    double mass;
//};

//struct tetrahedron{
//    Vector4i verts;
//    std::vector<Vector3d> normals; //index here maps to face opposite vert
//    std::vector<double> areas; // index here maps to face opposite vert
//    Matrix3d beta;
//};

//class Simulation
//{
//public:
//    Simulation();

//    void init();

//    void update(float seconds);

//    void draw(Shader *shader);

//    void toggleWire();
//private:
//    mpm_solver m_solver;

//    Shape m_shape;

//    Shape m_ground;
//    Shape m_sphere;

//    std::vector<std::shared_ptr<particle>> m_particles;
//    std::vector<std::shared_ptr<tetrahedron>> m_tets;
//    const double m_density = 1200.0f;
//    const double m_lambda = 1.0*pow(10, 3);
//    const double m_mu = 1.0*pow(10, 3);
//    const double m_phi = 2.0*pow(10, 1);
//    const double m_psi = 2.0*pow(10, 1);
//    const bool SPHERE_INTERSECT = true;
//    const float PENALTY_CONSTANT = 10000.f;
//    const Vector3d GRAVITY = Vector3d(0, -.1f, 0);
//    const std::string MESH_FILE = "example-meshes/sphere.mesh";

//    void initGround();
//    void initSphere();
//    void computeStressStrainForces(const std::vector<Vector3d> &positions, const std::vector<Vector3d> &velocities, const std::shared_ptr<tetrahedron> t, std::vector<Vector3d> &outForces);
//};

//#endif // SIMULATION_H
