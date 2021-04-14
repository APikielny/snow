//#include "simulation.h"

//#include <iostream>
//#include <unordered_set>
//#include <map>

//#include "graphics/MeshLoader.h"

//using namespace Eigen;

//static Vector3d computeNormal(Vector3d v1, Vector3d v2, Vector3d v3){
//    Vector3d n = (v2-v1).cross(v3-v1);
//    n.normalize();
//    return n;
//}

//static Vector3d getGroundIntersection(Vector3d pos){
//    pos = Affine3d(Translation3d(0, 2, 0))*pos;
//    double dist = 0;
//    if((abs(pos.x()) <= 5) && (abs(pos.z()) <= 5) && (pos.y() <= 0)){
//        dist = abs(0-pos.y());
//    }
//    return dist*Vector3d(0, 1, 0);
//}
//static Vector3d getSphereIntersection(Vector3d pos){
//    pos = Affine3d(Translation3d(0, 2, 0))*pos;
//    Vector3d center = Vector3d(1, 0, 0);
//    double r = 1;


//    double dist = (pos-center).norm() - r;
//    if (dist < 0){
//        Vector3d norm = pos-center;
//        norm.normalize();
//        return abs(dist)*norm;
//    } else{
//        return Vector3d(0, 0, 0);
//    }

//}



//Simulation::Simulation()
//{
//}

//void Simulation::init()
//{
//    m_solver.initialize();

//    std::vector<Vector3f> vertices;
//    std::vector<Vector3i> faces;
//    int j = 0;

//    //vertices of tetrahedron centered at origin with side length ~=offset
//    float offset = 0.1;
//    Vector3f v_1 = Vector3f(offset, offset, offset);
//    Vector3f v_2 = Vector3f(-offset, -offset, offset);
//    Vector3f v_3 = Vector3f(-offset, offset, -offset);
//    Vector3f v_4 = Vector3f(offset, -offset, -offset);

//    for(int i = 0; i < 1000; i++){
//        Vector2d p = m_solver.particles[i].x;
//        Vector3f pos = Vector3f(5*(float)p.x(), 5*(float)p.y(), 0);

//        vertices.push_back(v_1+pos);
//        vertices.push_back(v_2+pos);
//        vertices.push_back(v_3+pos);
//        vertices.push_back(v_4+pos);

//        Vector4i t = Vector4i(j, j+1, j+2, j+3);
//        Vector3i f_1 = Vector3i(t[3], t[1], t[2]); //opposite t[0]
//        Vector3i f_2 = Vector3i(t[2], t[0], t[3]); //opposite t[1]
//        Vector3i f_3 = Vector3i(t[3], t[0], t[1]); //opposite t[2]
//        Vector3i f_4 = Vector3i(t[1], t[0], t[2]); //opposite t[3]

//        faces.emplace_back(f_1);
//        faces.emplace_back(f_2);
//        faces.emplace_back(f_3);
//        faces.emplace_back(f_4);

//        j+=4;
//    }
//    m_shape.init(vertices, faces);
//}

//void Simulation::update(float seconds)
//{

////    //clear forces - ie set to just gravity
////    for(int i = 0; i < m_particles.size(); i++){
////        std::shared_ptr<particle> p = m_particles[i];
////        p->force = GRAVITY * p->mass;
////    }

////    //compute stess, strain and forces for each node
////    for(int i = 0; i < m_tets.size(); i++){
////        std::shared_ptr<tetrahedron> t = m_tets[i];

////        std::shared_ptr<particle> p_1 = m_particles[t->verts[0]];
////        std::shared_ptr<particle> p_2 = m_particles[t->verts[1]];
////        std::shared_ptr<particle> p_3 = m_particles[t->verts[2]];
////        std::shared_ptr<particle> p_4 = m_particles[t->verts[3]];

////        std::vector<Vector3d> posisions = {p_1->pos, p_2->pos, p_3->pos, p_4->pos};
////        std::vector<Vector3d> velocities = {p_1->vel, p_2->vel, p_3->vel, p_4->vel};
////        std::vector<Vector3d> forces;
////        forces.reserve(4);

////        computeStressStrainForces(posisions, velocities, t, forces);
////        p_1->force += forces[0];
////        p_2->force += forces[1];
////        p_3->force += forces[2];
////        p_4->force += forces[3];
////    }


//////    //do midpoint euler step
////    for(int i = 0; i < m_particles.size(); i++){
////        std::shared_ptr<particle> p = m_particles[i];
////        p->midpointPos = p->pos + p->vel*(0.5f*seconds);
////        Vector3d intersection = getGroundIntersection(p->midpointPos);
////        p->force += intersection*PENALTY_CONSTANT;
////        if(SPHERE_INTERSECT){
////            intersection = getSphereIntersection(p->midpointPos);
////            p->force += intersection*PENALTY_CONSTANT;
////        }
////        Vector3d a = p->force / (double) p->mass;
////        p->midpointVel = p->vel + a*(0.5f*seconds);
////    }

////    //recompute forces using midpoint positions and velocities
////    for(int i = 0; i < m_particles.size(); i++){
////        std::shared_ptr<particle> p = m_particles[i];
////        p->force = GRAVITY * p->mass;
////    }

////    //compute stess, strain and forces for each node
////    for(int i = 0; i < m_tets.size(); i++){
////        std::shared_ptr<tetrahedron> t = m_tets[i];

////        std::shared_ptr<particle> p_1 = m_particles[t->verts[0]];
////        std::shared_ptr<particle> p_2 = m_particles[t->verts[1]];
////        std::shared_ptr<particle> p_3 = m_particles[t->verts[2]];
////        std::shared_ptr<particle> p_4 = m_particles[t->verts[3]];

////        std::vector<Vector3d> posisions = {p_1->midpointPos, p_2->midpointPos, p_3->midpointPos, p_4->midpointPos};
////        std::vector<Vector3d> velocities = {p_1->midpointVel, p_2->midpointVel, p_3->midpointVel, p_4->midpointVel};
////        std::vector<Vector3d> forces;
////        forces.reserve(4);

////        computeStressStrainForces(posisions, velocities, t, forces);
////        p_1->force += forces[0];
////        p_2->force += forces[1];
////        p_3->force += forces[2];
////        p_4->force += forces[3];
////    }

////    //update position and velocity
////    for(int i = 0; i < m_particles.size(); i++){
////        std::shared_ptr<particle> p = m_particles[i];
////        p->pos = p->pos + p->vel*(seconds);
////        Vector3d intersection = getGroundIntersection(p->pos);
////        p->force += intersection*PENALTY_CONSTANT;
////        if(SPHERE_INTERSECT){
////            intersection = getSphereIntersection(p->pos);
////            p->force += intersection*PENALTY_CONSTANT;
////        }
////        Vector3d a = p->force / (double) p->mass;
////        p->vel = p->vel + a*(seconds);
////    }

////    std::vector<Vector3f> vertices;
////    for(int i = 0; i < m_particles.size(); i++){
////        Vector3d pPos = m_particles[i]->pos;
////        vertices.push_back(Vector3f(pPos.x(), pPos.y(), pPos.z()));
////    }

////    m_shape.setVertices(vertices);
//}

//void Simulation::computeStressStrainForces(const std::vector<Vector3d> &positions, const std::vector<Vector3d> &velocities, const std::shared_ptr<tetrahedron> t, std::vector<Vector3d> &outForces){
//    Vector3d p1 = positions[0];
//    Vector3d p2 = positions[1];
//    Vector3d p3 = positions[2];
//    Vector3d p4 = positions[3];

//    Vector3d v1 = velocities[0];
//    Vector3d v2 = velocities[1];
//    Vector3d v3 = velocities[2];
//    Vector3d v4 = velocities[3];

//    //compute P matrix
//    Matrix3d P;
//    P.col(0) = p1 - p4;
//    P.col(1) = p2 - p4;
//    P.col(2) = p3 - p4;

//    //compute V matrix
//    Matrix3d V;
//    V.col(0) = v1 - v4;
//    V.col(1) = v2 - v4;
//    V.col(2) = v3 - v4;

//    //compute derivatives
//    Matrix3d dx_du = P*t->beta;
//    Matrix3d dv_du = V*t->beta;

//    //compute strain + stress due to strain
//    Matrix3d epsilon = dx_du.transpose()*dx_du - Matrix3d::Identity();
//    Matrix3d sigma_epsilon = m_lambda*epsilon.trace()*Matrix3d::Identity() + 2.0f*m_mu*epsilon;

//    //compute strain rate + stress due to strain rate
//    Matrix3d epsilon_dot = dx_du.transpose()*dv_du + dv_du.transpose()*dx_du;
//    Matrix3d sigma_epsilon_dot = m_phi*epsilon_dot.trace()*Matrix3d::Identity() + 2.0f*m_psi*epsilon_dot;

//    //compute stress
//    Matrix3d sigma = sigma_epsilon+sigma_epsilon_dot;

//    double f_1 = t->areas[0]; //area of face opposite p_1
//    Vector3d n_1 = t->normals[0];
//    double f_2 = t->areas[1]; //area of face opposite p_2
//    Vector3d n_2 = t->normals[1];
//    double f_3 = t->areas[2]; //area of face opposite p_3
//    Vector3d n_3 = t->normals[2];
//    double f_4 = t->areas[3]; //area of face opposite p_4
//    Vector3d n_4 = t->normals[3];

//    outForces[0] = dx_du*sigma*(f_1*n_1);
//    outForces[1] = dx_du*sigma*(f_2*n_2);
//    outForces[2] = dx_du*sigma*(f_3*n_3);
//    outForces[3] = dx_du*sigma*(f_4*n_4);
//}


//void Simulation::draw(Shader *shader)
//{
//    m_shape.draw(shader);
//    m_ground.draw(shader);
//    if(SPHERE_INTERSECT){
//        m_sphere.draw(shader);
//    }
//}

//void Simulation::toggleWire()
//{
//    m_shape.toggleWireframe();
//}

//void Simulation::initGround()
//{
//    std::vector<Vector3f> groundVerts;
//    std::vector<Vector3i> groundFaces;
//    groundVerts.emplace_back(-5, 0, -5);
//    groundVerts.emplace_back(-5, 0, 5);
//    groundVerts.emplace_back(5, 0, 5);
//    groundVerts.emplace_back(5, 0, -5);
//    groundFaces.emplace_back(0, 1, 2);
//    groundFaces.emplace_back(0, 2, 3);
//    m_ground.init(groundVerts, groundFaces);
//}

//void Simulation::initSphere()
//{
//    std::vector<Vector3f> sphereVerts;
//    std::vector<Vector4i> sphereTets;
//    std::vector<Vector3i> faceSet;
//    std::vector<Vector3f> newVerts;
//    double offset = 0.1;
//    if(MeshLoader::loadTetMesh("example-meshes/sphere.mesh", sphereVerts, sphereTets)) {
//        int j = 0;
//        for(int i = 0; i < sphereVerts.size(); i++){
//            Vector3f v_1 = Vector3f(offset, offset, offset);
//            Vector3f v_2 = Vector3f(-offset, -offset, offset);
//            Vector3f v_3 = Vector3f(-offset, offset, -offset);
//            Vector3f v_4 = Vector3f(offset, -offset, -offset);
//            newVerts.push_back(v_1+sphereVerts[i]);
//            newVerts.push_back(v_2+sphereVerts[i]);
//            newVerts.push_back(v_3+sphereVerts[i]);
//            newVerts.push_back(v_4+sphereVerts[i]);

//            Vector4i t = Vector4i(j, j+1, j+2, j+3);
//            Vector3i f_1 = Vector3i(t[3], t[1], t[2]); //opposite t[0]
//            Vector3i f_2 = Vector3i(t[2], t[0], t[3]); //opposite t[1]
//            Vector3i f_3 = Vector3i(t[3], t[0], t[1]); //opposite t[2]
//            Vector3i f_4 = Vector3i(t[1], t[0], t[2]); //opposite t[3]

//            faceSet.emplace_back(f_1);
//            faceSet.emplace_back(f_2);
//            faceSet.emplace_back(f_3);
//            faceSet.emplace_back(f_4);
//            j+=4;

//        }
//    }
//    m_sphere.init(newVerts, faceSet);
//    m_sphere.setModelMatrix(Affine3f(Eigen::Translation3f(1, 0, 0)));


//}
