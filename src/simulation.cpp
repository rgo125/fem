#include "simulation.h"
#include "graphics/meshloader.h"

#include <iostream>
#include <fstream>

using namespace Eigen;

Simulation::Simulation() {}

void Simulation::init()
{
    // STUDENTS: This code loads up the tetrahedral mesh in 'example-meshes/single-tet.mesh'
    //    (note: your working directory must be set to the root directory of the starter code
    //    repo for this file to load correctly). You'll probably want to instead have this code
    //    load up a tet mesh based on e.g. a file path specified with a command line argument.
    std::vector<Vector3d> vertices;
    std::vector<Vector4i> tets;
    if (MeshLoader::loadTetMesh("/Users/rainergardner-olesen/Desktop/Spring_2025_Courses/Graphics/projects/fem-rgo125/example-meshes/ellipsoid.mesh", vertices, tets)) {
        // STUDENTS: This code computes the surface mesh of the loaded tet mesh, i.e. the faces
        //    of tetrahedra which are on the exterior surface of the object. Right now, this is
        //    hard-coded for the single-tet mesh. You'll need to implement surface mesh extraction
        //    for arbitrary tet meshes. Think about how you can identify which tetrahedron faces
        //    are surface faces...
        std::vector<Vector3i> faces;
        std::unordered_map<std::string, std::tuple<int, Vector3i>> check_dups;
        for(int i = 0; i < tets.size(); i++){

            std::vector<int> face1 = {tets[i].y(), tets[i].x(), tets[i].z()};
            put_faces_in_map(&check_dups, face1);

            std::vector<int> face2 = {tets[i].z(), tets[i].x(), tets[i].w()};
            put_faces_in_map(&check_dups, face2);

            std::vector<int> face3 = {tets[i].w(), tets[i].y(), tets[i].z()};
            put_faces_in_map(&check_dups, face3);

            std::vector<int> face4 = {tets[i].w(), tets[i].x(), tets[i].y()};
            put_faces_in_map(&check_dups, face4);

        }
        for(const auto& pair : check_dups){
            if(std::get<0>(pair.second) == 1){
                faces.push_back(std::get<1>(pair.second));
            }
        }
        m_faces = faces;
        m_shape.init(vertices, faces, tets);
    }

    //set material space vertices
    m_vertices_material = vertices;
    //set world space as well but this will change
    for(int i = 0; i < vertices.size(); i++){
        m_vertices_world.push_back(vertices[i]);
    }


    for(int i = 0; i < m_vertices_world.size(); i++){
        Vertex* v = new Vertex();
        v->position = &(m_vertices_world[i]);
        m_vertices.push_back(v);
    }

    for(int i = 0; i < tets.size(); i++){
        Tet* tet = new Tet();
        Matrix<double, 4, 4> mat = Matrix<double, 4, 4>::Zero();
        for(int j = 0; j < 4; j++){
            Vertex* vert = m_vertices[tets[i][j]];

            tet->vertices.push_back(vert);

            mat(0,j) = m_vertices_material[tets[i][j]].x();

            mat(1,j) = m_vertices_material[tets[i][j]].y();

            mat(2,j) = m_vertices_material[tets[i][j]].z();
        }
        mat.row(3).setOnes();
        tet->beta_mat = mat.inverse();

        find_tet_vol(tet);
        find_normals(tet);

        m_tets.push_back(tet);
    }


    m_shape.setModelMatrix(Affine3f(Eigen::Translation3f(0, 2, 0)));


    initGround();
}

void Simulation::find_normals(Tet* tet){
    Vertex* v1 = tet->vertices[0];
    Vertex* v2 = tet->vertices[1];
    Vertex* v3 = tet->vertices[2];
    Vertex* v4 = tet->vertices[3];

    Vector3d face_1_norm = Vector3d(*v1->position - *v2->position).cross(Vector3d(*v3->position - *v2->position));
    Vector3d face_2_norm = Vector3d(*v1->position - *v3->position).cross(Vector3d(*v4->position - *v3->position));
    Vector3d face_3_norm = Vector3d(*v2->position - *v4->position).cross(Vector3d(*v3->position - *v4->position));
    Vector3d face_4_norm = Vector3d(*v1->position - *v4->position).cross(Vector3d(*v2->position - *v4->position));

    v1->area_weighted_normals[tet] = (face_1_norm + face_2_norm + face_4_norm).normalized();
    v2->area_weighted_normals[tet] = (face_1_norm + face_3_norm + face_4_norm).normalized();
    v3->area_weighted_normals[tet] = (face_1_norm + face_2_norm + face_3_norm).normalized();
    v4->area_weighted_normals[tet] = (face_2_norm + face_3_norm + face_4_norm).normalized();

}

void Simulation::put_faces_in_map(std::unordered_map<std::string, std::tuple<int, Vector3i>>* map, std::vector<int> vertices){
    Vector3i face = {vertices[0], vertices[1], vertices[2]};
    std::sort(vertices.begin(), vertices.end());
    std::string key = "";
    for(int vert: vertices){
        key = key + std::to_string(vert) + ", ";
    }
    if(map->find(key) != map->end()){
        std::get<0>(map->at(key))++;
    }
    else{
        map->insert({key, std::tuple(1,face)});
    }
}

void Simulation::find_tet_vol(Tet* tet){
    Vertex* vert1 = tet->vertices[0];
    Vertex* vert2 = tet->vertices[1];
    Vertex* vert3 = tet->vertices[2];
    Vertex* vert4 = tet->vertices[3];

    Vector3d A = *vert1->position - *vert2->position;
    Vector3d B = *vert3->position - *vert2->position;
    Vector3d C = *vert4->position - *vert2->position;

    tet->volume = std::abs((1.0/6.0) * (A.cross(C)).dot(B));
    tet->mass = m_density * tet->volume;
    double v_mass = tet->mass/4.0;
    vert1->mass = v_mass;
    vert2->mass = v_mass;
    vert3->mass = v_mass;
    vert4->mass = v_mass;
}

void Simulation::update(double seconds)
{
    // STUDENTS: This method should contain all the time-stepping logic for your simulation.
    //   Specifically, the code you write here should compute new, updated vertex positions for your
    //   simulation mesh, and it should then call m_shape.setVertices to update the display with those
    //   newly-updated vertices.

    // STUDENTS: As currently written, the program will just continually compute simulation timesteps as long
    //    as the program is running (see View::tick in view.cpp) . You might want to e.g. add a hotkey for pausing
    //    the simulation, and perhaps start the simulation out in a paused state.

    // Note that the "seconds" parameter represents the amount of time that has passed since
    // the last update
}

void Simulation::draw(Shader *shader)
{
    m_shape.draw(shader);
    m_ground.draw(shader);
    export_obj();

    midpoint_method();

    m_shape.setVertices(m_vertices_world);
}

void Simulation::toggleWire()
{
    m_shape.toggleWireframe();
}

void Simulation::initGround()
{
    std::vector<Vector3d> groundVerts;
    std::vector<Vector3i> groundFaces;
    groundVerts.emplace_back(-5, 0, -5);
    groundVerts.emplace_back(-5, 0, 5);
    groundVerts.emplace_back(5, 0, 5);
    groundVerts.emplace_back(5, 0, -5);
    groundFaces.emplace_back(0, 1, 2);
    groundFaces.emplace_back(0, 2, 3);
    m_ground.init(groundVerts, groundFaces);
}


void Simulation::accumForces(){
    for(Vertex* v: m_vertices){
        v->force =Vector3d(0.0,0.0,0.0);
    }
}

void Simulation::update_velocity(){
    for(Vertex* v: m_vertices){
        Vector3d acceleration = v->force/v->mass;
        acceleration += m_gravity;
        v->velocity += acceleration;
    }
}

void Simulation::apply_velocity(){
    for(Vertex* v: m_vertices){
        *v->position += v->velocity;
    }
}

void Simulation::check_collis(){
    for(Vertex* vert: m_vertices){
        if(vert->position->y() < -2.0){
            vert->position->y() = -2.0;
            Vector3d coll_norm = {0.0, 1.0, 0.0};
            Vector3d norm_comp = {0.0, vert->velocity.y(), 0.0};
            Vector3d horiz_comp = {vert->velocity.x(), 0.0, vert->velocity.z()};

            vert->velocity = (-1.0 * norm_comp) + (0.8 * horiz_comp);
        }
    }
}

void Simulation::compute_strain(){
    for(Tet* tet: m_tets){
        tet->P.col(0) = *(tet->vertices[0]->position);
        tet->P.col(1) = *(tet->vertices[1]->position);
        tet->P.col(2) = *(tet->vertices[2]->position);
        tet->P.col(3) = *(tet->vertices[3]->position);

        tet->V.col(0) = (tet->vertices[0]->velocity);
        tet->V.col(1) = (tet->vertices[1]->velocity);
        tet->V.col(2) = (tet->vertices[2]->velocity);
        tet->V.col(3) = (tet->vertices[3]->velocity);

        tet->deformation_gradient = (tet->P * (tet->beta_mat).leftCols(3));
        tet->deformation_rate_gradient = (tet->V * (tet->beta_mat).leftCols(3));
        tet->strain_matrix = (((tet->deformation_gradient).transpose() * tet->deformation_gradient) - Matrix3d::Identity()) * (1.0/2.0);
        tet->strain_rate = (tet->deformation_gradient.transpose() * tet->deformation_rate_gradient) +
                           (tet->deformation_rate_gradient.transpose() * tet->deformation_gradient);

        compute_stress(tet);

    }
}

void Simulation::compute_stress(Tet* tet){
    tet->stress_matrix = (m_lambda * tet->strain_matrix.trace() * Matrix3d::Identity()) + (2.0 * m_mu * tet->strain_matrix) +
                         (m_phi * tet->strain_rate.trace() * Matrix3d::Identity()) + (2.0 * m_mu * tet->strain_matrix);

    if(!tet->stress_matrix.isZero()){
        add_stress(tet);
    }
}

void Simulation::add_stress(Tet* tet){
    for(Vertex* v : tet->vertices){
        Vector3d internal_force = (-1.0/3.0) * tet->deformation_gradient * tet->stress_matrix * (v->area_weighted_normals[tet]);
        v->force += (internal_force);
    }
}

void Simulation::midpoint_method(){
    accumForces();

    check_collis();
    compute_strain();
    update_velocity();
    apply_velocity();

    for(Vertex* v: m_vertices){
        //find midpoint
        *(v->position) = *(v->position) - (v->velocity/2.0);
    }
    accumForces();
    check_collis();
    compute_strain();

    for(Vertex* v: m_vertices){
        *(v->position) = *(v->position) - (v->velocity/2.0);
    }
    //check_collis();
    update_velocity();
    apply_velocity();

}

void Simulation::export_obj(){
    std::string filename = "/Users/rainergardner-olesen/Desktop/Spring_2025_Courses/Graphics/projects/fem-rgo125/objs/sphere/export" + std::to_string(m_frame) + ".obj";
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    // Write vertices
    for (Vertex* v: m_vertices) {
        file << "v " << v->position->x() << " " << v->position->y() << " " << v->position->z() << "\n";
    }

    // // Write texture coordinates
    // for (glm::vec2 texCoord : texture_coords) {
    //     file << "vt " << texCoord.x << " " << texCoord.y << "\n";
    // }

    // Write normals
    // for (glm::vec3 normal : normals) {
    //     file << "vn " << normal.x << " " << normal.y << " " << normal.z << "\n";
    // }

    // Write faces
    std::string curr_mtl;
    std::cout << "faces size: " << m_faces.size() << std::endl;
    for (Vector3i face : m_faces) {
        // if(curr_mtl != face.mtl){
        //     curr_mtl = face.mtl;
        //     file << "usemtl " << curr_mtl << "\n";
        // }
        file << "f";
        for (int i = 0; i < 3; i++) {
            file << " " << face[i] + 1;
            // file << "/";
            // file << face.all_refs[i][1] + 1;
            // file << "/" << face.all_refs[i][2] + 1;


        }
        file << "\n";

    }
    file.close();
    std::cout << "OBJ file written to " << filename << "\n";
    m_frame++;
}

