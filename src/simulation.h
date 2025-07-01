#pragma once

#include "graphics/shape.h"

class Shader;

class Simulation
{
public:
    Simulation();

    void init();

    void update(double seconds);

    void draw(Shader *shader);

    void toggleWire();
private:
    const double m_density = 5000.0;
    const double m_lambda = .08;
    const double m_mu = 0.6;
    const double m_phi = .1;
    const double m_psi = .15;

    int m_frame = 0;

    struct Tet;

    struct Vertex{
        Eigen::Vector3d *position;
        Eigen::Vector3d velocity = Eigen::Vector3d(0.0, 0.0, 0.0);
        Eigen::Vector3d force = Eigen::Vector3d(0.0,0.0,0.0);
        //std::vector<Tet*> tets;
        double mass;
        std::unordered_map<Tet*, Eigen::Vector3d> area_weighted_normals;
    };

    struct Tet{
        double mass;
        double volume;
        std::vector<Vertex*> vertices;
        Eigen::Matrix<double, 4, 4> beta_mat = Eigen::Matrix<double, 4, 4>::Zero();
        Eigen::Matrix<double, 3, 4> P = Eigen::Matrix<double, 3, 4>::Zero();
        Eigen::Matrix<double, 3, 4> V = Eigen::Matrix<double, 3, 4>::Zero();
        Eigen::Matrix<double, 3, 3> deformation_gradient = Eigen::Matrix<double, 3, 3>::Zero();
        Eigen::Matrix<double, 3, 3> deformation_rate_gradient = Eigen::Matrix<double, 3, 3>::Zero();
        Eigen::Matrix<double, 3, 3> strain_matrix = Eigen::Matrix<double, 3, 3>::Zero();
        Eigen::Matrix<double, 3, 3> strain_rate = Eigen::Matrix<double, 3, 3>::Zero();
        Eigen::Matrix3d stress_matrix = Eigen::Matrix<double, 3, 3>::Zero();
    };

    Shape m_shape;

    Shape m_ground;
    std::vector<Eigen::Vector3d> m_vertices_material;
    std::vector<Eigen::Vector3d> m_vertices_world;
    std::vector<Vertex*> m_vertices;
    std::vector<Tet*> m_tets;
    std::vector<Eigen::Vector3i> m_faces;
    Eigen::Vector3d m_gravity = Eigen::Vector3d(0.0, -0.001, 0.0);
    void initGround();
    void put_faces_in_map(std::unordered_map<std::string, std::tuple<int, Eigen::Vector3i>>* map, std::vector<int> vertices);
    void find_tet_vol(Tet* tet);
    void apply_velocity();
    void update_velocity();
    void accumForces();
    void applyGravity();
    void check_collis();
    void compute_strain();
    void compute_stress(Tet* tet);
    void add_stress(Tet* tet);
    void stress_per_v(Tet* tet, Vertex* v, Eigen::Vector3d norm1, Eigen::Vector3d norm2, Eigen::Vector3d norm3);
    void find_normals(Tet* tet);
    void midpoint_method();
    void export_obj();
};
