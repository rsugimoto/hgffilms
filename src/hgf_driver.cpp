#include <map>
#include <set>
#include <tuple>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include "hgf_driver.hpp"
#include "utils.hpp"
#include "trianglequality.h"
#include <igl/barycentric_coordinates.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>

using namespace LosTopos;

HGFDriver::HGFDriver(SurfTrack* st, int num_open_regions, double surface_tension_strength):
    NUM_OPEN_REGIONS(num_open_regions),
    st(st),
    beta(surface_tension_strength*2.0),
    mesh_update_event_handled(0)
{
    this->velocities = Eigen::MatrixXd::Zero(st->get_num_vertices(), 3);
    initial_region_volumes = get_volumes_by_region(st);
}

void HGFDriver::update_velocities_after_mesh_improvement() {
    for(; mesh_update_event_handled<st->m_mesh_change_history.size(); mesh_update_event_handled++) {
        const MeshUpdateEvent& event = st->m_mesh_change_history[mesh_update_event_handled];
        for ( const size_t & v: event.m_created_verts ){
            assert(v == velocities.rows());
            if(v >= velocities.rows()) velocities.conservativeResize(v+1, Eigen::NoChange);

            if( event.m_type == MeshUpdateEvent::EDGE_SPLIT) {
                assert(event.m_created_verts.size()==1);
                auto v0 = st->get_position(event.m_v0);
                auto v1 = st->get_position(event.m_v1);
                auto pos = event.m_vert_position;
                double w0 = mag(v0-pos);
                double w1 = mag(v1-pos);
                double w = w1/(w0+w1);
                velocities.row(event.m_created_verts[0]) = velocities.row(event.m_v0) * w + velocities.row(event.m_v1) * (1.0-w);
            }else if ( event.m_type == MeshUpdateEvent::VERTEX_POP) {
                //tentative implementation. Just taking average...
                std::set<size_t> adjacent_vertices_set;
                for (const Vec3st & tri: event.m_created_tri_data ){
                    if(tri[0] == v || tri[1] == v || tri[2] == v) {
                        if(tri[0] != v) adjacent_vertices_set.insert( tri[0] );
                        if(tri[1] != v) adjacent_vertices_set.insert( tri[1] );
                        if(tri[2] != v) adjacent_vertices_set.insert( tri[2] );
                    }
                }
                Eigen::RowVector3d res = Eigen::RowVector3d::Zero();
                int count = 0;
                for ( const size_t& adj: adjacent_vertices_set ){
                    if(adj>=velocities.rows())continue;
                    res += velocities.row(adj);
                    count++;
                }    
                velocities.row(v) = res/count;
            } else {
                std::set<size_t> adjacent_vertices_set;
                for (const Vec3st & tri: event.m_created_tri_data ){
                    if(tri[0] == v || tri[1] == v || tri[2] == v) {
                        if(tri[0] != v) adjacent_vertices_set.insert( tri[0] );
                        if(tri[1] != v) adjacent_vertices_set.insert( tri[1] );
                        if(tri[2] != v) adjacent_vertices_set.insert( tri[2] );
                    }
                }
                assert(adjacent_vertices_set.size() == 3);
                std::vector<size_t> adjacent_vertices(adjacent_vertices_set.begin(), adjacent_vertices_set.end());
                auto pos = st->get_position(v);
                auto v0 = st->get_position(adjacent_vertices[0]);
                auto v1 = st->get_position(adjacent_vertices[1]);
                auto v2 = st->get_position(adjacent_vertices[2]);
                Eigen::RowVector3d coord;
                igl::barycentric_coordinates(to_eigen_row_vector_const(pos), to_eigen_row_vector_const(v0), to_eigen_row_vector_const(v1), to_eigen_row_vector_const(v2), coord);
                velocities.row(v) = coord[0]*velocities.row(adjacent_vertices[0])
                                    + coord[1]*velocities.row(adjacent_vertices[1]) 
                                    + coord[2]*velocities.row(adjacent_vertices[2]);
            }
        }
    }
}

void HGFDriver::update_velocities_after_mesh_defrag(){
    Eigen::MatrixXd new_velocities(st->get_num_vertices(), 3);
    for (const Vec2st& pair: st->m_defragged_vertex_map) {
        const size_t &prev_index = pair[0];
        const size_t &curr_index = pair[1];
        new_velocities.row(curr_index) = velocities.row(prev_index);
    }
    velocities.swap(new_velocities);
}

void HGFDriver::step_evolution(double dt) {
    //Step 1: Evolution by dA/dx
    auto positions =  to_eigen(st->pm_positions);
    auto new_positions =  to_eigen(st->pm_newpositions);

    // std::vector<Eigen::Triplet<double>> L_elements;
    // for (size_t i=0; i<st->get_num_vertices(); ++i){
    //     double sum_element = 0.0;
    //     for ( const size_t &e : st->m_mesh.m_vertex_to_edge_map[i] ) {
    //         const Vec2st &edge = st->m_mesh.m_edges[e];
    //         size_t j = edge[0]==i? edge[1] : edge[0];
    //         double element = 0.0;
    //         for (const size_t &f : st->m_mesh.m_edge_to_triangle_map[e] ) {
    //             size_t k = st->m_mesh.get_third_vertex(e, f);
    //             const Vec3d &pi = st->get_position(i);
    //             const Vec3d &pj = st->get_position(j);
    //             const Vec3d &pk = st->get_position(k);
    //             // double R = circumcircle_radius(pi, pj, pk);
    //             Vec3d pipk = pi-pk;
    //             Vec3d pjpk = pj-pk;
    //             double cos_k = dot(pipk, pjpk)/(mag(pipk)*mag(pjpk));
    //             double sin_k2 = 1.0-cos_k*cos_k;
    //             element += cos_k/sqrt(sin_k2);//1.0/(2.0*R*R*sin_k2);
    //         }
    //         L_elements.push_back(Eigen::Triplet<double>(i, j, element));
    //         sum_element += element;
    //     }
    //     L_elements.push_back(Eigen::Triplet<double>(i, i, -sum_element));
    // }
    // Eigen::SparseMatrix<double> L(positions.rows(), positions.rows());
    // L.setFromTriplets(L_elements.begin(), L_elements.end());
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(positions, to_eigen_const(st->m_mesh.m_tris).cast<int>(), L);
    Eigen::SparseMatrix<double> M;
    igl::massmatrix(positions, to_eigen_const(st->m_mesh.m_tris).cast<int>(), igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(M);
    Eigen::MatrixXd dAdx = -solver.solve(L*positions);
    velocities -= dt*(beta*dAdx);
    new_positions = positions + dt*velocities;
}

void HGFDriver::step_volume_correction(double dt) {
    //Seems get_volumes_by_region is not working consistently.
    auto new_positions =  to_eigen(st->pm_newpositions);
    
    int num_closed_regions = 0;
    for (auto const& [region_id, volume]: initial_region_volumes){
        if (region_id >= NUM_OPEN_REGIONS) num_closed_regions++;
    }

    std::map<int, int> region_id_to_mat_index, mat_index_to_region_id;
    Eigen::VectorXd dV(num_closed_regions);
    {//calculate dV
        std::map<int, double> region_volumes = get_predicted_volumes_by_region(st);
        int i = 0;
        for (auto const& [region_id, volume]: region_volumes){
            std::cout<<"region"<<region_id<<": "<<volume<<std::endl;
            if (region_id < NUM_OPEN_REGIONS) continue;
            region_id_to_mat_index[region_id] = i;
            mat_index_to_region_id[i] = region_id;
            dV[i] = initial_region_volumes[region_id] - volume;
            i++;
        }
    }

    Eigen::SparseMatrix<double> A(num_closed_regions, num_closed_regions);
    {//calculate A
        std::map<std::tuple<int, int>, double> interface_areas = get_predicted_interface_areas(st);
        std::vector<Eigen::Triplet<double>> A_elements;
        for( auto const& [region_ids, area] : interface_areas ){
            auto [i, j] = region_ids;
            if (i >= NUM_OPEN_REGIONS && j >= NUM_OPEN_REGIONS){
                int i_index = region_id_to_mat_index[i];
                int j_index = region_id_to_mat_index[j];
                A_elements.push_back(Eigen::Triplet<double>(i_index, j_index, -area));
                A_elements.push_back(Eigen::Triplet<double>(j_index, i_index, -area));
                A_elements.push_back(Eigen::Triplet<double>(i_index, i_index, area));
                A_elements.push_back(Eigen::Triplet<double>(j_index, j_index, area));
            }
            else if (i >= NUM_OPEN_REGIONS && j < NUM_OPEN_REGIONS) {
                int i_index = region_id_to_mat_index[i];
                A_elements.push_back(Eigen::Triplet<double>(i_index, i_index, area));
            }
            else if (i < NUM_OPEN_REGIONS && j >= NUM_OPEN_REGIONS) {
                int j_index = region_id_to_mat_index[j];
                A_elements.push_back(Eigen::Triplet<double>(j_index, j_index, area));
            }  
        }
        A.setFromTriplets(A_elements.begin(), A_elements.end());
    }

    std::map<int, double> d;
    {//calculate d from A and dV
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        Eigen::VectorXd _d = solver.solve(dV);
        for ( size_t i = 0; i <_d.rows(); ++i ) d[mat_index_to_region_id[i]] = _d[i];
        for ( size_t i = 0; i < NUM_OPEN_REGIONS; ++i ) d[i] = 0.0;
    }

    Eigen::MatrixXd delta_d_normals(st->get_num_vertices(), 3);
    {//calculate delta_d_normals
        std::set<size_t> non_manifold_vertices;
        for ( size_t i = 0; i < st->get_num_vertices(); ++i ) {
            if(st->m_mesh.is_vertex_nonmanifold(i)){ //non-manifold vertex
                non_manifold_vertices.insert(i);
            }else { //manifold vertex
                auto incident_face = st->m_mesh.m_vertex_to_triangle_map[i][0];
                auto region_ids = st->m_mesh.m_triangle_labels[incident_face];
                double delta_d = d[region_ids[0]] - d[region_ids[1]];
                Vec3d normal = get_predicted_vertex_normal_max(st, i);
                delta_d_normals.row(i) = -delta_d*to_eigen_row_vector(normal);
            }
        }
        //For non-manifold vertices, take average of adjacent vertice's values
        for ( const size_t& i: non_manifold_vertices ) {
            Eigen::Vector3d sum_delta_d_normal = Eigen::Vector3d::Zero();
            std::vector<size_t> adjacent_vertices;
            st->m_mesh.get_adjacent_vertices(i, adjacent_vertices);
            size_t manifold_adjacent_vertices_count = 0;
            for (const size_t& j: adjacent_vertices) {
                if(non_manifold_vertices.find(j) == non_manifold_vertices.end()){
                    sum_delta_d_normal = delta_d_normals.row(j);
                    manifold_adjacent_vertices_count ++;
                }
            }
            delta_d_normals.row(i) = sum_delta_d_normal/manifold_adjacent_vertices_count;
        }
    }
    new_positions += delta_d_normals;
    velocities += (1./dt)*delta_d_normals;
}

void HGFDriver::step(double dt) {
    std::cout<<"HGFDriver::step"<<std::endl<<std::endl;
    step_evolution(dt);
    step_volume_correction(dt);
}

