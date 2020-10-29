#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <vector>
#include <map>
#include <Eigen/Core>
#include "surftrack.h"

inline LosTopos::Vec3d get_predicted_vertex_normal_max( const LosTopos::SurfTrack* const st,  size_t vertex_index );
inline std::map<int, double> get_volumes_by_region(const LosTopos::SurfTrack* const st);
inline std::map<int, double> get_predicted_volumes_by_region(const LosTopos::SurfTrack* const st);
inline std::set<std::tuple<int, int>> get_interfaces(const LosTopos::SurfTrack* const st);
inline std::map<int, double> get_surface_areas_by_region(const LosTopos::SurfTrack* const st);
inline std::map<int, double> get_predicted_surface_areas_by_region(const LosTopos::SurfTrack* const st);
inline std::map<std::tuple<int, int>, double> get_interface_areas(const LosTopos::SurfTrack* const st);
inline std::map<std::tuple<int, int>, double> get_predicted_interface_areas(const  LosTopos::SurfTrack* const st);

template<unsigned int N, typename T>
inline Eigen::Map<Eigen::Matrix<T, N, 1>> to_eigen(LosTopos::Vec<N, T>& vec){
    return Eigen::Map<Eigen::Matrix<T, N, 1>>(vec.v);
}

template<unsigned int N, typename T>
inline Eigen::Map<const Eigen::Matrix<T, N, 1>> to_eigen_const(LosTopos::Vec<N, T>& vec){
    return Eigen::Map<const Eigen::Matrix<T, N, 1>>(vec.v);
}

template<unsigned int N, typename T>
inline Eigen::Map<Eigen::Matrix<T, 1, N>> to_eigen_row_vector(LosTopos::Vec<N, T>& vec){
    return Eigen::Map<Eigen::Matrix<T, 1, N>>(vec.v);
}

template<unsigned int N, typename T>
inline Eigen::Map<const Eigen::Matrix<T, 1, N>> to_eigen_row_vector_const(LosTopos::Vec<N, T>& vec){
    return Eigen::Map<const Eigen::Matrix<T, 1, N>>(vec.v);
}

template<unsigned int N, typename T>
inline Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> to_eigen(std::vector<LosTopos::Vec<N, T>>& mat){
    return Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(mat[0].v, mat.size(), N);
}

template<unsigned int N, typename T>
inline Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> to_eigen_const(std::vector<LosTopos::Vec<N, T>>& mat){
    return Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(mat[0].v, mat.size(), N);
}

inline LosTopos::Vec3d get_predicted_vertex_normal_max( const LosTopos::SurfTrack* const st,  size_t vertex_index ) {
    const std::vector<size_t>& inc_tris = st->m_mesh.m_vertex_to_triangle_map[vertex_index];
    
    LosTopos::Vec3d sum_cross_products(0,0,0);
    
    for ( size_t i = 0; i < inc_tris.size(); ++i )
    {
        const LosTopos::Vec3st& curr_tri = st->m_mesh.get_triangle( inc_tris[i] );
        
        if ( curr_tri[0] == curr_tri[1] ) { continue; }
        
        LosTopos::Vec2st other_two;
        LosTopos::NonDestructiveTriMesh::index_in_triangle( curr_tri, vertex_index, other_two );
        
        size_t verti = (size_t)curr_tri[(int)other_two[0]];
        size_t vertnext = (size_t)curr_tri[(int)other_two[1]];
        
        LosTopos::Vec3d vi = st->get_newposition(verti) - st->get_newposition(vertex_index);
        LosTopos::Vec3d vnext = st->get_newposition(vertnext) - st->get_newposition(vertex_index);
        
        sum_cross_products += cross( vi, vnext ) / ( mag2(vi)*mag2(vnext) );
    }
    
    sum_cross_products /= mag( sum_cross_products );
    
    return sum_cross_products;
}

inline std::map<int, double> get_volumes_by_region(const LosTopos::SurfTrack* const st) {
    static const double inv_six = 1.0/6.0;
    std::map<int, double> volumes;

    const std::vector<LosTopos::Vec3st>& tris = st->m_mesh.get_triangles();
    const std::vector<LosTopos::Vec2i>& labels = st->m_mesh.get_triangle_labels();
    for(size_t t=0; t < tris.size(); ++t )
    {
        if ( tris[t][0] == tris[t][1] ) { continue; }
        const LosTopos::Vec3st& tri = tris[t];
        const LosTopos::Vec2i& lbl = labels[t];
        const int lbl0 = lbl[0];
        const int lbl1 = lbl[1];
        
        if(volumes.find(lbl0)==volumes.end()) volumes[lbl0]=0.0;
        if(volumes.find(lbl1)==volumes.end()) volumes[lbl1]=0.0;

        double local_volume = inv_six * triple(st->get_position(tri[0]), st->get_position(tri[1]), st->get_position(tri[2]));
        volumes[lbl0] += local_volume;
        volumes[lbl1] -= local_volume;
    }
    for(auto& [key, value] : volumes) {
        value = abs(value);
    }
    return volumes;
}

inline std::map<int, double> get_predicted_volumes_by_region(const LosTopos::SurfTrack* const st) {
    static const double inv_six = 1.0/6.0;
    std::map<int, double> volumes;

    const std::vector<LosTopos::Vec3st>& tris = st->m_mesh.get_triangles();
    const std::vector<LosTopos::Vec2i>& labels = st->m_mesh.get_triangle_labels();
    for(size_t t=0; t < tris.size(); ++t )
    {
        if ( tris[t][0] == tris[t][1] ) { continue; }
        const LosTopos::Vec3st& tri = tris[t];
        const LosTopos::Vec2i& lbl = labels[t];
        const int lbl0 = lbl[0];
        const int lbl1 = lbl[1];
        
        if(volumes.find(lbl0)==volumes.end()) volumes[lbl0]=0.0;
        if(volumes.find(lbl1)==volumes.end()) volumes[lbl1]=0.0;
        double local_volume = inv_six * triple(st->get_newposition(tri[0]), st->get_newposition(tri[1]), st->get_newposition(tri[2]));
        volumes[lbl0] += local_volume;
        volumes[lbl1] -= local_volume;
    }
    for(auto& [key, value] : volumes) {
        value = abs(value);
    }
    return volumes;
}

inline std::set<std::tuple<int, int>> get_interfaces(const LosTopos::SurfTrack* const st) {
    std::set<std::tuple<int, int>> interfaces;
    const std::vector<LosTopos::Vec3st>& tris = st->m_mesh.get_triangles();
    const std::vector<LosTopos::Vec2i>& labels = st->m_mesh.get_triangle_labels();

    for(size_t t=0; t < tris.size(); ++t )
    {
        if ( tris[t][0] == tris[t][1] ) { continue; }
        const LosTopos::Vec2i& lbl = labels[t]; 
        
        if(labels[t][0]<labels[t][1])interfaces.insert({labels[t][0], labels[t][1]});
        else interfaces.insert({labels[t][1], labels[t][0]});
    }
    return interfaces;
}

inline std::map<int, double> get_surface_areas_by_region(const LosTopos::SurfTrack* const st) {
    std::map<int, double> areas;

    const std::vector<LosTopos::Vec3st>& tris = st->m_mesh.get_triangles();
    const std::vector<LosTopos::Vec2i>& labels = st->m_mesh.get_triangle_labels();

    for(size_t t=0; t < tris.size(); ++t )
    {
        if ( tris[t][0] == tris[t][1] ) { continue; }
        const LosTopos::Vec3st& tri = tris[t];
        const LosTopos::Vec2i& lbl = labels[t]; 
        
        if(areas.find(lbl[0])==areas.end()) areas[lbl[0]]=0.0;
        if(areas.find(lbl[1])==areas.end()) areas[lbl[1]]=0.0;

        double local_area = st->get_triangle_area(t);
        areas[lbl[0]] += local_area;
        areas[lbl[1]] += local_area;
    }
    return areas;
}

inline std::map<int, double> get_predicted_surface_areas_by_region(const LosTopos::SurfTrack* const st) {
    std::map<int, double> areas;

    const std::vector<LosTopos::Vec3st>& tris = st->m_mesh.get_triangles();
    const std::vector<LosTopos::Vec2i>& labels = st->m_mesh.get_triangle_labels();

    for(size_t t=0; t < tris.size(); ++t )
    {
        if ( tris[t][0] == tris[t][1] ) { continue; }
        const LosTopos::Vec3st& tri = tris[t];
        const LosTopos::Vec2i& lbl = labels[t]; 
        
        if(areas.find(lbl[0])==areas.end()) areas[lbl[0]]=0.0;
        if(areas.find(lbl[1])==areas.end()) areas[lbl[1]]=0.0;

        const LosTopos::Vec3d &p0 = st->get_newposition(tris[t][0]);
        const LosTopos::Vec3d &p1 = st->get_newposition(tris[t][1]);
        const LosTopos::Vec3d &p2 = st->get_newposition(tris[t][2]);      
        double local_area = 0.5 * mag(cross(p1-p0, p2-p0));
        areas[lbl[0]] += local_area;
        areas[lbl[1]] += local_area;
    }
    return areas;
}

inline std::map<std::tuple<int, int>, double> get_interface_areas(const LosTopos::SurfTrack* const st) {
    std::map<std::tuple<int, int>, double> areas;

    const std::vector<LosTopos::Vec3st>& tris = st->m_mesh.get_triangles();
    const std::vector<LosTopos::Vec2i>& labels = st->m_mesh.get_triangle_labels();

    for(size_t t=0; t < tris.size(); ++t )
    {
        if ( tris[t][0] == tris[t][1] ) { continue; }
        const LosTopos::Vec2i& lbl = labels[t]; 
        std::tuple index =  (lbl[0]<lbl[1])? std::make_tuple(lbl[0], lbl[1]): std::make_tuple(lbl[0], lbl[1]);

        if(areas.find(index)==areas.end()) areas[index]=0.0;
        
        areas[index] += st->get_triangle_area(t);
    }
    return areas;
}

inline std::map<std::tuple<int, int>, double> get_predicted_interface_areas(const  LosTopos::SurfTrack* const st) {
    std::map<std::tuple<int, int>, double> areas;

    const std::vector<LosTopos::Vec3st>& tris = st->m_mesh.get_triangles();
    const std::vector<LosTopos::Vec2i>& labels = st->m_mesh.get_triangle_labels();

    for(size_t t=0; t < tris.size(); ++t )
    {
        if ( tris[t][0] == tris[t][1] ) { continue; }
        const LosTopos::Vec2i& lbl = labels[t]; 
        std::tuple index =  (lbl[0]<lbl[1])? std::make_tuple(lbl[0], lbl[1]): std::make_tuple(lbl[0], lbl[1]);

        if(areas.find(index)==areas.end()) areas[index]=0.0;
        
        const LosTopos::Vec3d &p0 = st->get_newposition(tris[t][0]);
        const LosTopos::Vec3d &p1 = st->get_newposition(tris[t][1]);
        const LosTopos::Vec3d &p2 = st->get_newposition(tris[t][2]);      
        areas[index] += 0.5 * mag(cross(p1-p0, p2-p0));
    }
    return areas;
}

#endif //__UTILS_HPP__