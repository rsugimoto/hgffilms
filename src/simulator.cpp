//
//  simulator.cpp
//
//  Christopher Batty, Fang Da 2014
//  Ryusuke Sugimoto 2020
//
//

#include <any>
#include <unordered_map>
#include <string>
#include "simulator.hpp"
#include "subdivisionscheme.h"
#include "hgf_driver.hpp"

using namespace LosTopos;
std::unordered_map<std::string, std::any> options;

Simulator::Simulator():
    m_verbose(true),
    m_bbwall(false),
    m_t1vel(false),
    m_dt(0),
    m_time(0),
    m_frameid(0),
    m_finished(false),
    m_st(NULL)
{

}

Simulator::~Simulator(){
    delete m_driver;
    delete m_st;
}

bool Simulator::init(){
    options["time-step"] = 0.01;
    options["simulation-time"] = 10.0;
    options["remeshing-resolution"] = 0.05;
    options["remeshing-iterations"] = 3;

    options["lostopos-collision-epsilon-fraction"] = 1e-4;       // lostopos collision epsilon (fraction of mean edge length)
    options["lostopos-merge-proximity-epsilon-fraction"] = 0.02; // lostopos merge proximity epsilon (fraction of mean edge length)
    options["lostopos-perform-smoothing"] = false;               // whether or not to perform smoothing
    options["lostopos-max-volume-change-fraction"] = 100.0;       // maximum allowed volume change during a remeshing operation (fraction of mean edge length cubed)
    options["lostopos-min-triangle-angle"] = 3.0;                // min triangle angle (in degrees)
    options["lostopos-max-triangle-angle"] = 177.0;              // max triangle angle (in degrees)
    options["lostopos-large-triangle-angle-to-split"] = 160.0;   // threshold for large angles to be split
    options["lostopos-min-triangle-area-fraction"] = 0.001;       // minimum allowed triangle area (fraction of mean edge length squared)
    options["lostopos-t1-transition-enabled"] = true;            // whether t1 is enabled
    options["lostopos-t1-pull-apart-distance-fraction"] = 0.1;   // t1 pull apart distance (fraction of mean edge legnth)
    options["lostopos-smooth-subdivision"] = false;              // whether to use smooth subdivision during remeshing
    options["lostopos-allow-non-manifold"] = true;               // whether to allow non-manifold geometry in the mesh
    options["lostopos-allow-topology-changes"] = true;           // whether to allow topology changes

    std::vector<Vec3d> vertices;
    std::vector<Vec3st> faces;
    std::vector<Vec2i> face_labels;
    std::vector<Vec3d> masses;
    size_t num_open_regions;
    load_scene(vertices, faces, face_labels, masses, num_open_regions);

    double mean_edge_len = std::any_cast<double>( options["remeshing-resolution"] );
    double min_edge_len = mean_edge_len * 0.5;
    double max_edge_len = mean_edge_len * 1.5;

    SurfTrackInitializationParameters params;
    params.m_proximity_epsilon = std::any_cast<double>( options["lostopos-collision-epsilon-fraction"] ) * mean_edge_len;
    params.m_merge_proximity_epsilon = std::any_cast<double>( options["lostopos-merge-proximity-epsilon-fraction"] ) * mean_edge_len;
    params.m_allow_vertex_movement_during_collapse = true;
    params.m_perform_smoothing = std::any_cast<bool>( options["lostopos-perform-smoothing"] );
    params.m_min_edge_length = min_edge_len;
    params.m_max_edge_length = max_edge_len;
    params.m_max_volume_change = std::any_cast<double>( options["lostopos-max-volume-change-fraction"] ) * pow(mean_edge_len, 3);
    params.m_min_triangle_angle = std::any_cast<double>( options["lostopos-min-triangle-angle"] );
    params.m_max_triangle_angle = std::any_cast<double>( options["lostopos-max-triangle-angle"] );
    params.m_large_triangle_angle_to_split = std::any_cast<double>( options["lostopos-large-triangle-angle-to-split"] );
    params.m_min_triangle_area = std::any_cast<double>( options["lostopos-min-triangle-area-fraction"] ) * pow(mean_edge_len, 2);
    params.m_verbose = false;
    params.m_allow_non_manifold = std::any_cast<bool>( options["lostopos-allow-non-manifold"] );
    params.m_allow_topology_changes = std::any_cast<bool>( options["lostopos-allow-topology-changes"] );
    params.m_collision_safety = true;
    params.m_remesh_boundaries = true;
    params.m_t1_transition_enabled = std::any_cast<bool>( options["lostopos-t1-transition-enabled"] );
    params.m_pull_apart_distance = std::any_cast<double>( options["lostopos-t1-pull-apart-distance-fraction"] ) * mean_edge_len;

    params.m_velocity_field_callback = NULL;
    if (m_t1vel)
        params.m_velocity_field_callback = this; // this is only turned on for specific scenes

    if ( std::any_cast<bool>( options["lostopos-smooth-subdivision"] ))
        params.m_subdivision_scheme = new ModifiedButterflyScheme();
    else
        params.m_subdivision_scheme = new MidpointScheme();

    params.m_use_curvature_when_collapsing = false;
    params.m_use_curvature_when_splitting = false;

    masses.resize(vertices.size(), Vec3d(1, 1, 1));
    m_st = new SurfTrack(vertices, faces, face_labels, masses, params);
    if (m_bbwall) m_st->m_solid_vertices_callback = this;

    m_driver = new HGFDriver(m_st, num_open_regions);

    // BB wall constraint: update the infinite masses
    updateBBWallConstraints();

    // prepare to start the simulation
    m_time = 0;
    m_dt = std::any_cast<double>( options["time-step"] );
    m_finished = false;

    return true;
}

void Simulator::step() {
    assert(m_st);
    if (m_verbose)
        std::cout << "Time stepping: t = " << m_time << ", dt = " << m_dt << std::endl;

    m_driver->step(m_dt);

    // handle collisions during the dynamics
    double actual_dt;
    m_st->integrate(m_dt, actual_dt);
    if (actual_dt != m_dt)
        std::cout << "Warning: SurfTrack::integrate() failed to step the full length of the time step!" << std::endl;
    
    
    // BB wall constraint: update the infinite masses
    updateBBWallConstraints();
    
    
    // remove faces completely inside the BB wall
    removeBBWallFaces();

    
    // mesh improvement
    for(int i = 0; i < std::any_cast<int>(options["remeshing-iterations"] ); i++)
    {
        m_st->topology_changes();
        m_st->improve_mesh();
    }

    m_driver->update_velocities_after_mesh_improvement();
    
    
    // remove faces completely inside the BB wall
    removeBBWallFaces();
    
    
    // defrag the mesh in the end, to ensure the next step starts with a clean mesh
    m_st->defrag_mesh();

    m_driver->update_velocities_after_mesh_defrag();

    // advance time
    m_frameid++;
    m_time += m_dt;
    if (m_time >= std::any_cast<double>( options["simulation-time"] ))
        m_finished = true;
}

void Simulator::load_scene(std::vector<Vec3d> & vs, std::vector<Vec3st> & fs, std::vector<Vec2i> & ls, std::vector<Vec3d> & ms, size_t& num_open_regions) {
    vs.clear(); vs.reserve(8);
    ms.clear(); ms.reserve(8);
    fs.clear(); fs.reserve(12);
    ls.clear(); ls.reserve(12);
    num_open_regions = 1;
    vs.push_back(Vec3d(0.0, 0.0, 0.0));
    vs.push_back(Vec3d(0.0, 0.0, 1.0));
    vs.push_back(Vec3d(0.0, 1.0, 0.0));
    vs.push_back(Vec3d(0.0, 1.0, 1.0));
    vs.push_back(Vec3d(1.0, 0.0, 0.0));
    vs.push_back(Vec3d(1.0, 0.0, 1.0));
    vs.push_back(Vec3d(1.0, 1.0, 0.0));
    vs.push_back(Vec3d(1.0, 1.0, 1.0));
    for(size_t i = 0; i<8; ++i){
        ms.push_back(Vec3d(1.0, 1.0, 1.0));
    }
    fs.push_back(Vec3st(2, 1, 0));
    fs.push_back(Vec3st(3, 1, 2));
    fs.push_back(Vec3st(6, 0, 4));
    fs.push_back(Vec3st(2, 0, 6));
    fs.push_back(Vec3st(7, 4, 5));
    fs.push_back(Vec3st(6, 4, 7));
    fs.push_back(Vec3st(3, 5, 1));
    fs.push_back(Vec3st(7, 5, 3));
    fs.push_back(Vec3st(7, 2, 6));
    fs.push_back(Vec3st(3, 2, 7));
    fs.push_back(Vec3st(4, 1, 5));
    fs.push_back(Vec3st(0, 1, 4));
    for(size_t i = 0; i<12; ++i){
        ls.push_back(Vec2i(0, 1));
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Utilities for maintaining a [0, 1]^3 bounding box
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Simulator::updateBBWallConstraints()
{
    if (m_bbwall)
    {
        for (size_t i = 0; i < m_st->m_mesh.nv(); i++)
        {
            Vec3d & x = m_st->pm_positions[i];
            if (x[0] > 1) x[0] = 1;
            if (x[0] < 0) x[0] = 0;
            if (x[1] > 1) x[1] = 1;
            if (x[1] < 0) x[1] = 0;
            if (x[2] > 1) x[2] = 1;
            if (x[2] < 0) x[2] = 0;
            
            int onwall = onBBWall(x);
            Vec3d & mass = m_st->m_masses[i];
            mass[0] = ((onwall & (1 << 0)) || (onwall & (1 << 3)) ? std::numeric_limits<double>::infinity() : 1);
            mass[1] = ((onwall & (1 << 1)) || (onwall & (1 << 4)) ? std::numeric_limits<double>::infinity() : 1);
            mass[2] = ((onwall & (1 << 2)) || (onwall & (1 << 5)) ? std::numeric_limits<double>::infinity() : 1);
        }
    }
}

void Simulator::removeBBWallFaces()
{
    if (!m_bbwall) // scene doesn't use the BB
        return;

    // remove all faces completely inside wall
    for (size_t i = 0; i < m_st->m_mesh.nt(); i++)
    {
        Vec3st & f = m_st->m_mesh.m_tris[i];
        int w0 = onBBWall(m_st->get_position(f[0]));
        int w1 = onBBWall(m_st->get_position(f[1]));
        int w2 = onBBWall(m_st->get_position(f[2]));
        
        if (w0 & w1 & w2)
            m_st->remove_triangle(i);
    }
}

int Simulator::onBBWall(const Vec3d & pos) const
{
    if (!m_bbwall) // scene doesn't use the BB
        return 0;
    
    static const double WALL_THRESHOLD = 1e-6;
    
    int walls = 0;
    if (pos[0] < 0 + WALL_THRESHOLD) walls |= (1 << 0);
    if (pos[1] < 0 + WALL_THRESHOLD) walls |= (1 << 1);
    if (pos[2] < 0 + WALL_THRESHOLD) walls |= (1 << 2);
    if (pos[0] > 1 - WALL_THRESHOLD) walls |= (1 << 3);
    if (pos[1] > 1 - WALL_THRESHOLD) walls |= (1 << 4);
    if (pos[2] > 1 - WALL_THRESHOLD) walls |= (1 << 5);
    
    return walls;
}

Vec3d Simulator::enforceBBWallConstraint(const Vec3d & input, int constraints) const
{
    if (!m_bbwall) // scene doesn't use the BB
        return input;
    
    Vec3d output = input;
    if (constraints & (1 << 0)) output[0] = 0;
    if (constraints & (1 << 1)) output[1] = 0;
    if (constraints & (1 << 2)) output[2] = 0;
    if (constraints & (1 << 3)) output[0] = 1;
    if (constraints & (1 << 4)) output[1] = 1;
    if (constraints & (1 << 5)) output[2] = 1;
    
    return output;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Callbacks
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Simulator::generate_collapsed_position(SurfTrack & st, size_t v0, size_t v1, Vec3d & pos)
{
    Vec3d x0 = st.get_position(v0);
    Vec3d x1 = st.get_position(v1);
    
    int label0 = onBBWall(Vec3d(x0[0], x0[1], x0[2]));
    int label1 = onBBWall(Vec3d(x1[0], x1[1], x1[2]));
    
    if (label0 == label1)
    {
        // on the same wall(s), prefer the one with higher max edge valence
        size_t maxedgevalence0 = 0;
        size_t maxedgevalence1 = 0;
        for (size_t i = 0; i < st.m_mesh.m_vertex_to_edge_map[v0].size(); i++)
            if (st.m_mesh.m_edge_to_triangle_map[st.m_mesh.m_vertex_to_edge_map[v0][i]].size() > maxedgevalence0)
                maxedgevalence0 = st.m_mesh.m_edge_to_triangle_map[st.m_mesh.m_vertex_to_edge_map[v0][i]].size();
        for (size_t i = 0; i < st.m_mesh.m_vertex_to_edge_map[v1].size(); i++)
            if (st.m_mesh.m_edge_to_triangle_map[st.m_mesh.m_vertex_to_edge_map[v1][i]].size() > maxedgevalence1)
                maxedgevalence1 = st.m_mesh.m_edge_to_triangle_map[st.m_mesh.m_vertex_to_edge_map[v1][i]].size();
        
        if (maxedgevalence0 == maxedgevalence1) // same max edge valence, use their midpoint
            pos = (x0 + x1) / 2;
        else if (maxedgevalence0 < maxedgevalence1)
            pos = x1;
        else
            pos = x0;
        
        return true;
        
    } else if ((label0 & ~label1) == 0)
    {
        // label0 is a proper subset of label1 (since label0 != label1)
        pos = x1;
        
        return true;
        
    } else if ((label1 & ~label0) == 0)
    {
        // label1 is a proper subset of label0
        pos = x0;
        
        return true;
        
    } else
    {
        // label0 and label1 are not subset of each other
        int newlabel = label0 | label1;
        assert(label0 != newlabel); // not subset of each other
        assert(label1 != newlabel);
        
        assert(!((label0 & (1 << 0)) != 0 && (label0 & (1 << 3)) != 0)); // can't have conflicting constraints in label0 and label1 already
        assert(!((label0 & (1 << 1)) != 0 && (label0 & (1 << 4)) != 0));
        assert(!((label0 & (1 << 2)) != 0 && (label0 & (1 << 5)) != 0));
        assert(!((label1 & (1 << 0)) != 0 && (label1 & (1 << 3)) != 0));
        assert(!((label1 & (1 << 1)) != 0 && (label1 & (1 << 4)) != 0));
        assert(!((label1 & (1 << 2)) != 0 && (label1 & (1 << 5)) != 0));
        
        bool conflict = false;
        if ((newlabel & (1 << 0)) != 0 && (newlabel & (1 << 3)) != 0) conflict = true;
        if ((newlabel & (1 << 1)) != 0 && (newlabel & (1 << 4)) != 0) conflict = true;
        if ((newlabel & (1 << 2)) != 0 && (newlabel & (1 << 5)) != 0) conflict = true;
        
        if (conflict)
        {
            // the two vertices are on opposite walls (conflicting constraints). Can't collapse this edge (which shouldn't have become a collapse candidate in the first place)
            return false;
        }
        
        pos = (x0 + x1) / 2;
        if (newlabel & (1 << 0))  pos[0] = 0;    // project the midpoint onto the constraint manifold (BB walls)
        if (newlabel & (1 << 1))  pos[1] = 0;
        if (newlabel & (1 << 2))  pos[2] = 0;
        if (newlabel & (1 << 3))  pos[0] = 1;
        if (newlabel & (1 << 4))  pos[1] = 1;
        if (newlabel & (1 << 5))  pos[2] = 1;
        
        return true;
    }
    
    return false;
}

bool Simulator::generate_split_position(SurfTrack & st, size_t v0, size_t v1, Vec3d & pos)
{
    pos = (st.get_position(v0) + st.get_position(v1)) / 2;
    
    return true;
}

Vec3c Simulator::generate_collapsed_solid_label(SurfTrack & st, size_t v0, size_t v1, const Vec3c & label0, const Vec3c & label1)
{
    Vec3d x0 = st.get_position(v0);
    Vec3d x1 = st.get_position(v1);
    
    int constraint0 = onBBWall(Vec3d(x0[0], x0[1], x0[2]));
    int constraint1 = onBBWall(Vec3d(x1[0], x1[1], x1[2]));
    
    assert(((constraint0 & (1 << 0)) || (constraint0 & (1 << 3))) == (bool)label0[0]);
    assert(((constraint0 & (1 << 1)) || (constraint0 & (1 << 4))) == (bool)label0[1]);
    assert(((constraint0 & (1 << 2)) || (constraint0 & (1 << 5))) == (bool)label0[2]);
    assert(((constraint1 & (1 << 0)) || (constraint1 & (1 << 3))) == (bool)label1[0]);
    assert(((constraint1 & (1 << 1)) || (constraint1 & (1 << 4))) == (bool)label1[1]);
    assert(((constraint1 & (1 << 2)) || (constraint1 & (1 << 5))) == (bool)label1[2]);
    
    Vec3c result;  // if either endpoint is constrained, the collapsed point shold be constrained. more specifically it should be on all the walls any of the two endpoints is on (implemented in generate_collapsed_position())
    int result_constraint = (constraint0 | constraint1);
    result[0] = ((result_constraint & (1 << 0)) || (result_constraint & (1 << 3)));
    result[1] = ((result_constraint & (1 << 1)) || (result_constraint & (1 << 4)));
    result[2] = ((result_constraint & (1 << 2)) || (result_constraint & (1 << 5)));
    
    return result;
}

Vec3c Simulator::generate_split_solid_label(SurfTrack & st, size_t v0, size_t v1, const Vec3c & label0, const Vec3c & label1)
{
    Vec3d x0 = st.get_position(v0);
    Vec3d x1 = st.get_position(v1);
    
    int constraint0 = onBBWall(Vec3d(x0[0], x0[1], x0[2]));
    int constraint1 = onBBWall(Vec3d(x1[0], x1[1], x1[2]));
    
    assert(((constraint0 & (1 << 0)) || (constraint0 & (1 << 3))) == (bool)label0[0]);
    assert(((constraint0 & (1 << 1)) || (constraint0 & (1 << 4))) == (bool)label0[1]);
    assert(((constraint0 & (1 << 2)) || (constraint0 & (1 << 5))) == (bool)label0[2]);
    assert(((constraint1 & (1 << 0)) || (constraint1 & (1 << 3))) == (bool)label1[0]);
    assert(((constraint1 & (1 << 1)) || (constraint1 & (1 << 4))) == (bool)label1[1]);
    assert(((constraint1 & (1 << 2)) || (constraint1 & (1 << 5))) == (bool)label1[2]);
  
    Vec3c result;  // the splitting midpoint has a positive constraint label only if the two endpoints are on a same wall (sharing a bit in their constraint bitfield representation)
    int result_constraint = (constraint0 & constraint1);
    result[0] = ((result_constraint & (1 << 0)) || (result_constraint & (1 << 3)));
    result[1] = ((result_constraint & (1 << 1)) || (result_constraint & (1 << 4)));
    result[2] = ((result_constraint & (1 << 2)) || (result_constraint & (1 << 5)));
    
    return result;
    
}

bool Simulator::generate_edge_popped_positions(SurfTrack & st, size_t oldv, const Vec2i & cut, Vec3d & pos_upper, Vec3d & pos_lower)
{
    Vec3d original_pos = st.get_position(oldv);
    int original_constraint = onBBWall(Vec3d(original_pos[0], original_pos[1], original_pos[2]));
    
    Vec3d new_pos_upper = enforceBBWallConstraint(Vec3d(pos_upper[0], pos_upper[1], pos_upper[2]), original_constraint);
    Vec3d new_pos_lower = enforceBBWallConstraint(Vec3d(pos_lower[0], pos_lower[1], pos_lower[2]), original_constraint);
    
    pos_upper = new_pos_upper;
    pos_lower = new_pos_lower;
    
    return true;
}

bool Simulator::generate_vertex_popped_positions(SurfTrack & st, size_t oldv, int A, int B, Vec3d & pos_a, Vec3d & pos_b)
{
    Vec3d original_pos = st.get_position(oldv);
    int original_constraint = onBBWall(Vec3d(original_pos[0], original_pos[1], original_pos[2]));
    
    Vec3d new_pos_a = enforceBBWallConstraint(Vec3d(pos_a[0], pos_a[1], pos_a[2]), original_constraint);
    Vec3d new_pos_b = enforceBBWallConstraint(Vec3d(pos_b[0], pos_b[1], pos_b[2]), original_constraint);
    
    pos_a = new_pos_a;
    pos_b = new_pos_b;
    
    return true;
}

bool Simulator::solid_edge_is_feature(const SurfTrack & st, size_t e)
{
    Vec3d x0 = st.get_position(st.m_mesh.m_edges[e][0]);
    Vec3d x1 = st.get_position(st.m_mesh.m_edges[e][1]);
    
    int constraint0 = onBBWall(Vec3d(x0[0], x0[1], x0[2]));
    int constraint1 = onBBWall(Vec3d(x1[0], x1[1], x1[2]));
    
    if (constraint0 & constraint1)  // edge is completely inside a wall
        return true;
    else
        return false;
}

Vec3d Simulator::sampleVelocity(Vec3d & pos)
{
    return Vec3d(0, 0, 0);
}

bool Simulator::sampleDirectionalDivergence(const Vec3d & pos, const Vec3d & dir, double & output)
{
    return false;
}