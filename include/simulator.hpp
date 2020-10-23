//
//  simulator.hpp
//
//  Christopher Batty, Fang Da 2014
//  Ryusuke Sugimoto 2020
//
//

#ifndef __SIMULATOR_HPP__
#define __SIMULATOR_HPP__

#include <string>
#include "surftrack.h"
#include "hgf_driver.hpp"
#include "utils.hpp"

class Simulator : public LosTopos::SurfTrack::SolidVerticesCallback, LosTopos::T1Transition::VelocityFieldCallback {
    public:
        Simulator();
        ~Simulator();
        bool init();
        void step();
        bool isFinished() const { return m_finished; }

    protected:
        void load_scene(std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<LosTopos::Vec3d> & ms, size_t& num_open_regions);
        void stepHGF(double dt);

    protected:
        void updateBBWallConstraints();
        void removeBBWallFaces();
        
        int onBBWall(const LosTopos::Vec3d & pos) const;
        LosTopos::Vec3d enforceBBWallConstraint(const LosTopos::Vec3d & input, int constraints) const;

    protected:
        // SurfTrack::SolidVerticesCallback method
        bool            generate_collapsed_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos);
        bool            generate_split_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos);
        LosTopos::Vec3c generate_collapsed_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1);
        LosTopos::Vec3c generate_split_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1);
        bool            generate_edge_popped_positions(LosTopos::SurfTrack & st, size_t oldv, const LosTopos::Vec2i & cut, LosTopos::Vec3d & pos_upper, LosTopos::Vec3d & pos_lower);
        bool            generate_vertex_popped_positions(LosTopos::SurfTrack & st, size_t oldv, int A, int B, LosTopos::Vec3d & pos_a, LosTopos::Vec3d & pos_b);
        bool            solid_edge_is_feature(const LosTopos::SurfTrack & st, size_t e);
        
        // T1Transition::VelocityFieldCallback methods
        LosTopos::Vec3d sampleVelocity(LosTopos::Vec3d & pos);
        bool sampleDirectionalDivergence(const LosTopos::Vec3d & pos, const LosTopos::Vec3d & dir, double & output);

    protected:
        bool m_verbose;

        bool m_bbwall;
        bool m_t1vel;
        
        double m_dt;
        double m_time;
        int m_frameid;
        bool m_finished;
    
        LosTopos::SurfTrack * m_st;
        HGFDriver* m_driver;

    public:
        double get_dt(){return m_dt;};
        auto V(){return to_eigen_const(m_st->pm_positions);};
        auto F(){return to_eigen_const(m_st->m_mesh.m_tris);};
};

#endif //__SIMULATOR_HPP__