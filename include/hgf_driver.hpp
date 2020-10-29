#ifndef __HGF_DRIVER_HPP__
#define __HGF_DRIVER_HPP__

#include <map>
#include <Eigen/Core>
#include "surftrack.h"

class HGFDriver {
    private:
        LosTopos::SurfTrack* st;
        size_t mesh_update_event_handled;
        std::map<int, double> initial_region_volumes;

        void step_evolution(double dt);
        void step_volume_correction(double dt);

    public:
        HGFDriver(LosTopos::SurfTrack* st, int num_open_regions=1, double surface_tension_strength=0.01);
        void step(double st);
        void update_velocities_after_mesh_defrag();
        void update_velocities_after_mesh_improvement();

        const int NUM_OPEN_REGIONS;
        const double beta;
        Eigen::MatrixXd velocities;
};

#endif //__HGF_DRIVER_HPP__