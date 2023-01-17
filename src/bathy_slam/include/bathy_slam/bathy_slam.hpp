#ifndef BATHY_SLAM_HPP
#define BATHY_SLAM_HPP

#include "submaps_tools/submaps.hpp"

#include "graph_optimization/graph_construction.hpp"
#include "graph_optimization/utils_g2o.hpp"

#include "registration/gicp_reg.hpp"

class BathySlam{

public:

    GraphConstructor* graph_obj_;
    SubmapRegistration* gicp_reg_;

    BathySlam(GraphConstructor& graph_obj, SubmapRegistration& gicp_reg);
    ~BathySlam();

    SubmapsVec runOffline(SubmapsVec &submaps_gt, GaussianGen &transSampler, GaussianGen &rotSampler, YAML::Node config);
    SubmapObj findLoopClosureByCombiningOverlappingSubmaps(
        SubmapObj& submap_i,
        ofstream& fileOutputStream,
        SubmapObj& submap_trg,
        SubmapsVec& submaps_reg,
        DRNoise& dr_noise,
        YAML::Node& config,
        GaussianGen& transSampler,
        GaussianGen& rotSampler,
        double info_thres=0.1);
};


#endif //BATHY_SLAM_HPP
