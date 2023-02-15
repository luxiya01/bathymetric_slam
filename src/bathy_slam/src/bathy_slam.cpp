#include "bathy_slam/bathy_slam.hpp"

BathySlam::BathySlam(GraphConstructor &graph_obj, SubmapRegistration &gicp_reg, benchmark::track_error_benchmark& benchmark):
    graph_obj_(&graph_obj), gicp_reg_(&gicp_reg), benchmark_(&benchmark) {

}

BathySlam::~BathySlam(){

}

void BathySlam::saveLCtoText(SubmapObj& submap_i, ofstream& fileOutputStream) {
    if(fileOutputStream.is_open()){
        fileOutputStream << submap_i.submap_id_;
        for(unsigned int j=0; j<submap_i.overlaps_idx_.size(); j++){
            fileOutputStream << " " << submap_i.overlaps_idx_.at(j);
        }
        fileOutputStream << "\n";
    }
}

SubmapObj BathySlam::findLoopClosureByCombiningOverlappingSubmaps(
    SubmapObj& submap_i,
    ofstream& fileOutputStream,
    SubmapObj& submap_trg,
    SubmapsVec& submaps_reg,
    DRNoise& dr_noise,
    YAML::Node& config,
    GaussianGen& transSampler,
    GaussianGen& rotSampler,
    double info_thres
    ) {
    // If potential loop closure detected
    SubmapObj submap_final = submap_i;

    if(!submap_i.overlaps_idx_.empty()){
        // Save loop closure to txt
        saveLCtoText(submap_i, fileOutputStream);

        // Register overlapping submaps
        submap_trg = gicp_reg_->constructTrgSubmap(submaps_reg, submap_i.overlaps_idx_, dr_noise);
        if (config["add_gaussian_noise"].as<bool>()) {
            addNoiseToSubmap(transSampler, rotSampler, submap_i); // Add disturbance to source submap
        }

        if(gicp_reg_->gicpSubmapRegistration(submap_trg, submap_i)){
            submap_final = submap_i;
        }
        submap_trg.submap_pcl_.clear();

        // Create loop closures
        graph_obj_->edge_covs_type_ = config["lc_edge_covs_type"].as<int>();
        graph_obj_->findLoopClosures(submap_final, submaps_reg, info_thres);
    }
    return submap_final;
}

Isometry3f BathySlam::loadLCFromFile(SubmapObj& submap, const string& filepath, const YAML::Node& config) {
    YAML::Node lc_file;
    Isometry3f identity = Isometry3f::Identity();

    try {
        lc_file = YAML::LoadFile(filepath);
    } catch (const std::exception& e) {
        cerr << "Failed to load LC file: " << filepath << endl;
        cerr << e.what() << endl;
        return identity;
    }

    // Do not transform the submap if the final RMS is larger than the given threshold
    double final_rms_error = lc_file["rms_error"]["final"].as<double>();
    if (final_rms_error >= config["offline_LC_max_rms"].as<double>()) {
        cout << "Ignore LC with large rms_error " << final_rms_error << ": " << filepath << endl;
        return identity;
    }
    auto tf = lc_file["transform"].as<std::vector<float>>();
    Matrix4f transform = Eigen::Map<Eigen::Matrix4f>(tf.data());
    Isometry3f rel_tf = Isometry3f (Isometry3f(Translation3f(transform.block<3,1>(0,3)))*
                                    Isometry3f(Quaternionf(transform.block<3,3>(0,0)).normalized()));


    submap.submap_tf_ = rel_tf * submap.submap_tf_;
    pcl::transformPointCloud(submap.submap_pcl_, submap.submap_pcl_, rel_tf);

    // Set GICP information matrix
    submap.submap_lc_info_.setZero();
    //auto info_diag = Eigen::Map<Eigen::Vector4d>(config["gicp_info_diag"].as<std::vector<double>>().data());
    auto gicp_info_diag = config["gicp_info_diag"].as<std::vector<double>>();
    auto info_diag = Eigen::Map<Eigen::Vector4d>(gicp_info_diag.data());
    submap.submap_lc_info_.bottomRightCorner(4, 4) = info_diag.asDiagonal();

    double xy_edge_info = 1/std::pow(final_rms_error, 2);
    Matrix2d lc_info_xy{
       {xy_edge_info, 0},
       {0, xy_edge_info}
    };
    submap.submap_lc_info_.topLeftCorner(2, 2) = lc_info_xy;
    return rel_tf;
}

SubmapsVec BathySlam::runOffline(SubmapsVec& submaps_gt, GaussianGen& transSampler, GaussianGen& rotSampler, YAML::Node config){
    DRNoise dr_noise = loadDRNoiseFromFile(config);
    SubmapObj submap_trg(dr_noise);
    SubmapsVec submaps_prev, submaps_reg;
    ofstream fileOutputStream;
    fileOutputStream.open("loop_closures.txt", std::ofstream::out);

    for(SubmapObj& submap_i: submaps_gt){
        std::cout << " ----------- Submap " << submap_i.submap_id_ << ", swath "
                  << submap_i.swath_id_ << " ------------"<< std::endl;

        // Look for loop closures
        for(SubmapObj& submap_k: submaps_reg){
            // Don't look for overlaps between submaps of the same swath or the prev submap
            if(submap_k.submap_id_ != submap_i.submap_id_ - 1){
                submaps_prev.push_back(submap_k);
            }
        }
        // Submaps in map_frame?
        bool submaps_in_map_tf = true;
            submap_i.findOverlaps(submaps_in_map_tf, submaps_prev, config["overlap_coverage"].as<double>());
            submaps_prev.clear();

    #if INTERACTIVE == 1
        // Update visualizer
        submaps_reg.push_back(submap_i); // Add submap_i to registered set (just for visualization here)
        visualizer->updateVisualizer(submaps_reg);
        while(!viewer.wasStopped ()){
            viewer.spinOnce ();
        }
        viewer.resetStoppedFlag();
        submaps_reg.pop_back();
    #endif
        // Create graph vertex i
        graph_obj_->createNewVertex(submap_i);

        // Create DR edge i and store (skip submap 0)
        if(submap_i.submap_id_ != 0 ){
            std::cout << "DR edge from " << submap_i.submap_id_ -1 << " to " << submap_i.submap_id_<< std::endl;
            graph_obj_->createDREdge(submap_i);
        }

        SubmapObj submap_final(dr_noise);
        submap_final = submap_i;
        if (config["online_LC_combine_overlapping_submaps"].as<bool>()) {
            submap_final = findLoopClosureByCombiningOverlappingSubmaps(
                submap_i, fileOutputStream, submap_trg, submaps_reg,
                dr_noise, config, transSampler, rotSampler);
        } else {
            if (!submap_i.overlaps_idx_.empty()) {
                saveLCtoText(submap_i, fileOutputStream);

                for (SubmapObj submap_j : submaps_reg) {
                    if(find(submap_i.overlaps_idx_.begin(), submap_i.overlaps_idx_.end(), submap_j.submap_id_)
                            != submap_i.overlaps_idx_.end()){
                        Isometry3f rel_tf = Isometry3f::Identity();
                        if (config["load_offline_LC"].as<bool>()) {
                            string lc_filepath = config["offline_LC_folder"].as<string>()
                                                + "submap_" + to_string(submap_i.submap_id_)
                                                + "-" + "submap_" + to_string(submap_j.submap_id_)
                                                + ".yaml";
                            rel_tf = loadLCFromFile(submap_final, lc_filepath, config);
                        } else {
                            // Compute GICP for each individual submap pairs online
                        }
                        graph_obj_->edge_covs_type_ = config["lc_edge_covs_type"].as<int>();
                        graph_obj_->createLCEdge(submap_final, submap_j);

                        // Transform the submap pose & point cloud TF back to before the LC (i.e. to the DR tf)
                        submap_final.submap_tf_ = rel_tf.inverse() * submap_final.submap_tf_;
                        pcl::transformPointCloud(submap_final.submap_pcl_, submap_final.submap_pcl_, rel_tf.inverse());
                        submap_final.submap_lc_info_.setZero();
                    }
                }
            }
        }
        submaps_reg.push_back(submap_final);    // Add registered submap_i

    #if INTERACTIVE == 1
        // Update visualizer
        visualizer->updateVisualizer(submaps_reg);
        while(!viewer.wasStopped ()){
            viewer.spinOnce ();
        }
        viewer.resetStoppedFlag();
    #endif
    }
    fileOutputStream.close();


    return submaps_reg;
}
