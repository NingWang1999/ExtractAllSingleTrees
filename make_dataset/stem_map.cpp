/*
 * Software License Agreement (Apache License)
 *
 *  Copyright (C) 2023, Xufei Wang (tjwangxufei@tongji.edu.cn),
 *                      Zexin Yang (zexinyang@tongji.edu.cn),
 *                      Liangliang Nan (liangliang.nan@gmail.com).
 *  All rights reserved.
 *
 *  This file is part of GlobalMatch (https://github.com/zexinyang/GlobalMatch),
 *  which implements the point cloud registration method described in the following paper:
 *  -----------------------------------------------------------------------------------------------------------
 *  GlobalMatch: Registration of forest terrestrial point clouds by global matching of relative stem positions.
 *  Xufei Wang, Zexin Yang, Xiaojun Cheng, Jantien Stoter, Wenbing Xu, Zhenlun Wu, and Liangliang Nan.
 *  ISPRS Journal of Photogrammetry and Remote Sensing. Vol. 197, 71-86, 2023.
 *  -----------------------------------------------------------------------------------------------------------
 *  We kindly ask you to cite the above paper if you use (part of) the code or ideas in your academic work.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 */
#include "stem_map.h"
 // pcl
#include <pcl/common/angles.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/io/pcd_io.h>
#include <pcl/surface/mls.h>
#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>
#include <pcl/point_cloud.h>
#include <pcl/filters/radius_outlier_removal.h>
// remove the leaf-size check
#include "voxel_grid_fix.h"
#include "octree_extract_clusters.h"
// csf
#include "CSF.h"
// cc
#include "jacobi.h"
// ICRA 2015 octree
#include "octree_unibn.hpp"
// vtk
#include <vtkOBBTree.h>
// savefile
#include <filesystem>

void
transformCylinderAxisToTwoPoints(const Coefficient& coeff,
    const double z_bottom,
    const double z_top,
    double(&pt_bottom)[3],
    double(&pt_top)[3]) {
    // transform the point-slope line coefficients to two 3D points
    const auto& x0 = coeff.values[0];
    const auto& y0 = coeff.values[1];
    const auto& z0 = coeff.values[2];
    const auto& dx = coeff.values[3];
    const auto& dy = coeff.values[4];
    const auto& dz = coeff.values[5];
    pt_bottom[0] = (z_bottom - z0) * dx / dz + x0; // x_bottom
    pt_bottom[1] = (z_bottom - z0) * dy / dz + y0; // y_bottom
    pt_bottom[2] = z_bottom;
    pt_top[0] = (z_top - z0) * dx / dz + x0; // x_top
    pt_top[1] = (z_top - z0) * dy / dz + y0; // y_top
    pt_top[2] = z_top;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
fitOneStem(const Cloud3D::Ptr& cloud_stem,
    int min_inr_per_cylinder,
    pcl::ModelCoefficients::Ptr& coefficients_cylinder) {

    CloudpointNormal::Ptr clouds_normal(new CloudpointNormal);
    pcl::search::KdTree<Point3D>::Ptr tree(new pcl::search::KdTree<Point3D>);
    tree->setInputCloud(cloud_stem);
    pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mls;
    mls.setComputeNormals(true);
    mls.setInputCloud(cloud_stem);
    mls.setPolynomialOrder(2);
    mls.setSearchMethod(tree);
    mls.setSearchRadius(0.03);
    mls.process(*clouds_normal);
  
    Cloud3D::Ptr cloud_trunk(new Cloud3D);
    CloudNormal::Ptr normal_trunk(new CloudNormal);
    pcl::copyPointCloud(*clouds_normal, *cloud_trunk);
    pcl::copyPointCloud(*clouds_normal, *normal_trunk);
    pcl::PointIndices::Ptr inliers(new pcl::PointIndices);

    pcl::SACSegmentationFromNormals<pcl::PointXYZ, pcl::Normal> seg;
    seg.setOptimizeCoefficients(true);
    seg.setModelType(pcl::SACMODEL_CYLINDER);
    seg.setMethodType(pcl::SAC_RANSAC);
    seg.setNormalDistanceWeight(0.1);
    seg.setMaxIterations(10000);
    seg.setDistanceThreshold(0.03);
    seg.setRadiusLimits(0.01, 0.2);
    seg.setInputCloud(cloud_trunk);
    seg.setInputNormals(normal_trunk);
    seg.segment(*inliers, *coefficients_cylinder);
    
    double inlier_ratio = (cloud_stem->size() > 0)
        ? static_cast<double>(inliers->indices.size()) / cloud_stem->size()
        : 0.0;
    if (min_inr_per_cylinder < inlier_ratio && inlier_ratio < 1 && inliers->indices.size() >= 10) {//FIXME
        return true;
    }
    else 
        return false;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<float>
calculateVerticality(
    const Cloud3D::Ptr& cloud_understory,
    float search_radius) {

    size_t cloud_size = cloud_understory->size();
    unibn::Octree<Point3D> octree;
    octree.initialize(*cloud_understory);
    std::vector<std::vector<uint32_t>> neighbors(cloud_size);
    std::vector<float> verticality(cloud_size, 0.0);

    // Query neighbors and calculate verticality in parallel
#if defined(_OPENMP)
#pragma omp parallel
#endif
    {
        // Query neighbors
#if defined(_OPENMP)
#pragma omp for
#endif
        for (int i = 0; i < cloud_size; ++i) {
            octree.radiusNeighbors<unibn::L2Distance<Point3D>>(
                cloud_understory->points[i], search_radius, neighbors[i]);
        }

        // Calculate verticality
#if defined(_OPENMP)
#pragma omp for
#endif
        for (int i = 0; i < neighbors.size(); ++i) {
            // Need at least 3 neighbors (including the point itself)
            if (neighbors[i].size() < 4)
                continue;

            // Calculate the gravity center of the neighbors
            const auto& count = neighbors[i].size();
            float gravity_center[3] = { 0, 0, 0 };
            for (auto& neighbor : neighbors[i]) {
                gravity_center[0] += cloud_understory->points[neighbor].x;
                gravity_center[1] += cloud_understory->points[neighbor].y;
                gravity_center[2] += cloud_understory->points[neighbor].z;
            }
            gravity_center[0] /= count;
            gravity_center[1] /= count;
            gravity_center[2] /= count;

            // Build the covariance matrix
            float mxx = 0.0, myy = 0.0, mzz = 0.0, mxy = 0.0, mxz = 0.0, myz = 0.0;
            for (auto& neighbor : neighbors[i]) {
                float dx = cloud_understory->points[neighbor].x - gravity_center[0];
                float dy = cloud_understory->points[neighbor].y - gravity_center[1];
                float dz = cloud_understory->points[neighbor].z - gravity_center[2];
                mxx += dx * dx;
                myy += dy * dy;
                mzz += dz * dz;
                mxy += dx * dy;
                mxz += dx * dz;
                myz += dy * dz;
            }

            CCCoreLib::SquareMatrixf mat_cov(3);
            mat_cov.m_values[0][0] = mxx / count;
            mat_cov.m_values[1][1] = myy / count;
            mat_cov.m_values[2][2] = mzz / count;
            mat_cov.m_values[1][0] = mat_cov.m_values[0][1] = mxy / count;
            mat_cov.m_values[2][0] = mat_cov.m_values[0][2] = mxz / count;
            mat_cov.m_values[2][1] = mat_cov.m_values[1][2] = myz / count;

            CCCoreLib::SquareMatrixf eigen_vectors;
            std::vector<float> eigen_values;
            if (!CCCoreLib::Jacobi<float>::ComputeEigenValuesAndVectors(
                mat_cov, eigen_vectors, eigen_values, true))
                continue;

            // Sort eigenvectors by eigenvalues in decreasing order
            CCCoreLib::Jacobi<float>::SortEigenValuesAndVectors(eigen_vectors, eigen_values);

            CCVector3f z(0, 0, 1);
            CCVector3f e3(z);
            CCCoreLib::Jacobi<float>::GetEigenVector(eigen_vectors, 2, e3.u);

            // Verticality is calculated based on the eigenvectors
            verticality[i] = 1.0f - std::abs(z.dot(e3));
        }
    }

    return verticality;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Mapping::extractUnderstory() {
    std::vector<int> indices_extracted, indices_remained;
    CSF csf;
    csf.params.bSloopSmooth = true;
    csf.params.dist_max = 0.6;//1.0, types on 0.5m
    csf.params.dist_min = 0.15;//0.2
    csf.params.cloth_resolution = 0.5;
    csf.params.interations = 500;
    csf.params.rigidness = 2;
    csf.setPointCloud(cloud_input_);
    csf.do_filtering(indices_extracted, indices_remained, true, mesh_ground_);
    indices_understory_->indices = indices_extracted;

    std::cout << "extract understory success" << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Mapping::extractStemPoints() {
    // downsampling
    Cloud3D::Ptr cloud_understory(new Cloud3D);
    pcl::VoxelGrid<Point3D> vg;
    vg.setInputCloud(cloud_input_);
    vg.setIndices(indices_understory_);
    vg.setLeafSize(leaf_size_, leaf_size_, leaf_size_);
    vg.filter(*cloud_understory);

    //statiscal_removal
    Cloud3D::Ptr statiscal_filtered1(new Cloud3D);
    pcl::StatisticalOutlierRemoval<pcl::PointXYZ> sor1;
    sor1.setInputCloud(cloud_understory);
    sor1.setMeanK(5);
    sor1.setStddevMulThresh(0.1);
    sor1.filter(*statiscal_filtered1);    
    cloud_understory = statiscal_filtered1;

    pcl::io::savePCDFileASCII("cloud_understory1.pcd", *cloud_understory);
    std::cout << "downsampling success" << std::endl;

    // Calculate verticality for each point
    std::vector<float> verticalities = calculateVerticality(cloud_understory, search_radius_);

    for (int i = 0; i < verticalities.size(); ++i) {
        if (verticalities[i] >= verticality_threshold_)
            cloud_stems_->push_back(cloud_understory->points[i]);
    }

    pcl::io::savePCDFileASCII("cloud_understory2.pcd", *cloud_stems_);

    // Calculate the normals of points and filter again 
    std::vector<float> verticalities_z(cloud_stems_->size(), 0.0f);

    pcl::NormalEstimation<Point3D, pcl::Normal> ne;
    pcl::PointCloud<pcl::Normal>::Ptr cloud_normals(new pcl::PointCloud<pcl::Normal>());
    ne.setInputCloud(cloud_stems_);
    ne.setKSearch(50);
    ne.compute(*cloud_normals);

    // Determine the angle between the normal and the Z axisz
    for (size_t i = 0; i < cloud_stems_->size(); ++i) {
        const auto& normal = cloud_normals->points[i];
        float verticality_z = std::abs(normal.normal_z);
        verticalities_z[i] = verticality_z;
    }

    Cloud3D::Ptr refined_cluster(new Cloud3D);

    for (size_t i = 0; i < cloud_stems_->size(); ++i) {
        if (verticalities_z[i] <= 0.5) {
            refined_cluster->push_back(cloud_stems_->points[i]);
        }
    }

    //Filter again according to the point spacing to eliminate the point cloud with low density
    pcl::RadiusOutlierRemoval<Point3D> radius_filter;
    radius_filter.setInputCloud(refined_cluster);
    radius_filter.setRadiusSearch(0.1);  
    radius_filter.setMinNeighborsInRadius(100); 
    radius_filter.filter(*refined_cluster);

    pcl::io::savePCDFileASCII("cloud_understory3.pcd", *refined_cluster);


    //statiscal_removal
    Cloud3D::Ptr statiscal_filtered2(new Cloud3D);
    pcl::StatisticalOutlierRemoval<pcl::PointXYZ> sor2;
    sor2.setInputCloud(refined_cluster);
    sor2.setMeanK(10);
    sor2.setStddevMulThresh(0.1);
    sor2.filter(*statiscal_filtered2);

    cloud_stems_->clear();
    cloud_stems_ = statiscal_filtered2;

    std::cout << "statiscal_removal success" << std::endl;
    pcl::io::savePCDFileASCII("cloud_understory4.pcd", *cloud_stems_);

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Mapping::extractTreePositions() {
    // Euclidean cluster extraction
    std::vector<pcl::PointIndices> ptid_clusters;
    OctreeEuclideanClusterExtraction<Point3D> oec;
    oec.setClusterTolerance(min_dist_between_points_);//0.05
    oec.setMinClusterSize(min_pts_per_cluster_);
    oec.setInputCloud(cloud_stems_);
    oec.extract(ptid_clusters);

    std::cout << "there are " << ptid_clusters.size() << " initial clusteres." << std::endl;

    // Cylinder fitting
    int j = 0;
    CoefficientVector coefficients(ptid_clusters.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < ptid_clusters.size(); ++i) {
        Cloud3D::Ptr cloud_cluster(new Cloud3D);
        for (int& pt_index : ptid_clusters[i].indices)
            cloud_cluster->push_back((*cloud_stems_)[pt_index]);

        pcl::ModelCoefficients::Ptr coefficients_cylinder(new pcl::ModelCoefficients);
        if (fitOneStem(cloud_cluster, 0.3, coefficients_cylinder)) {

            coefficients[i] = *coefficients_cylinder;
            j++;
        }
    }
    std::cout << "cylinder has " << j << std::endl;
    std::cout << "Cylinder fitting success" << std::endl;


    // Extract candidate positions
    // calculate two given z values
    pcl::PointCloud<pcl::PointXYZL>::Ptr candidate_positions(new pcl::PointCloud<pcl::PointXYZL>);
    double bbox[6]; // x_min, x_max, y_min, y_max, z_min, z_max
    mesh_ground_->GetBounds(bbox);
    double z_bottom = bbox[4] - 50.0; // z_min - rough guess
    double z_top = bbox[5] + 50.0; // z_max + rough guess
    vtkSmartPointer<vtkOBBTree> obbtree = vtkSmartPointer<vtkOBBTree>::New();
    obbtree->SetDataSet(mesh_ground_);
    obbtree->BuildLocator();
    for (int id_cluster = 0; id_cluster < ptid_clusters.size(); ++id_cluster) {
        if (coefficients[id_cluster].values.empty())
            continue;
        // transform the point-slope line coefficients to two 3D points
        double pt_bottom[3], pt_top[3];
        transformCylinderAxisToTwoPoints(coefficients[id_cluster], z_bottom, z_top, pt_bottom, pt_top);
        vtkSmartPointer<vtkPoints> intersections = vtkSmartPointer<vtkPoints>::New();
        obbtree->IntersectWithLine(pt_bottom, pt_top, intersections, nullptr);
        if (intersections->GetNumberOfPoints() == 1) {
            double intersection[3];
            intersections->GetPoint(0, intersection);
            pcl::PointXYZL position_with_cluster_id;
            position_with_cluster_id.x = static_cast<float>(intersection[0]);
            position_with_cluster_id.y = static_cast<float>(intersection[1]);
            position_with_cluster_id.z = static_cast<float>(intersection[2]);
            position_with_cluster_id.label = id_cluster;
            candidate_positions->push_back(position_with_cluster_id);
        }
    }

    std::cout << "extract " << candidate_positions->size() << " positions" << std::endl;

    // Optimize stem positions (to avoid generating several positions representing the same stem)
    std::vector<pcl::PointIndices> stemid_clusters;
    OctreeEuclideanClusterExtraction<pcl::PointXYZL> oec_stem;
    oec_stem.setClusterTolerance(min_dist_between_stems_);
    oec_stem.setInputCloud(candidate_positions);
    oec_stem.extract(stemid_clusters);

    for (auto& stemid_cluster : stemid_clusters) {
        if (stemid_cluster.indices.size() == 1) { // qualified stem positions
            const auto& stemid = stemid_cluster.indices[0];
            Point3D position;
            position.x = candidate_positions->points[stemid].x;
            position.y = candidate_positions->points[stemid].y;
            position.z = candidate_positions->points[stemid].z;
            tree_positions->push_back(position);
        }
        else { // unqualified stem positions (that are too close to each other)
            Cloud3D::Ptr merged_cluster(new Cloud3D);
            for (auto& stemid : stemid_cluster.indices) {
                const auto& cluster_id = candidate_positions->points[stemid].label;
                for (auto& ptid : ptid_clusters[cluster_id].indices)
                    merged_cluster->push_back((*cloud_stems_)[ptid]);
            }
            // re-estimate a stem position from the merged cluster
            pcl::ModelCoefficients::Ptr coefficients_cylinder(new pcl::ModelCoefficients);
            if (fitOneStem(merged_cluster, 0.1, coefficients_cylinder)) {
                // transform the point-slope line coefficients to two 3D points
                double pt_bottom[3], pt_top[3];
                transformCylinderAxisToTwoPoints(
                    *coefficients_cylinder, z_bottom, z_top, pt_bottom, pt_top);
                vtkSmartPointer<vtkPoints> intersections = vtkSmartPointer<vtkPoints>::New();
                obbtree->IntersectWithLine(pt_bottom, pt_top, intersections, nullptr);
                if (intersections->GetNumberOfPoints() == 1) {
                    double intersection[3];
                    intersections->GetPoint(0, intersection);
                    Point3D position;
                    position.x = static_cast<float>(intersection[0]);
                    position.y = static_cast<float>(intersection[1]);
                    position.z = static_cast<float>(intersection[2]);
                    tree_positions->push_back(position);
                }
            }
        }
    }

    std::cout << "extract position success, and there has " << tree_positions->size() << " positions" << std::endl;

}


// NingWang: extract every single tree
void Mapping::extractAllTree(const std::string& view_output_dir, const std::string& base_name) {
    // Check if the directory exists, if not create it
    if (!std::filesystem::exists(view_output_dir)) {
        std::filesystem::create_directories(view_output_dir);
    }

    // Set cylinder radius and height range
    double radius = 1.0;//FIXME
    double height_min = -2.0;
    double height_max = 5.0;

    pcl::KdTreeFLANN<Point3D> kdtree;
    kdtree.setInputCloud(cloud_input_);

    // Create a vector of local point clouds (one per thread)
    std::vector<Cloud3D::Ptr> local_trees(tree_positions->size());
    for (auto& cloud : local_trees) {
        cloud.reset(new Cloud3D);
    }

    std::vector<std::mutex> tree_mutexes(tree_positions->size());  // every tree has its mutex

#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for (int tree_index = 0; tree_index < tree_positions->size(); ++tree_index) {
        const auto& tree_position = tree_positions->points[tree_index];
        std::vector<int> point_indices;
        std::vector<float> point_distances;

        kdtree.radiusSearch(tree_position, height_max + 1.0, point_indices, point_distances);

#if defined(_OPENMP)
#pragma omp parallel for
#endif
        for (int i = 0; i < point_indices.size(); ++i) {
            int idx = point_indices[i];
            if (idx < 0 || idx >= cloud_input_->points.size()) continue;

            const auto& pt = cloud_input_->points[idx];
            double dx = pt.x - tree_position.x;
            double dy = pt.y - tree_position.y;
            double dz = pt.z - tree_position.z;
            double dist = sqrt(dx * dx + dy * dy);

            if (dist <= radius && dz >= height_min && dz <= height_max) {
                std::lock_guard<std::mutex> lock(tree_mutexes[tree_index]);
                local_trees[tree_index]->push_back(pt);
            }
        }
    }


    // Sort by number of points
    std::vector<std::pair<int, size_t>> tree_sizes;
    for (int i = 0; i < local_trees.size(); ++i) {
        tree_sizes.emplace_back(i, local_trees[i]->size());
    }

    std::sort(tree_sizes.begin(), tree_sizes.end(),
        [](const std::pair<int, size_t>& a, const std::pair<int, size_t>& b) {
            return a.second > b.second;
        });
   
    std::cout << "this step is no error" << std::endl;

    // Save each tree
#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for (int rank = 0; rank < tree_sizes.size(); ++rank) {
        int sorted_index = tree_sizes[rank].first;
        if (!local_trees[sorted_index]->empty()) {
            std::ostringstream ss;
            ss << view_output_dir << "/Tree_" << (rank + 1) << "_" << base_name << ".pcd";
            pcl::io::savePCDFileASCII(ss.str(), *local_trees[sorted_index]);
        }
    }

    std::cout << "Extract trees for " << base_name << " success!" << std::endl;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Mapping::extract(const std::string& view_output_dir, const std::string& base_name) {
    if (!cloud_input_) {
        PCL_WARN("[Mapping::extract] No input dataset given!\n");
        tree_positions->width = tree_positions->height = 0;
        tree_positions->points.clear();
        return;
    }

    // initialize
    tree_positions->clear();
    indices_understory_->indices.clear();
    cloud_stems_->clear();
    mesh_ground_->Initialize();

    Indices::Ptr indices_understory(new Indices);
    extractUnderstory();
    extractStemPoints();
    extractTreePositions();
    extractAllTree(view_output_dir, base_name);
}

