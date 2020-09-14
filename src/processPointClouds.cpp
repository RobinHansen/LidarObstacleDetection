// PCL lib Functions for processing point clouds 

#include "processPointClouds.h"


//constructor:
template<typename PointT>
ProcessPointClouds<PointT>::ProcessPointClouds() {}


//de-constructor:
template<typename PointT>
ProcessPointClouds<PointT>::~ProcessPointClouds() {}


template<typename PointT>
void ProcessPointClouds<PointT>::numPoints(typename pcl::PointCloud<PointT>::Ptr cloud)
{
    std::cout << cloud->points.size() << std::endl;
}


template<typename PointT>
typename pcl::PointCloud<PointT>::Ptr ProcessPointClouds<PointT>::FilterCloud(typename pcl::PointCloud<PointT>::Ptr cloud, float filterRes, Eigen::Vector4f minPoint, Eigen::Vector4f maxPoint)
{

    // Time segmentation process
    auto startTime = std::chrono::steady_clock::now();

    // 1. downsampling cloud into "Voxels" (cubes) of size 'filterRes' 
    typename pcl::PointCloud<PointT>::Ptr cloud_filtered (new pcl::PointCloud<PointT>());
    typename pcl::VoxelGrid<PointT> vg;
    vg.setInputCloud(cloud);
    vg.setLeafSize(filterRes, filterRes, filterRes);
    vg.filter(*cloud_filtered);

    // 2. filter out points outside of box (minPoint, maxPoint) - so too far away
    typename pcl::PointCloud<PointT>::Ptr cloud_cropped (new pcl::PointCloud<PointT>());
    pcl::CropBox<PointT> cropbox1(true);
    cropbox1.setMin(minPoint);
    cropbox1.setMax(maxPoint);
    cropbox1.setInputCloud(cloud_filtered);
    cropbox1.filter(*cloud_cropped);

    // 3. filter out points too near to the Lidar itself (roof reflections)
    pcl::PointIndices::Ptr roof_points(new pcl::PointIndices);
    pcl::CropBox<PointT> cropbox2(true);
    cropbox2.setMin(Eigen::Vector4f(-1.5, -1.7, -1, 1));
    cropbox2.setMax(Eigen::Vector4f(2.6, 1.7, -0.4, 1));
    cropbox2.setInputCloud(cloud_cropped);
    cropbox2.filter(roof_points->indices);

    pcl::ExtractIndices<PointT> extract;
    extract.setInputCloud(cloud_cropped);
    extract.setIndices(roof_points);
    extract.setNegative(true);
    extract.filter(*cloud_cropped);

    auto endTime = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    //std::cout << "filtering took " << elapsedTime.count() << " milliseconds" << std::endl;

    return cloud_cropped;
}


template<typename PointT>
std::pair<typename pcl::PointCloud<PointT>::Ptr, typename pcl::PointCloud<PointT>::Ptr> ProcessPointClouds<PointT>::SeparateClouds(pcl::PointIndices::Ptr inliers, typename pcl::PointCloud<PointT>::Ptr cloud) 
{
  // TODO: Create two new point clouds, one cloud with obstacles and other with segmented plane
    typename pcl::PointCloud<PointT>::Ptr plane(new pcl::PointCloud<PointT>()),
                                          obstacles(new pcl::PointCloud<PointT>());

    pcl::ExtractIndices<PointT> extract;
    extract.setInputCloud(cloud);
    extract.setIndices(inliers);
    extract.setNegative(false);
    extract.filter(*plane);

    extract.setNegative (true);
    extract.filter (*obstacles);

    return std::pair<typename pcl::PointCloud<PointT>::Ptr, typename pcl::PointCloud<PointT>::Ptr>(plane, obstacles);
}


template<typename PointT>
std::pair<typename pcl::PointCloud<PointT>::Ptr, typename pcl::PointCloud<PointT>::Ptr> ProcessPointClouds<PointT>::SegmentPlane(typename pcl::PointCloud<PointT>::Ptr cloud, int maxIterations, float distanceThreshold)
{
    // Time segmentation process
    auto startTime = std::chrono::steady_clock::now();

    pcl::ModelCoefficients::Ptr coefficients (new pcl::ModelCoefficients);  // coeffs are stored here
    pcl::PointIndices::Ptr inliers (new pcl::PointIndices);
    
    pcl::SACSegmentation<PointT> seg;  // segmentation object
    seg.setOptimizeCoefficients (true);       // Optional

    seg.setModelType (pcl::SACMODEL_PLANE);
    seg.setMethodType (pcl::SAC_RANSAC);
    seg.setDistanceThreshold (distanceThreshold);          // not sure... is this the distance in 3D.

    seg.setMaxIterations(maxIterations);
    seg.setInputCloud (cloud);
    seg.segment (*inliers, *coefficients);              // do work!

    auto endTime = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "plane segmentation took " << elapsedTime.count() << " milliseconds" << std::endl;
    //std::cout << "coeffs " << *coefficients << std::endl;

    std::pair<typename pcl::PointCloud<PointT>::Ptr, typename pcl::PointCloud<PointT>::Ptr> segResult = SeparateClouds(inliers,cloud);
    return segResult;
}


template<typename PointT>
std::vector<typename pcl::PointCloud<PointT>::Ptr> ProcessPointClouds<PointT>::Clustering(typename pcl::PointCloud<PointT>::Ptr cloud, float clusterTolerance, int minSize, int maxSize)
{

    // Time clustering process
    auto startTime = std::chrono::steady_clock::now();

    std::vector<typename pcl::PointCloud<PointT>::Ptr> clusters;


    // Creating the KdTree object for the search method of the extraction
    typename pcl::search::KdTree<PointT>::Ptr tree (new pcl::search::KdTree<PointT>);
    tree->setInputCloud (cloud);

    std::vector<pcl::PointIndices> cluster_indices;
    pcl::EuclideanClusterExtraction<PointT> ec;
    ec.setClusterTolerance(clusterTolerance);
    ec.setMinClusterSize(minSize);
    ec.setMaxClusterSize(maxSize);
    ec.setSearchMethod (tree);
    ec.setInputCloud(cloud);
    ec.extract (cluster_indices);

    // TODO:: Fill in the function to perform euclidean clustering to group detected obstacles

    auto endTime = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "clustering took " << elapsedTime.count() << " milliseconds and found " << clusters.size() << " clusters" << std::endl;

    for (const auto &ci : cluster_indices) {
        typename pcl::PointCloud<PointT>::Ptr cluster(new pcl::PointCloud<PointT>);
//        for (std::vector<int>::const_iterator pit = ci.indices.begin (); pit != ci.indices.end (); ++pit)
        for (const auto &e : ci.indices)
            cluster->push_back ((*cloud)[e]); 
        cluster->width = cluster->size ();
        cluster->height = 1;
        cluster->is_dense = true;
        clusters.push_back(cluster);
    }
        
    return clusters;
}


template<typename PointT>
Box ProcessPointClouds<PointT>::BoundingBox(typename pcl::PointCloud<PointT>::Ptr cluster)
{

    // Find bounding box for one of the clusters
    PointT minPoint, maxPoint;
    pcl::getMinMax3D(*cluster, minPoint, maxPoint);

    Box box;
    box.x_min = minPoint.x;
    box.y_min = minPoint.y;
    box.z_min = minPoint.z;
    box.x_max = maxPoint.x;
    box.y_max = maxPoint.y;
    box.z_max = maxPoint.z;

    return box;
}


template<typename PointT>
void ProcessPointClouds<PointT>::savePcd(typename pcl::PointCloud<PointT>::Ptr cloud, std::string file)
{
    pcl::io::savePCDFileASCII (file, *cloud);
    std::cerr << "Saved " << cloud->points.size () << " data points to "+file << std::endl;
}


template<typename PointT>
typename pcl::PointCloud<PointT>::Ptr ProcessPointClouds<PointT>::loadPcd(std::string file)
{

    typename pcl::PointCloud<PointT>::Ptr cloud (new pcl::PointCloud<PointT>);

    if (pcl::io::loadPCDFile<PointT> (file, *cloud) == -1) //* load the file
    {
        PCL_ERROR ("Couldn't read file \n");
    }
    std::cerr << "Loaded " << cloud->points.size () << " data points from "+file << std::endl;

    return cloud;
}


template<typename PointT>
std::vector<boost::filesystem::path> ProcessPointClouds<PointT>::streamPcd(std::string dataPath)
{

    std::vector<boost::filesystem::path> paths(boost::filesystem::directory_iterator{dataPath}, boost::filesystem::directory_iterator{});

    // sort files in accending order so playback is chronological
    sort(paths.begin(), paths.end());

    return paths;

}