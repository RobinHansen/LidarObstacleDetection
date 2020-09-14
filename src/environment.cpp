#include "sensors/lidar.h"
#include "render/render.h"
#include "processPointClouds.h"
// using templates for processPointClouds so also include .cpp to help linker
#include "processPointClouds.cpp"
#include <unordered_set>

struct KdTree
{
    struct Node
    {
        pcl::PointXYZI point;
        int id;
        Node* left = nullptr;
        Node* right = nullptr;

        Node() {}
        Node(pcl::PointXYZI &point, int id) : point(point), id(id) {}
    };
	
    Node* root;
    Node *nodes = nullptr;  // array of Nodes, consecutive in memory, allocated with one new, disposed of with one delete
    unsigned next_node_idx = 0;

    KdTree() : root(NULL) {}

	// initialize from a PCD
    KdTree(pcl::PointCloud<pcl::PointXYZI>::Ptr pcl) : root(NULL) {
        nodes = new Node[pcl->points.size()];
        for (int i = 0; i < pcl->points.size(); ++i)
            insert(pcl->points[i], i);
    }

    ~KdTree() {
        delete [] nodes;
    }

	void insert(pcl::PointXYZI &point, int id)
	{	
		//  obvious recursive algo is tail-recursive, hence an iterative variant is straight forward
		Node **n = &root;
		for( int dim = 0; *n != nullptr; dim = (dim+1)%3)
			n = &(point.data[dim] <= (*n)->point.data[dim] ? (*n)->left : (*n)->right);
		*n = new (&nodes[next_node_idx]) Node(point, id);		
        next_node_idx++;
	}

	void searchHelper(std::vector<int> &ids, const pcl::PointXYZI &target, float distanceTol, const Node *n, int depth) const {
		if (n != nullptr) {            
			float d_x = target.x - n->point.x,
                  d_y = target.y - n->point.y,
                  d_z = target.z - n->point.z;

			if ( sqrt(d_x*d_x + d_y*d_y + d_z*d_z) <= distanceTol )
				ids.push_back(n->id);

            // Note: target.x is an alias for target.data[0]  (and likewise for y and z)
            // (https://pcl.readthedocs.io/projects/tutorials/en/latest/adding_custom_ptype.html#adding-custom-ptype)
			int dim = depth % 3;
			if (n->point.data[dim] >= target.data[dim] - distanceTol)
				searchHelper(ids, target, distanceTol, n->left, depth+1);
			if (n->point.data[dim] <= target.data[dim] + distanceTol)
				searchHelper(ids, target, distanceTol, n->right, depth+1);
		}
	}

	// return a list of point ids in the tree that are within distance of target
	std::vector<int> search(const pcl::PointXYZI &target, float distanceTol) const
	{
		std::vector<int> ids;
		searchHelper(ids, target, distanceTol, root, 0);
		return ids;
	}
};

template<typename PointT>
std::unordered_set<int> Ransac_plane_inliers(typename pcl::PointCloud<PointT>::Ptr cloudptr, int maxIterations, float distance)
{
	srand(time(NULL));

	const pcl::PointCloud<PointT> &cloud = *cloudptr;

	int n = 0;
	double best_A, best_B, best_C, best_D;

	while(maxIterations-- > 0) {
		int index1, index2, index3;
		index1 = arc4random_uniform(cloud.size());
		do index2 = arc4random_uniform(cloud.size()); while (index2 == index1);
		do index3 = arc4random_uniform(cloud.size()); while (index3 == index2 || index3 == index1);

		const pcl::PointXYZI &p1 = cloud[index1],
		                     &p2 = cloud[index2],
		                     &p3 = cloud[index3];

		double x1 = p1.x, y1 = p1.y, z1 = p1.z,
		       x2 = p2.x, y2 = p2.y, z2 = p2.z,
			   x3 = p3.x, y3 = p3.y, z3 = p3.z;

		double i = (y2-y1)*(z3-z1) - (z2-z1)*(y3-y1),
               j = (z2-z1)*(x3-x1) - (x2-x1)*(z3-z1),
			   k = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1);

		double A = i,
		       B = j,
			   C = k,
			   D = -(i*x1 + j*y1 + k*z1);

		double sqrt_AA_BB_CC = sqrt(A*A + B*B + C*C);
		int inliers_count = 0;
		double distance_ = distance * sqrt_AA_BB_CC;

		for (const auto &pt : cloud)
			if (fabs(A*pt.x + B*pt.y + C * pt.z + D) <= distance_)
				++inliers_count;

        // store plane coefficients if this iteration provided more inliers.
		if (inliers_count > n) {
			best_A = A;
			best_B = B;
			best_C = C;
			best_D = D;
			n = inliers_count;
		}

	}

	std::unordered_set<int> inliersResult(n);
	
	double sqrt_AA_BB_CC = sqrt(best_A*best_A + best_B*best_B + best_C*best_C);
	double distance_ = distance * sqrt_AA_BB_CC; // take sqrt and division out of loop
	for (int i = 0; i < cloud.size(); ++i) {
		double d_ = fabs(best_A*cloud[i].x + best_B*cloud[i].y + best_C * cloud[i].z + best_D);
		if (d_ < distance_)
			inliersResult.insert(i);
	}
	
	return inliersResult;
}

class EuclideanCluster {
    // C++ flavour of a "closure" in order to avoid recursively passing invariant arguments
    // so in turn a minor runtime optimization
    // intended use:  std::vector<std::vector<int>> clusters = EuclideanCluster(params).find();
    
    const pcl::PointCloud<pcl::PointXYZI>::Ptr pcd;
    const KdTree &tree;
    float distanceTol;
	std::vector<bool> processed;

    void proximity(std::vector<int> &cluster, int i) {
        processed[i] = true;
        cluster.push_back(i);
        for (auto j : tree.search(pcd->points[i], distanceTol))
            if (!processed[j])
                proximity(cluster, j);
    }

public:
    EuclideanCluster(const pcl::PointCloud<pcl::PointXYZI>::Ptr pcd, const KdTree &tree, float distanceTol) :
        pcd(pcd), tree(tree), distanceTol(distanceTol), processed(pcd->points.size(), false) {}
    
    std::vector<std::vector<int>> find() {
    	std::vector<std::vector<int>> clusters;
        for (int i = 0; i < pcd->points.size(); ++i) {
            if (!processed[i]) {
                std::vector<int> cluster;
                proximity(cluster, i);
                clusters.push_back(cluster);
            }
        }
        return clusters;
    }
};

void cityBlock(pcl::visualization::PCLVisualizer::Ptr& viewer, ProcessPointClouds<pcl::PointXYZI> &pointProcessor,
                const pcl::PointCloud<pcl::PointXYZI>::Ptr& inputCloud)
{
    auto startTime = std::chrono::steady_clock::now();

    //std::cout << "point before filtering: " << inputCloud->size() << std::endl;
    typename pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_filtered = pointProcessor.FilterCloud(inputCloud, 0.1f, Eigen::Vector4f(-10,-5,-2,1), Eigen::Vector4f(20,7,5,1));
    //std::cout << "point after filtering: " << cloud_filtered->size() << std::endl;

    std::unordered_set<int> inliers = Ransac_plane_inliers<pcl::PointXYZI>(cloud_filtered, 50, 0.2);
  
    pcl::PointCloud<pcl::PointXYZI>::Ptr plane(new pcl::PointCloud<pcl::PointXYZI>);
    pcl::PointCloud<pcl::PointXYZI>::Ptr obstacles(new pcl::PointCloud<pcl::PointXYZI>);

	for(int index = 0; index < cloud_filtered->points.size(); ++index)
		if(inliers.count(index))
			plane->points.push_back(cloud_filtered->points[index]);
		else
			obstacles->points.push_back(cloud_filtered->points[index]);
   
    renderPointCloud(viewer, plane, "plane", Color(1,1,1));

    KdTree kdtree(obstacles);

    std::vector<std::vector<int>> clusters = EuclideanCluster(obstacles, kdtree, 0.25).find();

  	// Render clusters
  	int clusterId = 0;
	std::vector<Color> colors = {Color(1,0,0), Color(0,1,0), Color(0,0,1)};
  	for(std::vector<int> cluster : clusters)
  	{
  		pcl::PointCloud<pcl::PointXYZI>::Ptr one_obstacle(new pcl::PointCloud<pcl::PointXYZI>());
  		for(int i: cluster)
  			one_obstacle->points.push_back(obstacles->points[i]);
  		renderPointCloud(viewer, one_obstacle, "cluster"+std::to_string(clusterId),colors[clusterId % 3]);
        renderBox(viewer, pointProcessor.BoundingBox(one_obstacle), clusterId);
 		++clusterId;
  	}

    Box egocar;
    egocar.x_min = -2;
    egocar.y_min = -1;
    egocar.z_min = -1.5;
    egocar.x_max = 2;
    egocar.y_max = 1;
    egocar.z_max = 0;
    int egocar_id = 999999;
    renderBox(viewer, egocar, egocar_id, Color(0,0,1));

    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime);
    std::cout << "frame render time: " << elapsedTime.count() << " milliseconds" << std::endl;

}

//setAngle: SWITCH CAMERA ANGLE {XY, TopDown, Side, FPS}
void initCamera(CameraAngle setAngle, pcl::visualization::PCLVisualizer::Ptr& viewer)
{

    viewer->setBackgroundColor (0, 0, 0);
    
    // set camera position and angle
    viewer->initCameraParameters();
    // distance away in meters
    int distance = 16;
    
    switch(setAngle)
    {
        case XY : viewer->setCameraPosition(-distance, -distance, distance, 1, 1, 0); break;
        case TopDown : viewer->setCameraPosition(0, 0, distance, 1, 0, 1); break;
        case Side : viewer->setCameraPosition(0, -distance, 0, 0, 0, 1); break;
        case FPS : viewer->setCameraPosition(-10, 0, 0, 0, 0, 1);
    }

    if(setAngle!=FPS)
        viewer->addCoordinateSystem (1.0);
}


int main (int argc, char** argv)
{
    std::cout << "starting enviroment" << std::endl;
    pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer ("3D Viewer"));
    CameraAngle setAngle = XY;
    initCamera(setAngle, viewer);

    ProcessPointClouds<pcl::PointXYZI> pointProcessor;
    std::vector<boost::filesystem::path> stream = pointProcessor.streamPcd("../src/sensors/data/pcd/data_1");

    while (!viewer->wasStopped ())
        for ( auto const &path : stream)
        {
            viewer->removeAllPointClouds();
            viewer->removeAllShapes();
            cityBlock(viewer, pointProcessor, pointProcessor.loadPcd(path.string()));
            viewer->spinOnce ();
        } 
}