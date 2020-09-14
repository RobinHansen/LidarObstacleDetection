/* \author Aaron Brown */
// Quiz on implementing simple RANSAC line fitting

#include "../../render/render.h"
#include <unordered_set>
#include "../../processPointClouds.h"
// using templates for processPointClouds so also include .cpp to help linker
#include "../../processPointClouds.cpp"

pcl::PointCloud<pcl::PointXYZ>::Ptr CreateData()
{
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>());
  	// Add inliers
  	float scatter = 0.6;
  	for(int i = -5; i < 5; i++)
  	{
  		double rx = 2*(((double) rand() / (RAND_MAX))-0.5);
  		double ry = 2*(((double) rand() / (RAND_MAX))-0.5);
  		cloud->points.push_back(pcl::PointXYZ(i+scatter*rx, i+scatter*ry, 0));
  	}

	int numOutliers = 10;
  	while(numOutliers-- > 0)
  	{
  		double rx = 2*(((double) rand() / (RAND_MAX))-0.5);
  		double ry = 2*(((double) rand() / (RAND_MAX))-0.5);

  		cloud->points.push_back(pcl::PointXYZ(5*rx, 5*ry, 0));

  	}
  	cloud->width = cloud->points.size();
  	cloud->height = 1;

  	return cloud;

}

pcl::PointCloud<pcl::PointXYZ>::Ptr CreateData3D()
{
	ProcessPointClouds<pcl::PointXYZ> pointProcessor;
	return pointProcessor.loadPcd("../../../sensors/data/pcd/simpleHighway.pcd");
}


pcl::visualization::PCLVisualizer::Ptr initScene()
{
	pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer ("2D Viewer"));
	viewer->setBackgroundColor (0, 0, 0);
  	viewer->initCameraParameters();
  	viewer->setCameraPosition(0, 0, 15, 0, 1, 0);
  	viewer->addCoordinateSystem (1.0);
  	return viewer;
}

/** return random interger in range [r0;r1[ */
int randi(int r0, int r1) {
	return r0 + rand() % (r1-r0);
}

std::unordered_set<int> Ransac(pcl::PointCloud<pcl::PointXYZ>::Ptr cloudptr, int maxIterations, float distanceTol)
{
	std::unordered_set<int> inliersResult;
	srand(time(NULL));

	const pcl::PointCloud<pcl::PointXYZ> &cloud = *cloudptr;
	int half = cloud.size()/2;

	for (int i=0; i < maxIterations; ++i) {
		int index0 = randi(0,cloud.size());
		int index1 = index0 > half ? randi(0,half) : randi(half,cloud.size());

		const pcl::PointXYZ &p0 = cloud[index0],
		                    &p1 = cloud[index1];

		float A = p0.y - p1.y,
			  B = p1.x - p0.x,
			  C = p0.x*p1.y - p1.x*p0.y;
		
		std::unordered_set<int> inliersResultCandidate;
		double sqrt_AA_BB = sqrt(A*A + B*B);
		for (int i = 0; i < cloud.size(); ++i) 
			if (i != index0 && i != index1) {
				double d = fabs(A*cloud[i].x + B*cloud[i].y + C) / sqrt_AA_BB;
				if (d <= distanceTol)
					inliersResultCandidate.insert(i);
			}
		if (inliersResultCandidate.size() > inliersResult.size()) {
			std::cout << "new candidate with " << inliersResultCandidate.size() << " inliers; was " << inliersResult.size() << std::endl;
			inliersResult = inliersResultCandidate;
			inliersResult.insert(index0);
			inliersResult.insert(index1);
		}
	}
	
	return inliersResult;
}

std::unordered_set<int> Ransac3D(pcl::PointCloud<pcl::PointXYZ>::Ptr cloudptr, int maxIterations, float distance)
{
	srand(time(NULL));

	const pcl::PointCloud<pcl::PointXYZ> &cloud = *cloudptr;

	int n = 0;
	double best_A, best_B, best_C, best_D;

	while(maxIterations-- > 0) {
		int index1, index2, index3;
		index1 = arc4random_uniform(cloud.size());
		do index2 = arc4random_uniform(cloud.size()); while (index2 == index1);
		do index3 = arc4random_uniform(cloud.size()); while (index3 == index2 || index3 == index1);

		const pcl::PointXYZ &p1 = cloud[index1],
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

int main ()
{

	// Create viewer
	pcl::visualization::PCLVisualizer::Ptr viewer = initScene();

	// Create data
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = CreateData3D();
	

	// TODO: Change the max iteration and distance tolerance arguments for Ransac function
    auto startTime = std::chrono::steady_clock::now();
	std::unordered_set<int> inliers = Ransac3D(cloud, 500, 0.2);
    auto endTime = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "plane segmentation took " << elapsedTime.count() << " milliseconds" << std::endl;

	pcl::PointCloud<pcl::PointXYZ>::Ptr  cloudInliers(new pcl::PointCloud<pcl::PointXYZ>());
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloudOutliers(new pcl::PointCloud<pcl::PointXYZ>());

	for(int index = 0; index < cloud->points.size(); index++)
	{
		pcl::PointXYZ point = cloud->points[index];
		if(inliers.count(index))
			cloudInliers->points.push_back(point);
		else
			cloudOutliers->points.push_back(point);
	}


	// Render 2D point cloud with inliers and outliers
	if(inliers.size())
	{
		renderPointCloud(viewer,cloudInliers,"inliers",Color(0,1,0));
  		renderPointCloud(viewer,cloudOutliers,"outliers",Color(1,0,0));
	}
  	else
  	{
  		renderPointCloud(viewer,cloud,"data");
  	}
	
  	while (!viewer->wasStopped ())
  	{
  	  viewer->spinOnce ();
  	}
  	
}
