/* \author Aaron Brown */
// Quiz on implementing kd tree

#include "../../render/render.h"


// Structure to represent node of kd tree
struct Node
{
	std::vector<float> point;
	int id;
	Node* left = nullptr;
	Node* right = nullptr;

	Node(std::vector<float> point, int id) : point(point), id(id) {}
};

struct KdTree
{
	Node* root;

	KdTree()
	: root(NULL)
	{}

	void _insert_y(Node *&n, const std::vector<float> &point, int id) {
		if (n == NULL)
			n = new Node(point, id);
		else
			_insert_x(point[1] <= n->point[1] ? n->left : n->right, point, id);
	}

	void _insert_x(Node *&n, const std::vector<float> &point, int id) {
		if (n == NULL)
			n = new Node(point, id);
		else
			_insert_y(point[0] <= n->point[0] ? n->left : n->right, point, id);
	}

	void rec_insert(std::vector<float> point, int id)
	{
		_insert_x(root, point, id);
	}

	void insert(std::vector<float> point, int id)
	{	
		/*  try this, most simple, most efficient */
		Node **n = &root;
		for( int dim = 0; *n != nullptr; dim = (dim+1)%2)
			n = &(point[dim] <= (*n)->point[dim] ? (*n)->left : (*n)->right);
		*n = new Node(point, id);
		
	}

	void searchHelper(std::vector<int> &ids, const std::vector<float> &target, float distanceTol, const Node *n, int depth) {
		if (n != nullptr) {
			//std::cout << "considering " << n->id << " (" << n->point[0] << "," << n->point[1] << ") " << depth << std::endl;
			float d_x = target[0] - n->point[0], d_y = target[1] - n->point[1];
			if ( sqrt(d_x*d_x + d_y*d_y) <= distanceTol )
				ids.push_back(n->id);

			int dim = depth % 2;
			if (n->point[dim] >= target[dim] - distanceTol)
				searchHelper(ids, target, distanceTol, n->left, depth+1);
			if (n->point[dim] <= target[dim] + distanceTol)
				searchHelper(ids, target, distanceTol, n->right, depth+1);
		}
	}

	// return a list of point ids in the tree that are within distance of target
	std::vector<int> search(std::vector<float> target, float distanceTol)
	{
		std::vector<int> ids;
		searchHelper(ids, target, distanceTol, root, 0);
		return ids;
	}
	
	void dotOutputHelper(Node *n, int depth) {
		const char *color = depth % 2 == 0 ? "green" : "blue";

		std::cout << n->id << " [label=\"" << n->id << "\\n(" << n->point[0] << "," << n->point[1] << ")\"]\n";

		if (n->left)
			std::cout << n->id << "->" << n->left->id << " [label=ls color=" << color << "];" << std::endl;
		else if (n->right)
			std::cout << "N" << n->id << " [shape=point];\n"
					  << n->id << "->" << "N" << n->id << " [color=" << color << "];" << std::endl;
		if (n->right)
			std::cout << n->id << "->" << n->right->id << " [label=gt color=" << color << "];" << std::endl;
		else if (n->left)
			std::cout << "N" << n->id << " [shape=point];\n"
					  << n->id << "->" << "N" << n->id << " [color=" << color << "];" << std::endl;

		if (n->left)
			dotOutputHelper(n->left, depth+1);
		if (n->right)
			dotOutputHelper(n->right, depth+1);

	}

	void dotOutput() {
		std::cout << "digraph BT {\n";
		if (root)
			dotOutputHelper(root, 0);
		std::cout << "}" << std::endl;
	}

};
