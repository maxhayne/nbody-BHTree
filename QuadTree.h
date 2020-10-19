#ifndef __QUADTREE_H__
#define __QUADTREE_H__

#include "Body.h"
#include "Point.h"
#include "DoubleDouble.h"
#include <cstddef>

const double G = 6.674e-11;
const double pi = 3.14159265;

class QuadTree {
	public:
	Point* lowerLeftCorner;
	Body* body;
	Body com;
	double sideLength;
	bool subdivided;

	QuadTree* northwest;
	QuadTree* southwest;
	QuadTree* southeast;
	QuadTree* northeast;
	QuadTree* parent;

	double weightedX;
	double weightedY;
	double totalMass;
	double theta;

	QuadTree(Point* _lowerLeftCorner, double _sideLength) {
		lowerLeftCorner = _lowerLeftCorner;
		sideLength = _sideLength;
		northwest = NULL;
		southwest = NULL;
		southeast = NULL;
		northeast = NULL;
		parent = NULL;
		body = NULL;
		weightedX = 0;
		weightedY = 0;
		totalMass = 0;
		theta = 0.5;
		subdivided = false;
	}

	bool insert(Body*);
	bool withinBounds(Body*);
	void subdivide();
	void calculateCOM();
	int longestPath(int);
	DoubleDouble calculateNetForce(Body*);
	DoubleDouble calculateAcceleration(Body*);
	DoubleDouble accelerationHelper(Body*);
	int countNullQuadTrees();
	int countQuadTrees();
	static void deleteTree(QuadTree*);
};

#endif