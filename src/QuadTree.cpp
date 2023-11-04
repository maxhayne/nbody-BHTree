#include "QuadTree.h"
#include <cstddef>
#include <cmath>
#include <iostream>
#include <omp.h>
using namespace std;

bool QuadTree::insert(Body* b) {
	// If not within this QT's area, return false
	if (b == NULL || !withinBounds(b)) {
		return false;
	}

	// Updating values for calculating COM
	weightedX += b->x * b->mass;
	weightedY += b->y * b->mass;
	totalMass += b->mass;

	// If the QuadTree isn't already storing a body, store the body
	if (body == NULL && subdivided == false) {
		body = b;
		return true;
	}
	
	// If the QuadTree isn't subdivided, subdivide and add body to a subdivision
	if (subdivided == false) {
		subdivide();
		if (!northwest->insert(body)) {
			if (!southwest->insert(body)) {
				if (!southeast->insert(body)) {
					northeast->insert(body);
				}
			}
		}
		body = NULL; // Set pointer to NULL
	}

	// Insert new body into a subdivision
	if (!northwest->insert(b)) {
		if (!southwest->insert(b)) {
			if (!southeast->insert(b)) {
				northeast->insert(b);
			}
		}
	}

	return true;
}

bool QuadTree::withinBounds(Body* b) {
	return(b->x >= lowerLeftCorner->x && \
		b->x <= lowerLeftCorner->x + sideLength && \
		b->y >= lowerLeftCorner->y && \
		b->y <= lowerLeftCorner->y + sideLength);
}

void QuadTree::subdivide() {
	southwest = new QuadTree(new Point(lowerLeftCorner), sideLength/2);
	northwest = new QuadTree(new Point(lowerLeftCorner->x, lowerLeftCorner->y + sideLength/2), sideLength/2);
	southeast = new QuadTree(new Point(lowerLeftCorner->x + sideLength/2, lowerLeftCorner->y), sideLength/2);
	northeast = new QuadTree(new Point(lowerLeftCorner->x + sideLength/2, lowerLeftCorner->y + sideLength/2), sideLength/2);
	southwest->parent = this;
	northwest->parent = this;
	southeast->parent = this;
	northeast->parent = this;
	subdivided = true;
}

void QuadTree::calculateCOM() {
	com.x = weightedX/totalMass;
	com.y = weightedY/totalMass;
	com.mass = totalMass;
	return;
}

int QuadTree::longestPath(int x) {
	if (northwest == NULL && \
		southwest == NULL && \
		southeast == NULL && \
		northeast == NULL) {
		return(x);
	}

	int nw = -1, sw = -1, se = -1, ne = -1, longestPath = -1;
	if (northwest != NULL) {
		nw = northwest->longestPath(x+1);
		longestPath = nw;
	}
	if (southwest != NULL) {
		sw = southwest->longestPath(x+1);
		if (sw > longestPath) {
			longestPath = sw;
		}
	}
	if (southeast != NULL) {
		se = southeast->longestPath(x+1);
		if (se > longestPath) {
			longestPath = se;
		}
	}
	if (northeast != NULL) {
		ne = northeast->longestPath(x+1);
		if (ne > longestPath) {
			longestPath = ne;
		}
	}

	return(longestPath);
}

DoubleDouble QuadTree::calculateNetForce(Body *b) {
	DoubleDouble f = calculateAcceleration(b);
	f.x = f.x * b->mass;
	f.y = f.y * b->mass;
	return f;
}

DoubleDouble QuadTree::calculateAcceleration(Body *b) {
	DoubleDouble a = accelerationHelper(b);
	a.x = a.x * G;
	a.y = a.y * G;
	return a;
}

DoubleDouble QuadTree::accelerationHelper(Body *b) {
	DoubleDouble total(0,0);

	if (!subdivided && body == NULL) {
		return total;
	}

	// If the box contains the body we are calculating for, we must recurse
	// Also, if s/d is greater than theta, we must recurse

	// If the node we're in isn't subdivided, can just calculate the force straight out
	if (subdivided == false && body != NULL) {
		if (body != b) {
			calculateCOM();
			double distance = sqrt(pow(com.x-b->x, 2) + pow(com.y-b->y, 2));
			double accel = com.mass/(distance*distance);
			double dVec[2] = {com.x-b->x, com.y-b->y};
			total.x = accel * (dVec[0]/distance); 
			total.y = accel * (dVec[1]/distance);
			return total;
		} else {
			return total;
		}
	}

	calculateCOM(); // Calculate COM
	double distance = sqrt(pow(b->x-com.x, 2) + pow(b->y-com.y, 2));
	if (withinBounds(b) || sideLength/distance > theta) { // Recurse, theta is 0.5
		DoubleDouble nw(0,0);
		DoubleDouble sw(0,0);
		DoubleDouble se(0,0);
		DoubleDouble ne(0,0);
		if (northwest != NULL) {
			nw = northwest->accelerationHelper(b);
			total.x += nw.x;
			total.y += nw.y;
		}
		if (southwest != NULL) {
			sw = southwest->accelerationHelper(b);
			total.x += sw.x;
			total.y += sw.y;
		}
		if (southeast != NULL) {
			se = southeast->accelerationHelper(b);
			total.x += se.x;
			total.y += se.y;
		}
		if (northeast != NULL) {
			ne = northeast->accelerationHelper(b);
			total.x += ne.x;
			total.y += ne.y;
		}
	} else {
		double accel = com.mass/(distance*distance);
		double dVec[2] = {com.x-b->x, com.y-b->y};
		total.x = accel * (dVec[0]/distance); 
		total.y = accel * (dVec[1]/distance);
	}	
	return total;
}

int QuadTree::countNullQuadTrees() {
	int count = 0;
	
	if (northwest == NULL) {
		count++;
	} else {
		count += northwest->countNullQuadTrees();
	}

	if (southwest == NULL) {
		count++;
	} else {
		count += southwest->countNullQuadTrees();
	}

	if (southeast == NULL) {
		count++;
	} else {
		count += southeast->countNullQuadTrees();
	}

	if (northeast == NULL) {
		count++;
	} else {
		count += northeast->countNullQuadTrees();
	}

	return count;
}

int QuadTree::countQuadTrees() {
	int count = 0;
	
	if (northwest != NULL) {
		count++;
		count += northwest->countQuadTrees();
	}

	if (southwest != NULL) {
		count++;
		count += southwest->countQuadTrees();
	}

	if (southeast != NULL) {
		count++;
		count += southeast->countQuadTrees();
	}

	if (northeast != NULL) {
		count++;
		count += northeast->countQuadTrees();
	}

	return count;
}

// Won't delete body pointers. They may be kept in a list elsewhere.
void QuadTree::deleteTree(QuadTree* tree) {
	if (tree->northwest != NULL) deleteTree(tree->northwest);
	if (tree->southwest != NULL) deleteTree(tree->southwest);
	if (tree->southeast != NULL) deleteTree(tree->southeast);
	if (tree->northeast != NULL) deleteTree(tree->northeast);

	if (tree->northwest != NULL) delete tree->northwest;
	if (tree->southwest != NULL) delete tree->southwest;
	if (tree->southeast != NULL) delete tree->southeast;
	if (tree->northeast != NULL) delete tree->northeast;
	if (tree->lowerLeftCorner != NULL) delete tree->lowerLeftCorner;
	tree->northwest = NULL;
	tree->southwest = NULL;
	tree->southeast = NULL;
	tree->northeast = NULL;
	tree->parent = NULL;
	tree->body = NULL;
	tree->lowerLeftCorner = NULL;
}
