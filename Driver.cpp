#include "QuadTree.h"
#include <iostream>
#include <vector>
#include <ctime>
#include <list>
#include <chrono>
#include <cmath>
#include <fstream>
#include <string>
#include <omp.h>
using namespace std;

int getTime() {
	return(std::chrono::duration_cast< std::chrono::milliseconds >(
    std::chrono::system_clock::now().time_since_epoch()).count());
}

void leapfrog_part1(Body *bodies, int numberOfBodies, double dt) {
	//#pragma omp parallel for schedule(static)
	for (int i = 0; i < numberOfBodies; i++) {
		bodies[i].leapfrog_part1(dt);
	}
}

void leapfrog_part2(Body *bodies, int numberOfBodies, double dt) {
	//#pragma omp parallel for schedule(static)
	for (int i = 0; i < numberOfBodies; i++) {
		bodies[i].leapfrog_part2(dt);
	}
}

int main() {

	// First galaxy parameters
	int g1bodies = 5000;
	double g1range = 100000;
	double g1mass = 100000000;
	double g1centralmass = 100000000000000;

	// Second galaxy parameters
	int g2bodies = 1000;
	double g2range = 50000;
	double g2mass = 100000000;
	double g2centralmass = 50000000000000;
	double g2xoffset = 70000;
	double g2yoffset = 70000;
	double rg = sqrt((g2xoffset*g2xoffset) + (g2yoffset*g2yoffset));
	double vg = sqrt(G*g1centralmass/rg);
	double g2vX = vg * (g2xoffset/rg);
	double g2vY = vg * (g2yoffset/rg);

	// Simulation Parameters
	double dt = 1000; // in seconds
	double duration = dt*60*40; // in seconds
	int steps = duration/dt;

	int numberOfBodies = g1bodies+g2bodies;
	Body *bodies = new Body[numberOfBodies];
	
	// Keep track of QuadTree bounds
	double range = g1range+g2range;
	double xMin = -range/1.3;
	double yMin = -range/3;
	double xMax = range;
	double yMax = range;

	// To use for visualization
	double xMin1 = xMin*2;
	double yMin1 = yMin*2;
	double range1 = range*2.5;

	cout << "Driver is running.\n";
	cout << "Simulating " << numberOfBodies << " bodies.\n";
	srand(time(NULL)); // Initiate RNG

	int i;
	bodies[0] = new Body(0, 0, g1centralmass, 0, 0);
	for (i = 1; i < g1bodies; i++) {
		// generate random x within range
		double x = g1range/2 - rand() % (int)g1range;
		// generate random y within range of circular radius
		double yRange = 2*sqrt(((g1range/2)*(g1range/2)) - (x*x));
		double y = yRange/2 - rand() % (int)yRange;
		// find x and y components of velocity
		double r = sqrt((x*x) + (y*y));
		double velocity = sqrt(G*g1centralmass/r);
		double xV = velocity * (x/r);
		double yV = velocity * (y/r);
		bodies[i] = new Body(x, y, g1mass, -yV, xV);
	}

	bodies[i] = new Body(g2xoffset, g2yoffset, g2centralmass, -g2vY, g2vX);
	for (i = 1; i < g2bodies; i++) {
		// generate random x within range
		double x = g2range/2 - rand() % (int)g2range;
		// generate random y within range of circular radius
		double yRange = 2*sqrt(((g2range/2)*(g2range/2)) - (x*x));
		double y = yRange/2 - rand() % (int)yRange;
		// find x and y components of velocity
		double r = sqrt((x*x) + (y*y));
		double velocity = sqrt(G*g2centralmass/r);
		double xV = velocity * (x/r);
		double yV = velocity * (y/r);
		bodies[i] = new Body(x+g2xoffset, y+g2yoffset, g2mass, -yV-g2vY, xV+g2vX);
	}


	QuadTree *root;
	int start = getTime();
	for (int iterations = 0; iterations < steps; iterations++) {

		// First leapfrog drift
		leapfrog_part1(bodies, numberOfBodies, dt);

		// Update the locations of the bodies
		for (int i = 0; i < numberOfBodies; i++) {
			xMax = xMax < bodies[i].x ? bodies[i].x : xMax;
			xMin = xMin > bodies[i].x ? bodies[i].x : xMin;
			yMax = yMax < bodies[i].y ? bodies[i].y : yMax;
			yMin = yMin > bodies[i].y ? bodies[i].y : yMin;
		}
		range = ceil(xMax - xMin);
		range = range >= ceil(yMax - yMin) ? range : ceil(yMax - yMin);
		//cout << "xMin: " << xMin << ", yMin: " << yMin << ", Range: " << range << endl;

		// Build, then add the bodies to the tree
		root = new QuadTree(new Point(xMin,yMin), range);
		for (int i = 0; i < numberOfBodies; i++) {
			if(!root->insert(&bodies[i])) {
				cout << "Could not insert point at (" << bodies[i].x << ", " << bodies[i].y << ")\n";
			}
		}

		// Calculate Net Accelerations in x and y directions, AKA "kick"
		if (numberOfBodies > 65) {
			#pragma omp parallel for schedule(guided)
			for (int i = 0; i < numberOfBodies; i++) {
				DoubleDouble a = root->calculateAcceleration(&bodies[i]);
				bodies[i].aX = a.x;
				bodies[i].aY = a.y;
			}
		} else {
			for (int i = 0; i < numberOfBodies; i++) {
				DoubleDouble a = root->calculateAcceleration(&bodies[i]);
				bodies[i].aX = a.x;
				bodies[i].aY = a.y;
			}
		}

		// Second leapfrog drift
		leapfrog_part2(bodies, numberOfBodies, dt);

		// Delete the QuadTree
		QuadTree::deleteTree(root);
		root = NULL;
	}

	int bhtime = getTime() - start;
	cout << "Barnes-Hut Runtime for " << duration << " seconds with " << dt << "s steps: " << (float)bhtime/1000.0 << " seconds.\n";

	for (int i = 0; i < numberOfBodies; i++) {
		if ((int) bodies[i].locations.size() != steps) {
			cout << "Not equal.\n";
		}
	}

	string body;
	ofstream datafile;
	datafile.open("locations.data");
	if (!datafile) {
		cout << "Could not open the data file.\n";
	} else {
		datafile << numberOfBodies << "\n";
		datafile << "2\n"; // for the number of dimensions
		datafile << dt << "\n";
		datafile << duration << "\n";
		datafile << xMin1 << "," << yMin1 << "\n";
		datafile << range1 << "\n";
		// Writing all masses
		for (int i = 0; i < numberOfBodies; i++) {
			datafile << bodies[i].mass << "\n";
		}
		// Writing all positions
		for (int i = 0; i < steps; i++) {
			for (int j = 0; j < numberOfBodies; j++) {
				datafile << bodies[j].locations[i].x << "\t" << bodies[j].locations[i].y << "\n";
			}
		}
		cout << "Finished writing to 'locations.data'.\n";
	}

	delete [] bodies;
	cout << "Done.\n";
}

// This will make a visually appealing galaxy
/*
// Body Parameters
int numberOfBodies = 5000;
double range = 100000; // units undecided, m?
double massRange = 100000000; // units undecided, kg?
double centralMass = 100000000000000;
//int xVelRange = 1; // units undecided, m/s?
//int yVelRange = 1; // units undecided, m/s?

// Keep track of QuadTree bounds
double xMin = -range/2;
double yMin = -range/2;
double xMax = range/2;
double yMax = range/2;

// To use for visualization
double xMin1 = xMin;
double yMin1 = yMin;
double range1 = range;

// Simulation Parameters
double dt = 1000; // in seconds
double duration = dt*30*20; // in seconds
int steps = duration/dt;

cout << "Driver is running.\n";
cout << "Simulating " << numberOfBodies << " bodies.\n";
srand(time(NULL)); // Initiate RNG

// Create Array of bodies
Body *bodies = new Body[numberOfBodies];
bodies[0] = new Body(0, 0, centralMass, 0, 0);
for (int i = 1; i < numberOfBodies; i++) {
	// generate random x within range
	double x = range/2 - rand() % (int)range;
	// generate random y within range of circular radius
	double yRange = 2*sqrt(((range/2)*(range/2)) - (x*x));
	double y = yRange/2 - rand() % (int)yRange;
	// find x and y components of velocity
	double r = sqrt((x*x) + (y*y));
	double velocity = sqrt(G*centralMass/r);
	double xV = velocity * (x/r);
	double yV = velocity * (y/r);
	double mass = massRange;
	bodies[i] = new Body(x, y, mass, -yV, xV);
}*/