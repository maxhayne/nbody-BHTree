#ifndef __BODY_H__
#define __BODY_H__

#include "Point.h"
#include <vector>

class Body {
	public:
	double x;
	double y;
	double vX;
	double vY;
	double aX;
	double aY;
	double mass;

	std::vector<Point> locations;

	Body() {
		x = 0;
		y = 0;
		mass = 0;
		vX = 0;
		vY = 0;
		aX = 0;
		aY = 0;
	}

	Body(double _x, double _y, double _mass) {
		x = _x;
		y = _y;
		mass = _mass;
		vX = 0;
		vY = 0;
		aX = 0;
		aY = 0;
	}

	Body(Point* point, double _mass) {
		x = point->x;
		y = point->y;
		mass = _mass;
		vX = 0;
		vY = 0;
		aX = 0;
		aY = 0;
	}

	Body(double _x, double _y, double _mass, double _vX, double _vY) {
		x = _x;
		y = _y;
		mass = _mass;
		vX = _vX;
		vY = _vY;
		aX = 0;
		aY = 0;
	}

	Body(Point* point, double _mass, double _vX, double _vY) {
		x = point->x;
		y = point->y;
		mass = _mass;
		vX = _vX;
		vY = _vY;
		aX = 0;
		aY = 0;
	}

	Body(Body &b) {
		x = b.x;
		y = b.y;
		mass = b.mass;
		vX = b.vX;
		vY = b.vY;
		aX = b.aX;
		aY = b.aY;
	}

	Body(Body *b) {
		x = b->x;
		y = b->y;
		mass = b->mass;
		vX = b->vX;
		vY = b->vY;
		aX = b->aX;
		aY = b->aY;
	}

	// Naive step function
	void step(double dt) {
		locations.push_back(Point(x,y));
		vY += dt * aX;
	    vX += dt * aY;
	    x += dt * vX;
	    y += dt * vY;
	}

	// Call first, then calculate net accelerations
	void leapfrog_part1(double dt) {
		locations.push_back(Point(x,y));
		x  += 0.5 * dt * vX;
		y  += 0.5 * dt * vY;
	}

	// Call after calculating accelerations
	void leapfrog_part2(double dt) {
		vX += dt * aX;
		vY += dt * aY;
		x  += 0.5 * dt * vX;
		y  += 0.5 * dt * vY;
	}
};

#endif