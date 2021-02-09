#ifndef __POINT_H__
#define __POINT_H__

class Point {
	public:
	double x;
	double y;

	Point() {
		x = 0;
		y = 0;
	}

	Point(double _x, double _y) {
		x = _x;
		y = _y;
	}

	Point(Point* point) {
		x = point->x;
		y = point->y;
	}
};

#endif