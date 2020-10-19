#ifndef __DOUBLEDOUBLE_H__
#define __DOUBLEDOUBLE_H__

struct DoubleDouble {
	double x;
	double y;

	DoubleDouble() {
		x = 0;
		y = 0;
	}

	DoubleDouble(double _x, double _y) {
		x = _x;
		y = _y;
	}
};

#endif