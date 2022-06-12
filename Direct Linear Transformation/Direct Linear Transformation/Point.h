#pragma once

class Point2f
{
public:
	double mx = 0;
	double my = 0;

public:
	Point2f() {};
	Point2f(double x, double y) :mx(x), my(y) {};
};


class Point3f
{
public:
	double mx = 0;
	double my = 0;
	double mz = 0;

public:
	Point3f() {};
	Point3f(double x, double y, double z) :mx(x), my(y), mz(z) {};
};