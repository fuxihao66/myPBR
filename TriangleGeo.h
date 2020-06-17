#pragma once
#include "Geometry.h"
#include "utility.h"

class TriangleGeo : public Geometry
{
private:
	Point p1, p2, p3;
	Vec3<double> normal;
	Vec3<double> p2p1;
	Vec3<double> p2p3;
public:
	TriangleGeo(const Point& _p1, const Point& _p2, const Point& _p3);

	virtual bool Intersect(const Ray&, double&, Vec3<double>&)const ;
	virtual double SurfaceArea() const ;
	virtual Point RandomSampleOnSurf() const ;

	std::shared_ptr<Triangle> ToTriangle()const;
};

