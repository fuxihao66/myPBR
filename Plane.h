#pragma once
#include <algorithm>

#include "Ray.h"
#include "Geometry.h"

class Plane: public Geometry
{
private:
	Vec3<double> p1;  
	Vec3<double> p2;
	Vec3<double> p3;
	Vec3<double> p4;
	Vec3<double> p2p1;
	Vec3<double> p2p3;
	Vec3<double> p2p4;

	Vec3<double> normal;
	double w;
	double h;
	
public:
	//virtual const Vec3<double> GetNormal(const Vec3<double>& pos) const;
	virtual bool Intersect(const Ray&, double&, Vec3<double>&) const;
	//virtual bool Contain(const Vec3<double>& point) const;
	Plane();
	Plane(const Vec3<double>& _p1, const Vec3<double>& _p2, const Vec3<double>& _p3);
	virtual double SurfaceArea() const;
	virtual Point RandomSampleOnSurf() const;

	// ÎªÁËline sampling
	void SampleLine(Point& L0, Vec3<double>& Ldir, Vec3<double>& Lvert, double& length, double& pdf) const;

	Vec3<double> GetNormal() const;
};

