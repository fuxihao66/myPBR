#pragma once
#include "Geometry.h"

class Sphere : public Geometry
{
private:
	Vec3<double> position; 
	double radius;

	
public:
	Sphere();
	Sphere(double x, double y, double z, double r);

	//virtual const Vec3<double> GetNormal(const Vec3<double>& pos) const;
	virtual bool Intersect(const Ray&, double&, Vec3<double>&) const;
	//virtual bool Contain(const Vec3<double>& point)const;

};

