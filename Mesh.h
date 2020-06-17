#pragma once
#include <vector>
#include "Geometry.h"
#include "Ray.h"

class Mesh: public Geometry
{
private:
	UINT numVert;
	UINT numFace;
	std::vector<Vec3<size_t>> faceIndices;
	std::vector<Point> vertices;
	std::vector<Point> normals;// norm of every face

	double TriangleInsec(const Ray&,size_t faceId) const;
public:
	Mesh();
	Mesh(const std::string& path);
	virtual bool Intersect(const Ray&, double&, Vec3<double>&) const;

	//virtual const Vec3<double> GetNormal(const Vec3<double>& pos) const;
	//virtual bool Contain(const Vec3<double>& point) const;
	virtual double SurfaceArea() const;
	virtual Point RandomSampleOnSurf() const;
};

