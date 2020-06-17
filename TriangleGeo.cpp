#include "TriangleGeo.h"

TriangleGeo::TriangleGeo(const Point& _p1, const Point& _p2, const Point& _p3) {
	p1 = _p1;
	p2 = _p2;
	p3 = _p3;
	p2p1 = p1 - p2;
	p2p3 = p3 - p2;
	normal = p2p1.OuterProd(p2p3);
	normal.normalize();
}

std::shared_ptr<Triangle> TriangleGeo::ToTriangle()const {
	return std::make_shared<Triangle>(p1, p2, p3);
}

bool TriangleGeo::Intersect(const Ray& r, double& t, Vec3<double>& norm)const {
	if (normal.dot(r.dir) == 0)
		return -1;

	t = (p1 - r.ori).dot(normal) / r.dir.dot(normal);

	if (t <= EPSILON)
		return false;

	norm = normal;

	if (PointInsideTriangle(r.ori + t * r.dir, p1, p2, p3))
		return true;
	else
		return false;
}


// TODO: delete in later version
double TriangleGeo::SurfaceArea() const {
	return 0.0;
}

// TODO: delete in later version
Point TriangleGeo::RandomSampleOnSurf() const {
	return Point();
}