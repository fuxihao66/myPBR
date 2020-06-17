#include "Sphere.h"

Sphere::Sphere() {

}
Sphere::Sphere(double x, double y, double z, double r) {
	position.x = x;
	position.y = y;
	position.z = z;
	radius = r;
}
//const Vec3<double> Sphere::GetNormal(const Vec3<double>& hitPos) const{
//	Vec3<double> m = hitPos - position;
//	return m.normalize();
//}
bool Sphere::Intersect(const Ray& r, double& t, Vec3<double>& norm)const {
	double closestPoint = (position.dot(r.dir) - r.ori.dot(r.dir)) / r.dir.dot(r.dir);

	double distance = (r.ori + closestPoint * r.dir - position).length();

	double offset = sqrt(radius * radius - distance * distance);


	Vec3<double> ray2SphereCenter = position - r.ori;

	if (ray2SphereCenter.length() > radius) {
		if (distance > radius)
			return false;
		else if (ray2SphereCenter.dot(r.dir) < 0) 
			return false;
		else {
			t = closestPoint - offset;
		}
	}
	else {
		t = closestPoint + offset;
	}
	
	if (t < 0 || t < EPSILON)
		return false;

	norm = r.ori + t * r.dir - position;
	norm.normalize();

	return true;

}

//bool Sphere::Contain(const Vec3<double>& point)const {
//	if (distance(point, position) <= radius)
//		return true;
//	else
//		return false;
//}