#include "Mesh.h"


Mesh::Mesh() {
	numVert = 0;
	numFace = 0;
}
Mesh::Mesh(const std::string& path) {
	LoadMeshFromFile(path, faceIndices, vertices, normals);
	numVert = vertices.size();
	numFace = faceIndices.size();
}

double Mesh::TriangleInsec(const Ray& r, size_t faceId) const {
	const Point& norm = normals[faceId];
	if (norm.dot(r.dir) == 0)
		return -1;
	const Point& a = vertices[faceIndices[faceId].x];
	const Point& b = vertices[faceIndices[faceId].y];
	const Point& c = vertices[faceIndices[faceId].z];

	double t = (a - r.ori).dot(norm) / r.dir.dot(norm);

	if (t <= 0)
		return t;

	if (PointInsideTriangle(r.ori + t * r.dir, a, b, c))
		return t;
	else
		return -1;

}


//const Vec3<double> Mesh::GetNormal(const Vec3<double>& pos) const {
//
//}
bool Mesh::Intersect(const Ray& r, double& t, Vec3<double>& norm) const {
	double local_t = 0xFFFFF;
	size_t selectedFace = -1;


	for (size_t i = 0; i < numFace; i++) {
		double temp_t = TriangleInsec(r, i);

		if (temp_t < EPSILON)  // 在发出点相交或不相交
			continue;

		if (local_t > temp_t) {
			selectedFace = i;
			local_t = temp_t;
		}
	}

	if (selectedFace == -1)
		return false;
	else {
		t = local_t;
		norm = normals[selectedFace];
	}
}
//bool Mesh::Contain(const Vec3<double>& point) const {
//
//}

double Mesh::SurfaceArea() const {
	return 0.0;
}
Point Mesh::RandomSampleOnSurf() const {
	return Point();
}