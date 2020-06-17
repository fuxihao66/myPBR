#pragma once
#include <cmath>
#include <ctime>
#include <algorithm>
#include <cstring>
#include "happly.h"

typedef unsigned int UINT;
#define PI 3.1415926
#define EPSILON (10e-7)

struct Config {
	UINT WIDTH;
	UINT HEIGHT;
	UINT SSP;
	UINT THREAD;
	Config() {
		WIDTH = 480;
		HEIGHT = 480;
		SSP = 100;
		THREAD = 1;
	}
	Config(UINT w, UINT h, UINT s, UINT t) {
		WIDTH = w;
		HEIGHT = h;
		SSP = s;
		THREAD = t;
	}
};

template <class T>
struct Vec2 {
	T x;
	T y;
	Vec2() {
		x = static_cast<T>(0);
		y = static_cast<T>(0);;
	}
	Vec2(const T& _x, const T& _y){
		x = _x;
		y = _y;
	}
	bool operator<(const Vec2<T>& v) {
		if (x > v.x)
			return true;
		else
			return false;
	}
};


template <class T>
struct Vec3 {
	T x;
	T y;
	T z;
	Vec3() {
		x = static_cast<T>(0);
		y = static_cast<T>(0);
		z = static_cast<T>(0);
	}
	Vec3(const T& _x, const T& _y, const T& _z) {
		x = _x;
		y = _y;
		z = _z;
	}	

	bool Parallel(const Vec3<double>& v) const{
		auto mul = (*this).OuterProd(v);
		if (abs(mul.x) < EPSILON && abs(mul.y) < EPSILON && abs(mul.z) < EPSILON)
			return true;
		else
			return false;
	}

	 Vec3<T>& normalize() {
		double norm = static_cast<double>(sqrt(x * x + y * y + z * z));
		x = static_cast<double>(x / norm);
		y = static_cast<double>(y / norm);
		z = static_cast<double>(z / norm);

		return *this;
	}
	T dot(const Vec3<T>& anotherVec) const{
		return x * anotherVec.x + y * anotherVec.y + z * anotherVec.z;
	}

	 Vec3<T> OuterProd(const Vec3<T>& v2) const {
		return Vec3<T>(y * v2.z - v2.y * z, v2.x * z - x * v2.z, x * v2.y - v2.x * y);
	}

	T length() const{
		return sqrt(x * x + y * y + z * z);
	}

	friend  Vec3<T> operator *(T para, const Vec3<T>& v2) {
		return Vec3<T>(para * v2.x, para * v2.y, para * v2.z);
	}

	friend  Vec3<T> operator +(const Vec3<T>& v1, const Vec3<T>& v2) {
		return Vec3<T>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
	}

	friend  Vec3<T> operator -(const Vec3<T>& v1, const Vec3<T>& v2) {
		return Vec3<T>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
	}

	friend  Vec3<T> operator *(const Vec3<T>& v1, const Vec3<T>& v2) {
		return Vec3<T>(v1.x * v2.x, v1.y*v2.y, v1.z*v2.z);
	}
};

typedef Vec3<double> Color;
typedef Vec3<double> Point;
typedef UINT Mat;

inline double distance(const Vec3<double>& p1, const Vec3<double>& p2) {
	return (p1 - p2).length();
}

inline double clamp(double x) { return x < 0 ? 0 : x>1 ? 1 : x; }
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }

inline double cos2sin(double cosValue) {
	//assert(cosValue >= 0);
	return sqrt(1 - cosValue * cosValue);
}
inline double sin2cos(double sinValue) {
	//assert(sinValue >= 0);
	return sqrt(1 - sinValue * sinValue);
}

inline double SampleRandomNumOverride(double start, double end, UINT seed) {
	std::srand(seed*(unsigned)time(0));
	return (double)rand() / RAND_MAX * (end - start) + start;
}

inline double SampleRandomNum(double start, double end) {
	//std::srand((unsigned)time(0));
	return (double)rand() / RAND_MAX * (end - start) + start;
}
/*
inline Vec3<double> SampleOnSemiSphere(const Vec3<double>& normal) {
	//std::srand((unsigned)time(0));
	//double phi = (double)rand() / RAND_MAX * PI / 2;

	//std::srand((unsigned)time(0));
	//double theta = (double)rand() / RAND_MAX * PI * 2;

	double r1 = 2 * PI * SampleRandomNum(0, 1.0f);
	double r2 = SampleRandomNum(0, 1.0f);
	double r2s = sqrt(r2);
	Vec3<double> w = normal;
	Vec3<double> u = ( (fabs(w.x) > .1 ? Vec3<double>(0, 1, 0) : Vec3<double>(1, 0, 0)).OuterProd( w )  ).normalize();
	auto v = w.OuterProd(u);
	Vec3<double> d = (cos(r1) * r2s * u + sin(r1) * r2s * v + sqrt(1 - r2) * w).normalize();

	return d;
}
*/
inline void resetRandomGen() {
	std::srand((unsigned)time(0));
}

// 半球上的均匀采样
inline Vec3<double> SampleOnSemiSphereUniform(const Vec3<double>& normal, const Vec3<double>& incidentDir) {
	double xi1 = SampleRandomNum(0, 1);
	double xi2 = SampleRandomNum(0, 1);
	double x = cos(2 * PI * xi2) * sqrt(1 - xi1 * xi1);
	double y = sin(2 * PI * xi2) * sqrt(1 - xi1 * xi1);
	double z = xi1;

	Vec3<double> newZAxis = normal;
	Vec3<double> newXAxis = incidentDir.OuterProd(normal).normalize();
	Vec3<double> newYAxis = newZAxis.OuterProd(newXAxis).normalize();


	// 坐标不变 但是基底改变

	double x_prime = x * newXAxis.x + y * newYAxis.x + z * newZAxis.x;
	double y_prime = x * newXAxis.y + y * newYAxis.y + z * newZAxis.y;
	double z_prime = x * newXAxis.z + y * newYAxis.z + z * newZAxis.z;

	return Vec3<double>(x_prime, y_prime, z_prime);
}

// phong材质的采样
inline Vec3<double> SamplePhongRefl(const Vec3<double>& normal, const Vec3<double>& mirrReflDir, UINT coeff) {
	double xi1 = SampleRandomNum(0, 1);
	double xi2 = SampleRandomNum(0, 1);
	double x = cos(2 * PI * xi2) * sqrt(1 - pow(xi1, 2.0 / (coeff + 1)));
	double y = sin(2 * PI * xi2) * sqrt(1 - pow(xi1, 2.0 / (coeff + 1)));
	double z = pow(xi1, 1.0 / (coeff + 1));


	Vec3<double> newZAxis = mirrReflDir;
	Vec3<double> newYAxis = normal.OuterProd(mirrReflDir).normalize();
	
	Vec3<double> newXAxis = newYAxis.OuterProd(newZAxis).normalize();
	double x_prime = x * newXAxis.x + y * newYAxis.x + z * newZAxis.x;
	double y_prime = x * newXAxis.y + y * newYAxis.y + z * newZAxis.y;
	double z_prime = x * newXAxis.z + y * newYAxis.z + z * newZAxis.z;

	return Vec3<double>(x_prime, y_prime, z_prime);
}

//inline Vec3<double> SampleOnSemiSphere(const Vec3<double>& normal) {
//	int index = 0;
//	while (true) {
//		index++;
//		double phi = SampleRandomNum(0.0, PI);
//		double theta = SampleRandomNum(0.0, 2 * PI);
//		Vec3<double> d(sin(phi)*cos(theta), sin(phi)*sin(theta),cos(phi));
//		if (d.dot(normal) >= 0)
//		return d;
//	}
//
//}


//inline Vec3<double> SampleScatterDir() {
//	double phi = SampleRandomNum(0.0, PI);
//	double theta = SampleRandomNum(0.0, 2 * PI);
//	Vec3<double> d(sin(phi) * cos(theta), sin(phi) * sin(theta), cos(phi));
//	return d;
//}

enum MatType {
	Reflect,
	//Specular,
	Transmit,
	Both,
	Emission,
	Subsurface
};

struct Photon {
	Vec3<double> position;
	Vec3<double> spectrum;
	Vec3<double> direction;

};


inline void LoadMeshFromFile(const std::string& path, std::vector<Vec3<size_t>>& face, std::vector<Point>& vert, std::vector<Point>& norms) {


	// read
	happly::PLYData plyIn(path);
	std::vector<std::array<double, 3>> vPos = plyIn.getVertexPositions();
	std::vector<std::vector<size_t>> fInd = plyIn.getFaceIndices<size_t>();


	// copy data
	face.resize(fInd.size());
	vert.resize(vPos.size());
	
	memcpy(&vert[0], &vPos[0], vPos.size() * 3 * sizeof(double));
	for (size_t i = 0; i < fInd.size(); i++) {
		face[i].x = fInd[i][0];
		face[i].y = fInd[i][1];
		face[i].z = fInd[i][2];
	}

	// compute normal
	norms.resize(face.size());
	
	for (size_t i = 0; i < face.size(); i++) {
		Point& p1 = vert[face[i].x];
		Point& p2 = vert[face[i].y];
		Point& p3 = vert[face[i].z];

		norms[i] = (p2 - p1).OuterProd(p3 - p1).normalize();
	}
}



inline bool PointInsideTriangle(const Point& testPoint, const Point& p1, const Point& p2, const Point& p3) {
	Vec3<double> p1p3 = p3 - p1;
	Vec3<double> p1p2 = p2 - p1;
	Vec3<double> p1p = testPoint - p1;

	double u = (p1p2.dot(p1p2)*p1p.dot(p1p3)-p1p3.dot(p1p2)*p1p.dot(p1p2))/(p1p3.dot(p1p3)*p1p2.dot(p1p2)-p1p3.dot(p1p2)*p1p2.dot(p1p3));
	double v = (p1p3.dot(p1p3)*p1p.dot(p1p2)-p1p3.dot(p1p2)*p1p.dot(p1p3)) / (p1p3.dot(p1p3) * p1p2.dot(p1p2) - p1p3.dot(p1p2) * p1p2.dot(p1p3));

	if (u >= 0 && u <= 1 && v >= 0 && v <= 1 && u + v <= 1) 
		return true;
	else
		return false;
}


inline bool NumInRange(double num, double top, double bot) {
	if (top < bot) {
		double temp = bot;
		bot = top;
		top = temp;
	}
	if (num <= top && num >= bot)
		return true;
	else
		return false;
}

inline double GetLineRecInsec(double p1u, double spp1dotLdir, double p2u, double spp2dotLdir) {
	if (p1u < EPSILON)
		return spp1dotLdir;
	if (p2u < EPSILON)
		return spp2dotLdir;
	double biggerOne, smallerOne;
	double biggerU,  smallerU;
	if (spp1dotLdir > spp2dotLdir) {
		biggerOne = spp1dotLdir;
		biggerU = p1u;
		smallerOne = spp2dotLdir;
		smallerU = p2u;
	}
	else {
		biggerOne = spp2dotLdir;
		biggerU = p2u;
		smallerOne = spp1dotLdir;
		smallerU = p1u;
	}
	return (biggerOne * smallerU + biggerU * smallerOne) / (smallerU + biggerU);
}

// only used for line lighting
// TODO: 这个算法有问题  l=0时结果不是0
//inline double DiffuseLineIntegrate(double A, double B, double C, double D, double F, double l) {
//	return (atan(D / sqrt(C - D * D)) - atan((D + l) / sqrt(C - D * D))) * (A - B * D) / (2 * pow(C - D * D, 1.5)) - F * (l * B * C * (C + l) + A * (C - D * l - 2 * D * D)) / (2 * C * (C - D * D) * (l * l + 2 * D * l + C));
//}

inline double DiffuseLineIntegrate(double A, double B, double C, double D, double F, double l) {
	double result = (atan(D / sqrt(C - D * D)) - atan((D + l) / sqrt(C - D * D))) * (A - B * D) / (2 * pow(C - D * D, 1.5)) - ((B * C - A * D) * l * l + (A * C + B * C * D - 2 * A * D * D) * l) / (2 * C * (C - D * D) * (l * l + 2 * D * l + C));
	return abs(result * F);
}

struct Triangle {
	Point p1, p2, p3;
	Vec3<double> normal;
	Vec3<double> p2p1;
	Vec3<double> p2p3;
	Triangle(const Point& _p1, const Point& _p2, const Point& _p3) {
		p1 = _p1;
		p2 = _p2;
		p3 = _p3;
		p2p1 = p1 - p2;
		p2p3 = p3 - p2;

		normal = p2p1.OuterProd(p2p3);
		normal.normalize();
	}
	
};







inline double PointPlaneRelation(const Point& pOnPlane, const Vec3<double>& planeNorm, const Point& p) {
	// 平面方程
	// nx*x + ny*y + nz*z - D = 0
	double D = planeNorm.x * pOnPlane.x + planeNorm.y * pOnPlane.y + planeNorm.z * pOnPlane.z;

	return planeNorm.x * p.x + planeNorm.y * p.y + planeNorm.z * p.z - D;
}


struct Pair {
	double rela;
	Point p;
	Pair(double _rela, const Point& _p) {
		rela = _rela;
		p = _p;
	}
	bool operator<(const Pair& p1)
	{
		return rela < p1.rela;
	}
};


inline Point LinePlaneInsec(const Vec3<double>& dir, const Point& pLine, const Vec3<double>& norm, const Point& pOnPlane) {
	double t = (pOnPlane.dot(norm) - pLine.dot(norm)) / dir.dot(norm);
	return pLine + t * dir;
}

inline bool TrianglePlaneIntersection(const Triangle& t, const Point& pOnPlane, const Vec3<double>& planeNorm, Point& p1, Point& p2) {
	std::vector<Pair> triPointWithPlaneRela{
		Pair(PointPlaneRelation(pOnPlane, planeNorm, t.p1), t.p1),
		Pair(PointPlaneRelation(pOnPlane, planeNorm, t.p2), t.p2),
		Pair(PointPlaneRelation(pOnPlane, planeNorm, t.p3), t.p3)
	};



	std::sort(triPointWithPlaneRela.begin(), triPointWithPlaneRela.end());

	
	if (triPointWithPlaneRela[0].rela >= EPSILON || triPointWithPlaneRela[2].rela <= 0)
		return false;
	
	Vec3<double> dir1, dir2;
	if (triPointWithPlaneRela[1].rela > EPSILON) {
		// <0 >0 >0
		dir1 = triPointWithPlaneRela[1].p - triPointWithPlaneRela[0].p;
		dir2 = triPointWithPlaneRela[2].p - triPointWithPlaneRela[0].p ;
		p1 = LinePlaneInsec(dir1, triPointWithPlaneRela[0].p, planeNorm, pOnPlane);
		p2 = LinePlaneInsec(dir2, triPointWithPlaneRela[0].p, planeNorm, pOnPlane);

	}
	else {
		// <0 <0 >0
		dir1 = triPointWithPlaneRela[0].p - triPointWithPlaneRela[2].p;
		dir2 = triPointWithPlaneRela[1].p - triPointWithPlaneRela[2].p;
		p1 = LinePlaneInsec(dir1, triPointWithPlaneRela[2].p, planeNorm, pOnPlane);
		p2 = LinePlaneInsec(dir2, triPointWithPlaneRela[2].p, planeNorm, pOnPlane);
	}
	return true;
}
// only used for line lighting
inline bool TriangleIntersection(const Triangle& t1, const Triangle& t2, Point& p1, Point& p2) {
	// 转化成三角形和平面交 平面和三角形交  最后两条线交
	Point t1p1, t1p2, t2p1, t2p2;

	auto t1InsecP2 = TrianglePlaneIntersection(t1, t2.p1, t2.normal, t1p1, t1p2);
	if (!t1InsecP2)
		return false;
	auto t2InsecP1 = TrianglePlaneIntersection(t2, t1.p1, t1.normal, t2p1, t2p2);
	if (!t2InsecP1)
		return false;
	// t1p1为原点
	Vec3<double> dir = t1p2 - t1p1;
	double len = dir.length();
	dir.normalize();

	double coort2p1 = (t2p1 - t1p1).dot(dir);
	double coort2p2 = (t2p2 - t1p1).dot(dir);
	
	double minT = std::max(0.0, std::min(coort2p1, coort2p2));
	double maxT = std::min(len, std::max(coort2p1, coort2p2));

	if (minT >= len || maxT <= 0.0) {
		return false;
	}
	else if (maxT - minT < EPSILON) {
		return false;
	}
	else {
		p1 = t1p1 + minT * dir;
		p2 = t1p1 + maxT * dir;
		

		return true;
	}
}

inline double ClampBeforeACosAndASin(double beforeClamp) {
	if (beforeClamp > 1.0)
		return 1.0;
	if (beforeClamp < -1.0)
		return - 1.0;
	return beforeClamp;
}

inline double IntegrateCosExpSin(size_t exponent, double from, double to) {
	return (pow(cos(from), exponent + 1) - pow(cos(to), exponent + 1)) / (exponent + 1);
}

// 计算 cos^n
inline double IntegrateCosExp(size_t exponent, double from, double to) {
	switch (exponent)
	{
	case 0:
		return to - from;
		break;
	case 1:
		return sin(to) - sin(from);
		break;
	case -1:
		std::cout << "illegal exponent when doing integral" << std::endl;
		return -0xFFFFF;
		break;
	default: // >= 2
		return (pow(cos(to), exponent - 1) * sin(to) - pow(cos(from), exponent - 1) * sin(from)) / exponent + double(exponent - 1) / exponent * IntegrateCosExp(exponent - 2, from, to);
	}
}


inline double PointLineDistance(const Point& p, const Point& L0, const Vec3<double>& dir) {
	double t = dir.dot(L0) - dir.dot(p);
	return (p - L0 - t * dir).length();
}

inline double Sign(double t) {
	if (t > 0)
		return 1;
	else
		return -1;
}

inline void Swap(double& t1, double& t2) {
	double temp = t1;
	t1 = t2;
	t2 = temp;
}

// matrix
