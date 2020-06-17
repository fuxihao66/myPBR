#include "Plane.h"

Plane::Plane() {

}

Plane::Plane(const Vec3<double>& _p1, const Vec3<double>& _p2, const Vec3<double>& _p3) {
	p1 = _p1;
	p2 = _p2;
	p3 = _p3;
	p2p1 = p1 - p2;
	p2p3 = p3 - p2;
	p2p4 = p2p1 + p2p3;
	p4 = p2 + p2p4;
	w = p2p1.length();
	h = p2p3.length();
	normal = p2p3.OuterProd(p2p1); //
	normal.normalize();
}
//
Vec3<double> Plane::GetNormal() const{
	return normal;
}


bool Plane::Intersect(const Ray& r, double& t, Vec3<double>& norm)const {
	if (normal.dot(r.dir) == 0)
		return false;

	t = (normal.dot(p2) - normal.dot(r.ori)) / normal.dot(r.dir);
	
	if (t < 0)
		return false;

	if (abs(t) < EPSILON)
		return false;

	norm = normal;

	Vec3<double> p2InsecPoint = r.ori + t * r.dir - p2;


	double innerProduct1 = p2InsecPoint.dot(p2p1);
	double innerProduct2 = p2InsecPoint.dot(p2p3);


	if (innerProduct1 > w* w || innerProduct1 < 0 || innerProduct2 > h * h || innerProduct2 < 0)
		return false;
	else
		return true;

}

//bool Plane::Contain(const Vec3<double>& point) const {
//	return false;
//}

double Plane::SurfaceArea() const{

	return w*h;
}

Point Plane::RandomSampleOnSurf() const {
	double coorX = SampleRandomNum(0.0, 1.0);
	double coorY = SampleRandomNum(0.0, 1.0);

	return p2 + coorX * p2p1 + coorY * p2p3;
}


// 注意返回的L0是世界坐标的
// diffuse处理是相交点为原点的坐标
void Plane::SampleLine(Point& L0, Vec3<double>& Ldir, Vec3<double>& Lvert, double& length, double& pdf) const {
	// 随机选取坐标系
	double theta = SampleRandomNum(0, 2 * PI);
	double cosTheta = cos(theta);
	double sinTheta = sin(theta);

	Vec3<double> xDir = p2p1;   xDir.normalize();
	Vec3<double> yDir = p2p3;   yDir.normalize();
	Ldir = cosTheta * xDir + sinTheta * yDir;

	Lvert = Ldir.OuterProd(normal);
	Lvert.normalize();
	// p2 为原点   Lvert Ldir分别为 x y

	// a < b < c < d
	double a, b, c, d, l;

	// l是最长的line
	if (abs(cosTheta) < EPSILON)
		l = h / abs(sinTheta);
	else if (abs(sinTheta) < EPSILON)
		l = w / abs(cosTheta);
	else
		l = std::min(w / abs(cosTheta), h / abs(sinTheta));

	double p2u = 0;
	double p1u = p2p1.dot(Lvert);
	double p3u = p2p3.dot(Lvert);
	double p4u = p2p4.dot(Lvert);
	std::vector<double> fourPoints{p2u, p1u, p3u, p4u};
	std::sort(fourPoints.begin(), fourPoints.end());
	a = fourPoints[0]; b = fourPoints[1]; c = fourPoints[2]; d = fourPoints[3];

	///// 重要性采样

	double area = (c - b + d - a) * l / 2.0;


	double t;
	// axis-aligned
	if (abs(a-b) < EPSILON) { 
		t = SampleRandomNum(a, d);
		L0 = p2 + t * Ldir;
		length = l;

		pdf = l / area;

		return;
	}

	// 否则
	double k = l / (b - a); // 斜率

	while (true) {
		double xi = SampleRandomNum(a, d);
		double lamb = SampleRandomNum(0, l);
		if (a <= xi && xi <= b) {
			if (lamb > (xi - a) * k)
				continue;
		}
		else if (c <= xi && xi <= d) {
			if (lamb > (d - xi) * k)
				continue;
		}
		t = xi;
		break;
	}

	Point startPoint = p2 + t * Lvert; // 原点在p2上
	double scale = 0xFFFFF;
	

	// 判断和矩形哪条边有交点
	// TODO: 修改
	if (NumInRange(t, p2u, p1u)) {
		Vec3<double> spp2 = p2 - startPoint;
		Vec3<double> spp1 = p1 - startPoint;
		scale = std::min(scale, GetLineRecInsec(abs(p1u-t), spp1.dot(Ldir), abs(p2u-t), spp2.dot(Ldir)) );
	}
	if (NumInRange(t, p2u, p3u)) {
		Vec3<double> spp2 = p2 - startPoint;
		Vec3<double> spp3 = p3 - startPoint;
		scale = std::min(scale, GetLineRecInsec(abs(p2u - t), spp2.dot(Ldir), abs(p3u - t), spp3.dot(Ldir)));
	}
	if (NumInRange(t, p4u, p3u)) {
		Vec3<double> spp4 = p4 - startPoint;
		Vec3<double> spp3 = p3 - startPoint;
		scale = std::min(scale, GetLineRecInsec(abs(p4u - t), spp4.dot(Ldir), abs(p3u - t), spp3.dot(Ldir)));
	}
	if (NumInRange(t, p4u, p1u)) {
		Vec3<double> spp4 = p4 - startPoint;
		Vec3<double> spp1 = p1 - startPoint;
		scale = std::min(scale, GetLineRecInsec(abs(p1u - t), spp1.dot(Ldir), abs(p4u - t), spp4.dot(Ldir)));
	}
	L0 = startPoint + scale * Ldir;

	if (a <= t && t <= b) {
		length = (t - a) * k;
	}
	else if (c <= t && t <= d) {
		length = (d - t) * k;
	}
	else
		length = l;
	pdf = length / area;

}