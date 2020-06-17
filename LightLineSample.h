#pragma once
#include "utility.h"
/* 相互引用 可能有问题*/



class LightLineSample
{
public:
	Point L0;
	Vec3<double> Ldir;
	Vec3<double> Lvertical;
	double length;
	std::vector<Vec2<double>> lineSeg;
	Color illumination;
	Vec3<double> normal;


	LightLineSample(const Point& l0, const Vec3<double>& ldir, const Vec3<double>& lvert, double leng, const Color& ill, const Vec3<double>& norm) {
		L0 = l0;
		Ldir = ldir;
		Lvertical = lvert;
		length = leng;
		illumination = ill;

		normal = norm;

		lineSeg.push_back(Vec2<double>(0.0, length));
	}

	// generate line segments
	// TODO: 效率改进
	void VisibilityTest(const std::vector<std::shared_ptr<Triangle const>>& triangleList, const Point& inscPoint) {
		
		for (auto pTri : triangleList) {
			std::vector<Vec2<double>> newSegs;
			for (auto seg : lineSeg) {
				Point p1, p2;
				Triangle t1(inscPoint, L0 + seg.x * Ldir, L0 + seg.y * Ldir);
				if (TriangleIntersection(t1, *pTri, p1, p2)) {
					// 投影
					Vec3<double> xp1 = p1 - inscPoint; double lenXP1 = xp1.length(); xp1.normalize();
					Vec3<double> xp2 = p2 - inscPoint; double lenXP2 = xp2.length(); xp2.normalize();

					double p1LineDis = PointLineDistance(p1, L0, Ldir);
					double p2LineDis = PointLineDistance(p2, L0, Ldir);
					double xLineDis = PointLineDistance(inscPoint, L0, Ldir);

					double p1LineProjPointDis = (lenXP1 * p1LineDis) / (xLineDis - p1LineDis);
					double p2LineProjPointDis = (lenXP2 * p2LineDis) / (xLineDis - p2LineDis);

					Point p1ProjOnLine = p1 + p1LineProjPointDis * xp1;
					Point p2ProjOnLine = p2 + p2LineProjPointDis * xp2;

					// L0 + t* dir = p1ProjOnLine
					Vec3<double> vp1 = (p1ProjOnLine - L0);
					Vec3<double> vp2 = (p2ProjOnLine - L0);


					double t1 = vp1.dot(Ldir);
					double t2 = vp2.dot(Ldir);



					if (t1 > t2) {
						Swap(t1, t2);
					}


					// seg.x seg.y  t1 t2 (t1 < t2)相交
					if (seg.y <= t1 || seg.x >= t2)
						newSegs.push_back(Vec2<double>(seg.x, seg.y));
					else if (seg.x < t1 && seg.y > t2) {
						newSegs.push_back(Vec2<double>(seg.x, t1));
						newSegs.push_back(Vec2<double>(t2, seg.y));
					}
					else if (seg.x >= t1 && seg.y <= t2) {
						continue;
					}
					else if (seg.x >= t1 && seg.x < t2) {
						newSegs.push_back(Vec2<double>(t2, seg.y));
					}
					else if (seg.y > t1 && seg.y <= t2) {
						newSegs.push_back(Vec2<double>(seg.x, t1));
					}
				}
				else {
					newSegs.push_back(seg);
				}
			}
			
			lineSeg = std::move(newSegs);
		}

		std::sort(lineSeg.begin(), lineSeg.end());
	}

	std::vector<Vec2<double>>& LineSegments() {
		return lineSeg;
	}
};