#pragma once
#include <vector>
#include <ctime>
#include <iostream>

#include "RenderableObject.h"
#include "LightLineSample.h"

/*
 2020-5-28  改成Intersect只返回相交点的信息 不负责之后的采样等处理
*/

class ObjWarehouse
{
private:
	std::vector<std::shared_ptr<RenderableObject const >> objList;
public:
	std::shared_ptr<RenderableObject const> Intersect(const Ray& r, Vec3<double>& norm, Point& inscPoint)const {

		double closestPoint = 0xFFFFF;
		
		Vec3<double> temp_norm;

		std::shared_ptr<RenderableObject const> pClosestObject;
		for (auto it = objList.begin(); it != objList.end(); it++) {
			
			double t = 0xFFFFF;

			bool intersected = (*it)->Intersect(r, t, temp_norm);
			
			if (intersected && closestPoint > t) {
				closestPoint = t;
				pClosestObject = *it;
				norm = temp_norm;
			}
		}
		if (closestPoint == 0xFFFFF){
			
			return nullptr;
		}
		inscPoint = r.ori + closestPoint * r.dir;

		return pClosestObject;
		
	}

	// for shadow ray 
	bool Intersect(const Ray& r, double& t)const {

		double closestPoint = 0xFFFFF;


		for (auto it = objList.begin(); it != objList.end(); it++) {


			bool intersected = (*it)->Intersect(r, t);

			if (intersected && closestPoint > t) {
				closestPoint = t;
			}
		}
		if (closestPoint == 0xFFFFF) {

			return false;
		}

		t = closestPoint;

		return true;

	}

	/*void VisibilityTest(const LightLineSample& lls, const Point& inscPoint, std::vector<Vec2<double>>& lineSegs) const {


	}*/

	void AddObject(std::shared_ptr<RenderableObject const> obj) {
		objList.push_back(obj);
	}
};

