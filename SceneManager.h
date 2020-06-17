#pragma once

#include <cmath>
#include <iostream>
#include <thread>
#include "utility.h"
#include "Image.h"
#include "Ray.h"
#include "RenderableObject.h"
#include "ObjWarehouse.h"
#include "Camera.h"
#include "Material.h"
#include "PerfectDiffuse.h"
#include "Mirror.h"
#include "Light.h"
#include "Plane.h"
#include "Sphere.h"
#include "Phong.h"
#include "LightLineSample.h"

class SceneManager
{
private:
	ObjWarehouse _objList;

	// only used for line sampling
	std::vector<std::shared_ptr<Triangle const>> _triList;
	

	std::vector<std::shared_ptr<Light>> _lightList; // 直接光照
	Config _f;
	std::shared_ptr<Camera const> _pCam;

	// path tracing 
	Color DoRenderPathTracing(const Ray&, int depth) const;

	// 增加shadow ray 处理直接光照
	Color DoRenderPathTracingWithDirectIllumination(const Ray&, int depth) const;

	//std::shared_ptr<Image> DoRenderPhotonMapping() const;
	Color DirectIlluminationPoint(const Ray& r, const Point& inscPoint, const Vec3<double>& inscNormal, std::shared_ptr<Material const> mat, size_t sampleNum) const;
	Color DirectIlluminationLine(const Ray& r, const Point& inscPoint, const Vec3<double>& inscNormal, std::shared_ptr<Material const> mat) const;

	// 不包括end
	static void DoRenderHeightRange(const UINT heightStart, const UINT heightEnd, const UINT width, const UINT SSP, const Point& camPos, double pixelWidth, double pixelHeight, std::shared_ptr<Image> img, const SceneManager* psm);
public:
	SceneManager();

	bool BuildDefaultScene();  // 

	void AddTriangle(std::shared_ptr<Triangle const> obj);

	void AddObject(std::shared_ptr<const RenderableObject>);
	void SetConfig(const Config& f);
	void SetCamera(std::shared_ptr<const Camera>);
	std::shared_ptr<Image> DoRender() const;

};

