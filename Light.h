#pragma once
#include "RenderableObject.h"
#include "Plane.h"
#include "LightLineSample.h"

class Light : public RenderableObject
{
public:
	Light();
	Light(std::shared_ptr<Material const> mat, std::shared_ptr<Geometry const> geo);
	Ray SamplePoint(const Point& , double& tMax);
	LightLineSample SampleLine(double& samplePDF);
	double FuncG(double cosTheta, double cosAlpha, double pxLeng); // render equation 面积分形式中的G
	double PointSamplePDF();
};

