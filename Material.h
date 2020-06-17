#pragma once
#include <iostream>
#include "utility.h"
#include "Ray.h"
#include "LightLineSample.h"

// 使用继承的方式重载
// TODO: 为emission专门派生一种material 发光体同时也可以反射
// TODO: make Material pure virtual

class Material
{
protected:
	Color albedo;
	Color emit;
	MatType matType;
public:
	Material();
	Material(MatType m, const Vec3<double>& );

	Color GetEmission() const;
	Color GetAlbedo() const;
	void SetAlbedo(const Color& v);
	void SetEmission(const Color& c);
	MatType GetMatType() const;


	virtual Color AnalyticLineIllumination(const std::vector<Vec2<double>>& lineSeg, const LightLineSample& lls, const Point& inscPoint, const Vec3<double>& inscNorm, const Vec3<double>& lightNorm, const Vec3<double>& incidentDir)const;
	// for material imp sampling
	virtual Ray ReflImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf, Color& albedo) const;

	virtual Ray RefrImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf, Color& albedo) const;

	virtual Ray ReflImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf) const;

	virtual Ray RefrImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf) const;

	// can be called to compute direct illumination
	virtual double BTDF(const Vec3<double>& wi, const Vec3<double>& wo, const Vec3<double>& norm) const;
	virtual double BRDF(const Vec3<double>& wi, const Vec3<double>& wo, const Vec3<double>& norm) const;
};

