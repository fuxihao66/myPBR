#pragma once
#include "Material.h"

class PerfectDiffuse : public Material
{
public:
	PerfectDiffuse(MatType m, const Vec3<double>& _albedo);

	virtual Ray ReflImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf, Color& albedo) const;

	virtual Ray RefrImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf, Color& albedo) const;

	virtual Ray ReflImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf) const;

	virtual Ray RefrImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf) const;

	virtual double BTDF(const Vec3<double>& wi, const Vec3<double>& wo, const Vec3<double>& norm) const;
	virtual double BRDF(const Vec3<double>& wi, const Vec3<double>& wo, const Vec3<double>& norm) const;


	virtual Color AnalyticLineIllumination(const std::vector<Vec2<double>>& lineSeg, const LightLineSample& lls, const Point& inscPoint, const Vec3<double>& inscNorm, const Vec3<double>& lightNorm, const Vec3<double>& incidentDir)const;
};

