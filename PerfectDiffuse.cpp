#include "PerfectDiffuse.h"


PerfectDiffuse::PerfectDiffuse(MatType m, const Vec3<double>& _albedo) : Material(m,  _albedo) {

}

Ray PerfectDiffuse::ReflImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf, Color& albedo) const {
	albedo = GetAlbedo();

	return ReflImpSampling(r, pos, normal, pdf, brdf);
}

Ray PerfectDiffuse::RefrImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf, Color& albedo) const {

	return Ray();
}

Ray PerfectDiffuse::ReflImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf) const {
	Vec3<double> modifiedNorm = normal;
	if (normal.dot(r.dir) > 0)
		modifiedNorm = -1*normal;

	Ray newRay(pos, SampleOnSemiSphereUniform(modifiedNorm, r.dir));

	// return by reference
	brdf = BRDF(r.dir, newRay.dir, modifiedNorm);
	
	pdf = 1 / (2 * PI);

	return newRay;
}

Ray PerfectDiffuse::RefrImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf) const {

	return Ray();
}

// no transmission
double PerfectDiffuse::BTDF(const Vec3<double>& wi, const Vec3<double>& wo, const Vec3<double>& norm) const {
	return 0.0;
}

double PerfectDiffuse::BRDF(const Vec3<double>& wi, const Vec3<double>& wo, const Vec3<double>& norm) const {
	return 1 / PI;
}


Color PerfectDiffuse::AnalyticLineIllumination(const std::vector<Vec2<double>>& lineSeg, const LightLineSample& lls, const Point& inscPoint, const Vec3<double>& inscNorm, const Vec3<double>& lightNorm, const Vec3<double>& incidentDir)const {
	Point newL0 = lls.L0 - inscPoint; // 变换坐标

	// norm 方向

	double A = newL0.dot(inscNorm);
	double B = lls.Ldir.dot(inscNorm);
	double C = newL0.dot(newL0);
	double D = newL0.dot(lls.Ldir);
	double E = 1;
	double F = newL0.dot(lightNorm);

	Color dirIllu;
	// 假设seg是降序排序  seg.x < seg.y
	for (auto seg : lineSeg) {
		dirIllu = dirIllu + DiffuseLineIntegrate(A, B, C, D, F, seg.y) * (1.0 / PI) * albedo * lls.illumination;
		dirIllu = dirIllu - DiffuseLineIntegrate(A, B, C, D, F, seg.x) * (1.0 / PI) * albedo * lls.illumination;
	}

	return dirIllu;
}
