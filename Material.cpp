#include "Material.h"

/*
 * 2020-5-28  修改material接口；
			  对于相交点处 使用ReflImpSampling方法  以reference返回相交点处的：
			  1）albedo
			  2）brdf
			  3）根据brdf采样的方向
			  4) pdf
	material使用MatType来表示需要如何处理表面材质
			1）both
			2）transmit
			3）reflect
			4）emit
			5）subsurface
 */

Material::Material() {
	//matType = MatType::Matt;
}

// can be used to init emission or albedo
// no reflection on light
Material::Material(MatType m, const Vec3<double>& _albedoOrEmission) {
	
	matType = m;
	if (m == MatType::Emission)
		emit = _albedoOrEmission;
	else
		albedo = _albedoOrEmission;
}
Color Material::GetEmission() const{
	return emit;
}
Color Material::GetAlbedo()const{
	return albedo;
}

MatType Material::GetMatType()const {
	return matType;
}

void Material::SetAlbedo(const Color& v) {
	albedo = v;
}
void Material::SetEmission(const Color& c) {
	emit = c;
}




// makes no sense
Ray Material::RefrImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf) const {

	//Vec3<double> modified_norm;

	//if (normal.dot(r.dir) > 0)  //
	//	modified_norm = -1 * normal;
	//else
	//	modified_norm = normal;

	//Ray reflRay;
	//reflRay.ori = pos;
	//switch (matType) {
	//case MatType::Specular:
	//	reflRay.dir = 2 * (r.dir - r.dir.dot(modified_norm) * modified_norm) - r.dir;
	//	reflRay.dir.normalize();
	//	break;
	//case MatType::Matt:

	//	reflRay.dir = SampleOnSemiSphere(modified_norm);

	//	break;
	//case MatType::Glass:
	//	reflRay.dir = 2 * (r.dir - r.dir.dot(modified_norm) * modified_norm) - r.dir;
	//	reflRay.dir.normalize();

	//	break;
	//case MatType::Emission:
	//	reflRay.SetValidValue(false);
	//	break;
	//}

	//return reflRay;

	return Ray();
}

Ray Material::ReflImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf) const {

	//Vec3<double> modified_norm;

	//if (normal.dot(r.dir) > 0) {  //
	//	modified_norm = -1 * normal;

	//}
	//else
	//	modified_norm = normal;


	//double realRefrRatio = refrRatio;
	//if (inside) {
	//	realRefrRatio = 1 / realRefrRatio;
	//}


	//Ray refrRay;
	//refrRay.ori = pos;
	//switch (matType) {
	//case MatType::Glass:
	//{
	//	double cosTheta = -1 * modified_norm.dot(r.dir);
	//	double refrSin = cos2sin(cosTheta) / realRefrRatio;
	//	if (refrSin >= 1)
	//		refrRay.SetValidValue(false);
	//	else {
	//		double refrCos = sin2cos(refrSin);
	//		Vec3<double> horizontal = refrSin * (r.dir + cosTheta * modified_norm).normalize();
	//		refrRay.dir = horizontal - refrCos * modified_norm;
	//	}

	//	break;
	//}
	//default:
	//	refrRay.SetValidValue(false);
	//}

	//return refrRay;

	return Ray();
}

// makes no sense
Ray Material::RefrImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf, Color& albedo) const {

	return Ray();
}

Ray Material::ReflImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf, Color& albedo) const {

	return Ray();
}

double Material::BTDF(const Vec3<double>& wi, const Vec3<double>& wo, const Vec3<double>& norm) const {
	return 0.0;
}
double Material::BRDF(const Vec3<double>& wi, const Vec3<double>& wo, const Vec3<double>& norm) const {
	return 0.0;
}



////double Material::IntegrateOverLightLineSample(double l)
Color Material::AnalyticLineIllumination(const std::vector<Vec2<double>>& lineSeg, const LightLineSample& lls, const Point& inscPoint, const Vec3<double>& inscNorm, const Vec3<double>& lightNorm, const Vec3<double>& incidentDir)const {
	return Color(0.0, 0.0, 0.0);
}

