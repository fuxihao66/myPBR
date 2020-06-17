#include "Light.h"


Light::Light() : RenderableObject(){
}

Light::Light(std::shared_ptr<Material const> mat, std::shared_ptr<Geometry const> geo) : RenderableObject(mat, geo) {

}

// new ray from _x to a sample point
Ray Light::SamplePoint(const Point& _x, double& tMax) {
	Point p = _geo->RandomSampleOnSurf();

	Vec3<double> dir = p - _x;

	tMax = dir.length();

	return Ray(_x, dir.normalize());
}


double Light::PointSamplePDF() {
	return 1.0/ _geo->SurfaceArea();
}


double Light::FuncG(double cosTheta, double cosAlpha, double pxLeng) {
	/*Vec3<double> px = x - PointOnLight;
	double pxLeng;
	Vec3<double> normLight;
	GetGeo()->Intersect(Ray(x, px.normalize()), pxLeng, normLight);
	double cosAlpha = px.dot(normLight);*/

	return abs(cosTheta * cosAlpha) / (pxLeng * pxLeng);
}



/*----------------仅用于line sampling-----------------*/

// TODO: 如何处理不同形状的光源  geometry增加虚函数？
// 
LightLineSample Light::SampleLine(double& samplePDF) {
	
	auto plane = std::dynamic_pointer_cast<Plane const>(GetGeo());

	Point L0;
	Vec3<double> Ldir, Lvert;
	double length;

	plane->SampleLine(L0, Ldir, Lvert, length, samplePDF);

	

	LightLineSample lls(L0, Ldir, Lvert, length, GetMat()->GetEmission(), plane->GetNormal());


	return lls;
}

