#include "RenderableObject.h"


RenderableObject::RenderableObject() {

}
RenderableObject::RenderableObject(std::shared_ptr<Material const> m, std::shared_ptr<Geometry const> g) {
	_mat = m;
	_geo = g;
}
bool RenderableObject::Intersect(const Ray& r, double& t, Vec3<double>& norm) const{
	
	return _geo->Intersect(r, t, norm);
}

bool RenderableObject::Intersect(const Ray& r, double& t) const {
	Vec3<double> norm;
	return _geo->Intersect(r, t, norm);
}

//
//MatType RenderableObject::GetMatType() const{
//	return _mat->GetMatType();
//}
//
//Color RenderableObject::GetAlbedo() const{
//	return _mat->GetAlbedo();
//}
//
//Color RenderableObject::GetEmission() const {
//	return _mat->GetEmission();
//}

void RenderableObject::SetGeo(std::shared_ptr<Geometry const> g) {
	_geo = g;
}
void RenderableObject::SetMat(std::shared_ptr<Material const> m) {
	_mat = m;
}

//Ray RenderableObject::GetRefl(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf) const{
//
//	return _mat->ReflImpSampling(r, pos, normal, pdf);
//}
//Ray RenderableObject::GetRefr(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf) const{
//
//	return _mat->RefrImpSampling(r, pos, normal, pdf);
//	
//}

std::shared_ptr<Material const> RenderableObject::GetMat() const {
	return _mat;
}
std::shared_ptr<Geometry const> RenderableObject::GetGeo() const {
	return _geo;
}