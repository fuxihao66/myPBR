#pragma once
#include <memory>
#include "Geometry.h"
#include "Material.h"

/*
* 2020-5-31 modify  RenderableObject不再代理Material和Geometry的接口
	直接提供获取material和geometry的接口
*/

class RenderableObject
{
protected:
	std::shared_ptr<Material const> _mat;
	std::shared_ptr<Geometry const> _geo;
public:
	RenderableObject();
	RenderableObject(std::shared_ptr<Material const> m, std::shared_ptr<Geometry const> g);
	/*Ray GetRefl(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf) const;
	Ray GetRefr(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf) const;*/
	void SetGeo(std::shared_ptr<Geometry const> );
	void SetMat(std::shared_ptr<Material const>);
	/*Color GetAlbedo() const;
	Color GetEmission() const;
	MatType GetMatType() const;*/
	bool Intersect(const Ray& r, double& t, Vec3<double>&) const;
	bool Intersect(const Ray& r, double& t) const;
	std::shared_ptr<Material const> GetMat() const;
	std::shared_ptr<Geometry const> GetGeo() const;

};

