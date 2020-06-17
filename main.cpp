#include <memory>

#include "Camera.h"
#include "Material.h"
#include "SceneManager.h"
#include "Sphere.h"
//#include "Glass.h"
//#include "Mirror.h"
#include "Phong.h"
#include "ObjWarehouse.h"
#include "Mesh.h"
#include "TriangleGeo.h"



//int main1() {
//	// test intersection
//	ObjWarehouse _objList;
//	Ray r(Vec3<double>(-5, 4.99962, -3.51913), Vec3<double>(0.616031, 0.773086, -0.151143));
//
//	auto matt1Mat = std::make_shared<Material>(MatType::Matt, Vec3<double>(0.2f, 0.4f, 0.2f));
//	auto matt2Mat = std::make_shared<Material>(MatType::Matt, Vec3<double>(0.4f, 0.2f, 0.2f));
//	auto matt3Mat = std::make_shared<Material>(MatType::Matt, Vec3<double>(0.2f, 0.2f, 0.4f));
//	auto mattMat = std::make_shared<Material>(MatType::Matt, Vec3<double>(0.1f, 0.1f, 0.1f));
//
//	auto planeLeft = std::make_shared<RenderableObject>(matt1Mat, std::make_shared<Plane>(Vec3<double>(-5, -5, -5), Vec3<double>(-5, 5, -5), Vec3<double>(-5, -5, 5)));
//	auto planeRight = std::make_shared<RenderableObject>(matt1Mat, std::make_shared<Plane>(Vec3<double>(5, -5, -5), Vec3<double>(5, 5, -5), Vec3<double>(5, -5, 5)));
//	auto planeFar = std::make_shared<RenderableObject>(matt2Mat, std::make_shared<Plane>(Vec3<double>(-5, -5, -5), Vec3<double>(5, -5, -5), Vec3<double>(5, 5, -5)));
//
//	auto planeNear = std::make_shared<RenderableObject>(mattMat, std::make_shared<Plane>(Vec3<double>(-5, -5, 5), Vec3<double>(5, -5, 5), Vec3<double>(5, 5, 5)));
//
//	auto planeTop = std::make_shared<RenderableObject>(matt3Mat, std::make_shared<Plane>(Vec3<double>(-5, 5, -5), Vec3<double>(5, 5, -5), Vec3<double>(5, 5, 5)));
//	auto planeBottom = std::make_shared<RenderableObject>(matt3Mat, std::make_shared<Plane>(Vec3<double>(5, -5, 5), Vec3<double>(5, -5, -5), Vec3<double>(-5, -5, -5)));
//
//	_objList.AddObject(planeLeft);
//	_objList.AddObject(planeRight);
//	_objList.AddObject(planeFar);
//	_objList.AddObject(planeNear);
//	_objList.AddObject(planeTop);
//	_objList.AddObject(planeBottom);
//
//	_objList.Intersect(r, 1);
//
//	return 0;
//}


int main() {
	resetRandomGen();

	auto mainCam = std::make_shared<Camera>(Vec3<double>(0, 0, 10));
	auto sm = std::make_shared<SceneManager>();
	
	
	sm->BuildDefaultScene();  
	/*auto glassMat = std::make_shared<Material>(MatType::Glass, Color(0.9f, 0.9f, 0.9f));
	auto glassBall = std::make_shared<RenderableObject>();
	glassBall->SetMat(glassMat);
	glassBall->SetGeo(std::make_shared<Sphere>(0.0f, -3.0f, 0.0f, 2.0f));
	sm->AddObject(glassBall);*/



	//  可以增加折射
	std::shared_ptr<Material> mirror = std::make_shared<Phong>(MatType::Reflect, Color(0.9f, 0.9f, 0.9f), 10);  
	std::shared_ptr<Material> matt0 = std::make_shared<PerfectDiffuse>(MatType::Reflect, Color(0.3f, 0.5f, 0.7f));
	std::shared_ptr<Material> matt1 = std::make_shared<PerfectDiffuse>(MatType::Reflect, Color(0.7f, 0.5f, 0.3f));

	auto triangleGeo = std::make_shared<TriangleGeo>(Vec3<double>(5, 0, 0), Vec3<double>(-5, -5, 5), Vec3<double>(-5, 5, -5));

	auto triangle = std::make_shared<RenderableObject>(matt1, triangleGeo);

	auto triangleBlock = triangleGeo->ToTriangle();

	/*sm->AddObject(triangle);
	sm->AddTriangle(triangleBlock);*/

	/*auto materialTestBall0 = std::make_shared<RenderableObject>();
	materialTestBall0->SetMat(matt0);
	materialTestBall0->SetGeo(std::make_shared<Mesh>("models/Mesh000.ply"));

	auto materialTestBall1 = std::make_shared<RenderableObject>();
	materialTestBall1->SetMat(matt1);
	materialTestBall1->SetGeo(std::make_shared<Mesh>("models/Mesh001.ply"));

	auto materialTestBall2 = std::make_shared<RenderableObject>();
	materialTestBall2->SetMat(mirror);
	materialTestBall2->SetGeo(std::make_shared<Mesh>("models/Mesh002.ply"));*/

	/*auto teapot00 = std::make_shared<RenderableObject>();
	teapot00->SetMat(matt0);
	teapot00->SetGeo(std::make_shared<Mesh>("teapot/Mesh000.ply"));

	auto teapot01 = std::make_shared<RenderableObject>();
	teapot01->SetMat(matt1);
	teapot01->SetGeo(std::make_shared<Mesh>("teapot/Mesh001.ply"));*/

	std::cout << "after loaded" << std::endl;

	/*sm->AddObject(materialTestBall0);
	sm->AddObject(materialTestBall1);
	sm->AddObject(materialTestBall2);*/
	/*sm->AddObject(teapot00);
	sm->AddObject(teapot01);*/
	std::cout << "add to sm" << std::endl;

	// resolution x, resolution y, samples, threads
	Config f = Config(300,300, 20, 4);
	sm->SetConfig(f);
	sm->SetCamera(mainCam);

	std::cout << "start to render";
	auto image = sm->DoRender();

	if (image != nullptr) {
		image->Save("test.png");
	}

	return 0;
}