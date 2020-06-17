#include "SceneManager.h"


SceneManager::SceneManager() {

}
void SceneManager::SetCamera(std::shared_ptr<Camera const > cam) {
	_pCam = cam;
}

bool SceneManager::BuildDefaultScene() {

	std::shared_ptr<Material> matt1Mat = std::make_shared<PerfectDiffuse>(MatType::Reflect, Vec3<double>(0.75f,0.75f,0.75f));
	std::shared_ptr<Material> matt2Mat = std::make_shared<PerfectDiffuse>(MatType::Reflect, Vec3<double>(0.75f, 0.25f, 0.25f));
	std::shared_ptr<Material> matt3Mat = std::make_shared<PerfectDiffuse>(MatType::Reflect, Vec3<double>(0.25f, 0.25f, 0.75f));
	std::shared_ptr<Material> mattMat = std::make_shared<PerfectDiffuse>(MatType::Reflect, Vec3<double>(0.1f, 0.1f, 0.1f));
	std::shared_ptr<Material> matt4Mat = std::make_shared<PerfectDiffuse>(MatType::Reflect, Vec3<double>(0.25f, 0.75f, 0.25f));

	std::shared_ptr<Material> phongMat = std::make_shared<Phong>(MatType::Reflect, Vec3<double>(0.9f, 0.9f, 0.9f), 50);

	auto planeLeft = std::make_shared<RenderableObject>(matt2Mat, std::make_shared<Plane>(Vec3<double>(-5,-5,-5), Vec3<double>(-5,5,-5), Vec3<double>(-5, -5, 5)));
	auto planeRight = std::make_shared<RenderableObject>(matt3Mat, std::make_shared<Plane>(Vec3<double>(5, -5, -5), Vec3<double>(5, 5, -5), Vec3<double>(5, -5, 5)));

	auto planeFar = std::make_shared<RenderableObject>(phongMat, std::make_shared<Plane>(Vec3<double>(-5,-5,-5), Vec3<double>(5,-5,-5), Vec3<double>(5,5,-5)));

	auto planeNear = std::make_shared<RenderableObject>(mattMat, std::make_shared<Plane>(Vec3<double>(-5, -5, 5), Vec3<double>(5, -5, 5), Vec3<double>(5, 5, 5)));

	auto planeTop = std::make_shared<RenderableObject>(matt1Mat, std::make_shared<Plane>(Vec3<double>(-5,5,-5), Vec3<double>(5,5,-5), Vec3<double>(5,5,5)));

	auto planeBottom = std::make_shared<RenderableObject>(matt4Mat, std::make_shared<Plane>(Vec3<double>(5,-5,5), Vec3<double>(5,-5,-5), Vec3<double>(-5,-5,-5)));

	auto emissionMat = std::make_shared<Material>(MatType::Emission, Vec3<double>(18.0f, 18.0f, 18.0f));
	auto emissionMat2 = std::make_shared<Material>(MatType::Emission, Vec3<double>(9.0f, 10.0f, 10.0f));
	
	auto mainLight = std::make_shared<Light>(emissionMat, std::make_shared<Plane>(Vec3<double>(-2, 4.9, 2), Vec3<double>(-2, 4.9, -2),Vec3<double>(2, 4.9, -2))  );
	/*auto mainLight2 = std::make_shared<Light>(emissionMat2, std::make_shared<Plane>(Vec3<double>(4.9999, 5, -5),Vec3<double>(4.9999, -5, -5), Vec3<double>(4.9999, -5, 5)));
	auto mainLight3 = std::make_shared<Light>(emissionMat, std::make_shared<Plane>(Vec3<double>(-4.9999, 5, -5), Vec3<double>(-4.9999, -5, -5), Vec3<double>(-4.9999, -5, 5)));
	auto backLight = std::make_shared<Light>(emissionMat, std::make_shared<Plane>(Vec3<double>(-5, -5, -4.99), Vec3<double>(5, -5, -4.99), Vec3<double>(5, 5, -4.99)));*/


	_objList.AddObject(planeLeft);
	_objList.AddObject(planeRight);
	_objList.AddObject(planeFar);
	_objList.AddObject(planeNear);
	_objList.AddObject(planeTop);
	_objList.AddObject(planeBottom);
	_objList.AddObject(mainLight);
	
	_lightList.push_back(mainLight);

	return true;
}

void SceneManager::AddObject(std::shared_ptr<RenderableObject const> obj) {
	_objList.AddObject(obj);
}

// only used for line sampling
void SceneManager::AddTriangle(std::shared_ptr<Triangle const> obj) {
	_triList.push_back(obj);
}

void SceneManager::SetConfig(const Config& f) {
	_f = f;
}

void SceneManager::DoRenderHeightRange(const UINT heightStart, const UINT heightEnd, const UINT width, const UINT SSP, const Point& camPos,  double pixelWidth, double pixelHeight, std::shared_ptr<Image> img, const SceneManager* psm)  {
	for (UINT col = heightStart; col < heightEnd; col++) {
		for (UINT row = 0; row < width; row++) {
			// subpixel 
			Color finalColor;
			for (int sy = 0; sy < 2; sy++)
				for (int sx = 0; sx < 2; sx++) {
					Color color;
					for (UINT sampIndex = 0; sampIndex < SSP; sampIndex++) {

						//stratified 
						double randomWidthOffset = SampleRandomNum(0, pixelWidth / 2) + pixelWidth / 2 * sx;
						double randomHeightOffset = SampleRandomNum(0, pixelHeight / 2) + pixelHeight / 2 * sy;
						/*Vec3<double> randomPoint(-5.0f + row * widthPixel + widthPixel / 2, -(-5.0f + col * heightPixel + heightPixel / 2), 5.0f);*/
						Vec3<double> randomPoint(-5.0f + row * pixelWidth + randomWidthOffset, -(-5.0f + col * pixelHeight + randomHeightOffset), 4.9f);

						auto dir = randomPoint - camPos;
						dir.normalize();
						Ray r(randomPoint, dir);

						//color = color + DoRenderPathTracing(r, 1);
						color = color + psm->DoRenderPathTracingWithDirectIllumination(r, 1);
					}
					color = (1.0 / SSP) * color;

					finalColor = finalColor + 0.25 * Vec3<double>(clamp(color.x), clamp(color.y), clamp(color.z));
				}
			img->SetColor(row, col, finalColor);
			//count++;
			//std::cout << 100 * double(count) / _f.HEIGHT / _f.WIDTH << "\%" << " over" << std::endl;
		}
	}
	std::cout << "thread exit" << std::endl;
}

std::shared_ptr<Image> SceneManager::DoRender() const {
	auto img = std::make_shared<Image>(_f.WIDTH, _f.HEIGHT);

	//auto eachThreadRow = ceil((double)(_f.WIDTH) / _f.THREAD);

	// plane is on z=5
	// TODO: camera module
	double widthPixel = 10 / (double)_f.WIDTH;
	double heightPixel = 10 / (double)_f.HEIGHT;
	//size_t count = 0;

	std::vector<std::thread> threadGroup(_f.THREAD);

	UINT heightGap = ceil(double(_f.HEIGHT) / _f.THREAD);

	for (UINT tIndex = 0; tIndex < threadGroup.size(); tIndex++) {
		threadGroup[tIndex] = std::thread(SceneManager::DoRenderHeightRange, tIndex*heightGap, std::min(tIndex * heightGap+heightGap, _f.HEIGHT), _f.WIDTH, _f.SSP, _pCam->position, widthPixel, heightPixel ,img, this);
	}
	
	for (size_t tIndex = 0; tIndex < threadGroup.size(); tIndex++) {
		threadGroup[tIndex].join();
	}

	//for (UINT col = 0*_f.HEIGHT; col < _f.HEIGHT; col++) {
	//	for (UINT row = 0; row < _f.WIDTH; row++) {
	//		// subpixel 
	//		Color finalColor;
	//		for (int sy = 0; sy < 2; sy++)      
	//			for (int sx = 0; sx < 2; sx++) {      
	//				Color color;
	//				for (UINT sampIndex = 0; sampIndex < _f.SSP; sampIndex++) {

	//					//stratified 
	//					double randomWidthOffset = SampleRandomNum(0, widthPixel / 2) + widthPixel / 2 * sx;
	//					double randomHeightOffset = SampleRandomNum(0, heightPixel / 2) + heightPixel / 2 * sy;
	//					/*Vec3<double> randomPoint(-5.0f + row * widthPixel + widthPixel / 2, -(-5.0f + col * heightPixel + heightPixel / 2), 5.0f);*/
	//					Vec3<double> randomPoint(-5.0f + row * widthPixel + randomWidthOffset, -(-5.0f + col * heightPixel + randomHeightOffset), 4.9f);

	//					auto dir = randomPoint - _pCam->position;
	//					dir.normalize();
	//					Ray r(randomPoint, dir);

	//					//color = color + DoRenderPathTracing(r, 1);
	//					color = color + DoRenderPathTracingWithDirectIllumination(r, 1);
	//				}
	//				color = (1.0 / _f.SSP) * color;

	//				finalColor = finalColor + 0.25 * Vec3<double>(clamp(color.x), clamp(color.y), clamp(color.z));
	//			}
	//		img->SetColor(row, col, finalColor);
	//		count++;
	//		std::cout <<  100*double(count) / _f.HEIGHT / _f.WIDTH << "\%" << " over" << std::endl;
	//	}
	//}

	return img;
}




Color SceneManager::DoRenderPathTracing(const Ray& r, int depth) const {
	if (!r.Valid())
		return Color(0.0, 0.0, 0.0);

	Vec3<double>	inscNormal;
	Point			inscPosition;

	// 获取相交点的信息: 相交物体 normal inscPoint
	auto intersecObj = _objList.Intersect(r, inscNormal, inscPosition);

	if (intersecObj == nullptr)  // 没有相交
		return Color(0.0, 0.0, 0.0);

	auto intersecMat = intersecObj->GetMat();

	Color	albedo		= intersecMat->GetAlbedo();
	Color	emission	= intersecMat->GetEmission();
	auto	matType		= intersecMat->GetMatType();

	//RR
	double p = 0.6;
	double rnd = SampleRandomNum(0, 1);
	if (rnd > p)
		return emission;   // 终止trace

	albedo = (1 / p) * albedo;


	double reflPDF;
	double refrPDF;
	double brdf;
	double btdf;
	Ray reflRay;
	Ray refrRay;

	switch (matType) {
	case MatType::Emission:

		return emission;
		break;
	case MatType::Reflect:

		reflRay = intersecMat->ReflImpSampling(r, inscPosition, inscNormal, reflPDF, brdf);
		return emission + abs(inscNormal.dot(reflRay.dir) / reflPDF) * brdf * albedo * DoRenderPathTracing(reflRay, depth + 1);
		break;
	case MatType::Transmit:
		// 处理折射 （根据法线和光线的夹角判断在里在外）
		refrRay = intersecMat->RefrImpSampling(r, inscPosition, inscNormal, refrPDF, btdf);

		return emission + abs(inscNormal.dot(refrRay.dir) / refrPDF) * btdf * albedo * DoRenderPathTracing(refrRay, depth + 1);
		break;


	//TODO: 待修改  
	//case MatType::Both:
	//	// 处理折射 （根据法线和光线的夹角判断在里在外）
	//	// 解frenel方程
	//	Ray reflRay = intersecMat->GetRefl(r, inscPosition, inscNormal, reflPDF);
	//	Ray refrRay = intersecMat->GetRefr(r, inscPosition, inscNormal, refrPDF);

	//	double reflRatio = 0.5;
	//	double refrRatio = 0.5;

	//	return (emission + reflRatio * albedo * DoRenderPathTracing(reflRay, depth + 1) + refrRatio * albedo * DoRenderPathTracing(refrRay, depth + 1));
	//	break;
	//TODO: 待增加  次表面散射
	//case MatType::Subsurface:
	}

}


Color SceneManager::DirectIlluminationPoint(const Ray& r, const Point& inscPoint, const Vec3<double>& inscNormal, std::shared_ptr<Material const> mat, size_t sampleNum) const{
	Color directLight(0.0, 0.0, 0.0);
	for (auto pLight : _lightList) {
		Color perLight(0.0, 0.0, 0.0);
		for (size_t index = 0; index < sampleNum; index++) {
			double tMax, t;
			Vec3<double> lightNorm;
			Ray shadowRay = pLight->SamplePoint(inscPoint, tMax); 
			auto insecIndicator = _objList.Intersect(shadowRay, t);
			if (t > EPSILON && t < tMax-EPSILON){ 
				// in shadow
				continue;
			}
			pLight->Intersect(shadowRay, t, lightNorm);

			perLight = perLight +  mat->BRDF(r.dir, shadowRay.dir, inscNormal) * pLight->FuncG(shadowRay.dir.dot(inscNormal), shadowRay.dir.dot(lightNorm), tMax) / pLight->PointSamplePDF() * pLight->GetMat()->GetEmission() * mat->GetAlbedo();
		}
		directLight = directLight + (1.0 / sampleNum) * perLight;
	}

	return Color(clamp(directLight.x), clamp(directLight.y), clamp(directLight.z));
}

Color SceneManager::DirectIlluminationLine(const Ray& r, const Point& inscPoint, const Vec3<double>& inscNormal, std::shared_ptr<Material const> mat) const{
	 /*norm 方向
	 光线norm方向*/
	Vec3<double> modifiedNormal = inscNormal;
	if (inscNormal.dot(r.dir) > 0){
		modifiedNormal = -1*inscNormal;
	}
	Color directLight(0.0, 0.0, 0.0);

	for (auto pLight : _lightList) {
		double pdf;
		LightLineSample lls = pLight->SampleLine(pdf);

		Vec3<double> reflDir = lls.L0 + lls.length * lls.Ldir;
		Vec3<double> lightNorm = lls.normal;
		if (lls.normal.dot(reflDir) < 0) { // TODO:
			lightNorm = -1 * lls.normal;
		}

		
		
		lls.VisibilityTest(_triList, inscPoint);
		auto lineSegs = lls.LineSegments();

		/*if (lineSegs.size() > 0)
			std::cout << "line segs is " << lineSegs.size() << std::endl;*/
		
		directLight = directLight + (1.0 / pdf) * mat->AnalyticLineIllumination(lineSegs, lls, inscPoint, modifiedNormal, lightNorm, r.dir);
	}
	
	return directLight;
}

Color SceneManager::DoRenderPathTracingWithDirectIllumination(const Ray& r, int depth) const {
	if (!r.Valid())
		return Color(0.0, 0.0, 0.0);

	Vec3<double>	inscNormal;
	Point			inscPosition;

	// 获取相交点的信息: 相交物体 normal inscPoint
	auto intersecObj = _objList.Intersect(r, inscNormal, inscPosition);

	if (intersecObj == nullptr)  // 没有相交
		return Color(0.0, 0.0, 0.0);

	auto intersecMat = intersecObj->GetMat();

	Color	albedo = intersecMat->GetAlbedo();
	Color	emission = intersecMat->GetEmission();
	auto	matType = intersecMat->GetMatType();


	// 10 samples on light
	Color directLighting;
	
	
	//directLighting = DirectIlluminationPoint(r, inscPosition, inscNormal, intersecMat, 20);
	directLighting = DirectIlluminationLine(r, inscPosition, inscNormal, intersecMat);


	// only direct illumination
	if (matType == MatType::Emission) {
		return emission;
	}
	else
		return directLighting;

	
	//RR
	double p = 0.6 / sqrt(double(depth));


	double rnd = SampleRandomNum(0, 1);
	if (depth > 5 && rnd > p)
		return directLighting;   // 终止trace
	else if (depth > 5) // TODO: 
		albedo = (1 / p) * albedo;


	double reflPDF;
	double refrPDF;
	double brdf;
	double btdf;
	Ray reflRay;
	Ray refrRay;


	switch (matType) {
	case MatType::Emission:
		// 如果直接看到光源  则返回光源emission
		if (depth == 1)
			return emission;
		// 否则不考虑光源发光 防止计算两次
		else
			return Color(0.0, 0.0, 0.0);
		break;
	case MatType::Reflect:

		reflRay = intersecMat->ReflImpSampling(r, inscPosition, inscNormal, reflPDF, brdf);
		return directLighting + abs(inscNormal.dot(reflRay.dir) / reflPDF) * brdf * albedo * DoRenderPathTracingWithDirectIllumination(reflRay, depth + 1);
		break;
	case MatType::Transmit:
		// 处理折射 （根据法线和光线的夹角判断在里在外）
		refrRay = intersecMat->RefrImpSampling(r, inscPosition, inscNormal, refrPDF, btdf);

		return directLighting + abs(inscNormal.dot(refrRay.dir) / refrPDF) * btdf * albedo * DoRenderPathTracingWithDirectIllumination(refrRay, depth + 1);
		break;
	}

}

