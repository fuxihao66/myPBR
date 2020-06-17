#include "Phong.h"


Phong::Phong(MatType m, const Vec3<double>& _albedo, UINT n) 
	: Material(m, _albedo)
{
	_phongCoeff = n;
}

Ray Phong::ReflImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf, Color& albedo) const {

	albedo = GetAlbedo();

	return ReflImpSampling(r, pos, normal, pdf, brdf);


}

Ray Phong::RefrImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf, Color& albedo) const {
	return Ray();
}


Ray Phong::ReflImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf) const {

	Vec3<double> modifiedNorm = normal;
	if (normal.dot(r.dir) > 0)
		modifiedNorm = -1 * normal;


	Vec3<double> mirrorReflVec = sqrt(2) * r.dir + 2.0 * modifiedNorm;
	mirrorReflVec.normalize();



	Vec3<double> newDir = SamplePhongRefl(modifiedNorm, mirrorReflVec, _phongCoeff);
	Ray newRay(pos, newDir);

	// return by reference
	brdf = BRDF(r.dir, newRay.dir, modifiedNorm);
	pdf = brdf * (_phongCoeff + 1) / (_phongCoeff + 2);

	return newRay;


}

Ray Phong::RefrImpSampling(const Ray& r, const Vec3<double>& pos, const Vec3<double>& normal, double& pdf, double& brdf) const {
	return Ray();
}

// no transmission
double Phong::BTDF(const Vec3<double>& wi, const Vec3<double>& wo, const Vec3<double>& norm) const {
	return 0.0;
}

double Phong::BRDF(const Vec3<double>& wi, const Vec3<double>& wo, const Vec3<double>& norm) const {
	
	Vec3<double> modifiedNorm = norm;
	if (norm.dot(wi) > 0)
		modifiedNorm = -1 * norm;

	Vec3<double> mirrorReflVec = sqrt(2) * wi + 2.0 * modifiedNorm;
	mirrorReflVec.normalize();

	return  (_phongCoeff +2.0)/(2*PI)* pow(wo.dot(mirrorReflVec), _phongCoeff);
}

Color Phong::AnalyticLineIllumination(const std::vector<Vec2<double>>& lineSeg, const LightLineSample& lls, const Point& inscPoint, const Vec3<double>& inscNorm, const Vec3<double>& lightNorm, const Vec3<double>& incidentDir)const {
	// TODO:  norm 方向

	if (lineSeg.size() == 0)
		return Color(0.0, 0.0, 0.0);

	/*std::cout << "line sample num is " << lineSeg.size() << std::endl;
	std::cout << "t1 " << lineSeg[0].x << " t2: " << lineSeg[0].y << std::endl;*/
	//int k;
	//std::cin >> k;

	double tMin = lineSeg[lineSeg.size() - 1].x;
	double tMax = lineSeg[0].y;
	Point pMin = lls.L0 + tMin * lls.Ldir;
	Point pMax = lls.L0 + tMax * lls.Ldir;
	Vec3<double> dir1 = pMin - inscPoint; dir1.normalize();
	Vec3<double> dir2 = pMax - inscPoint; dir2.normalize();

	Vec3<double> ZAxis = dir2.OuterProd(dir1).normalize();
	Vec3<double> YAxis = dir1;
	Vec3<double> XAxis = YAxis.OuterProd(ZAxis).normalize();

	// 镜面反射的反射方向
	Vec3<double> mirrorReflectDir = sqrt(2) * incidentDir + 2.0 * inscNorm;
	mirrorReflectDir.normalize();

	double thetaL = atan2(lls.Ldir.dot(YAxis), lls.Ldir.dot(XAxis));
	double phiNx = atan2(inscNorm.dot(ZAxis), inscNorm.dot(XAxis));
	double phiNy = atan2(lightNorm.dot(ZAxis), lightNorm.dot(XAxis));
	Point newL0 = lls.L0 - inscPoint; // 变换到inscPoint为原点
	double L0x = newL0.dot(XAxis);
	double L0y = newL0.dot(YAxis);

	double commonFactor = abs( sin(phiNx) * sin(phiNy) / (L0x * sin(thetaL) - L0y * cos(thetaL)) );
	

	double thetaR = atan2(mirrorReflectDir.dot(XAxis), mirrorReflectDir.dot(YAxis));
	double thetaNx = atan2(inscNorm.dot(XAxis), inscNorm.dot(YAxis));
	double thetaNy = atan2(lightNorm.dot(XAxis), lightNorm.dot(YAxis));

	double integral = 0.0;

	// 假设seg是降序排序  seg.x < seg.y
	for (auto seg : lineSeg) {
		Vec3<double> dirMin = newL0 + seg.x * lls.Ldir; dirMin.normalize();
		Vec3<double> dirMax = newL0 + seg.y * lls.Ldir; dirMax.normalize();

		double thetaMin = acos(ClampBeforeACosAndASin(dirMin.dot(dir1)));
		double thetaMax = acos(ClampBeforeACosAndASin(dirMax.dot(dir1)));

		double uMin = thetaMin - thetaR;
		double uMax = thetaMax - thetaR;
		// TODO: abs?????
		integral += cos(thetaR - thetaNx) * cos(thetaR - thetaNy) * IntegrateCosExp(_phongCoeff + 2, uMin, uMax) ;
		integral -= cos(thetaR - thetaNx) * sin(thetaR - thetaNy) * IntegrateCosExpSin(_phongCoeff + 1, uMin, uMax);
		integral -= sin(thetaR - thetaNx) * cos(thetaR - thetaNy) * IntegrateCosExpSin(_phongCoeff + 1, uMin, uMax);
		integral += sin(thetaR - thetaNx) * sin(thetaR - thetaNy) * (IntegrateCosExp(_phongCoeff, uMin, uMax) - IntegrateCosExp(_phongCoeff + 2, uMin, uMax));

		

	}

	return commonFactor * (_phongCoeff + 2.0) / (2.0 * PI) * abs(integral) * albedo * lls.illumination;
}

