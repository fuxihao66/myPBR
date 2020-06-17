#pragma once
#include "utility.h"

struct Ray
{
	Vec3<double> ori;
	Vec3<double> dir;
	Ray();
	Ray(const Vec3<double>& ori, const Vec3<double>& dir);

	bool Valid() const;  // whether ray is valid
	void SetValidValue(bool v);
private:
	bool isValid;
};

