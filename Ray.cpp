#include "Ray.h"
Ray::Ray() {
	isValid = true;
}
Ray::Ray(const Vec3<double>& _ori, const Vec3<double>& _dir) {
	ori = _ori;
	dir = _dir;
	if (dir.length() > 1)
		dir.normalize();
	isValid = true;
}

bool Ray::Valid() const{  // whether ray is valid
	return isValid;
}

void Ray::SetValidValue(bool v) {
	isValid = v;
}