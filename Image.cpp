#include "Image.h"


Image::Image() {
	_img = new std::vector<Color>(0);
}
Image::Image(UINT _w, UINT _h) {
	w = _w;
	h = _h;
	_img = new std::vector<Color>(w*h);
}

bool Image::SetColor(UINT i, UINT j, const Color& c) {
	if (i >= w || j >= h) {  // Ô½½ç
		return false;
	}
	else {
		(*_img)[j*w+i] = c;
		return true;
	}
}

bool Image::Save(const std::string& path) {
	//FILE* f = fopen(path.c_str(), "w");         // Write image to PPM file. 
	//fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
	//for (int i = 0; i < w * h; i++)
	//	fprintf(f, "%d %d %d ", toInt((*_img)[i].x), toInt((*_img)[i].y), toInt((*_img)[i].z));


	unsigned char* image = new unsigned char[w * h * 4];

	for (size_t i = 0; i < w * h; i++) {
		image[4 * i + 0] = toInt((*_img)[i].x);
		image[4 * i + 1] = toInt((*_img)[i].y);
		image[4 * i + 2] = toInt((*_img)[i].z);
		image[4 * i + 3] = 0xFF;
	}

	lodepng_encode32_file(path.c_str(), image, w, h);

	return true;
}	