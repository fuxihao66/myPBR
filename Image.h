#pragma once

#include <string>
#include <vector>
#include "utility.h"
#include "lodepng/lodepng.h"

class Image
{
private:
	// ´æ´¢½á¹¹
	std::vector<Color>* _img;
	UINT w, h;
public:
	Image();
	Image(UINT w, UINT h);

	bool SetColor(UINT i, UINT j, const Color& c);

	bool Save(const std::string& path);

};

