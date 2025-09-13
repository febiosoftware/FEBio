#pragma once
#include "feimglib_api.h"

class Image;

enum class BlurMethod
{
	BLUR_AVERAGE,
	BLUR_FFT,
};

FEIMGLIB_API void blur_image_2d(Image& trg, Image& src, float d, BlurMethod blurMethod = BlurMethod::BLUR_AVERAGE);
FEIMGLIB_API void blur_image_3d(Image& trg, Image& src, float d, BlurMethod blurMethod = BlurMethod::BLUR_AVERAGE);
