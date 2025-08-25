#pragma once
#include <FECore/FECoreClass.h>
#include "feimglib_api.h"

class Image;

FEIMGLIB_API void blur_image_2d(Image& trg, Image& src, float d);
FEIMGLIB_API void blur_image(Image& trg, Image& src, float d);

FEIMGLIB_API void fftblur_2d(Image& trg, Image& src, float d);
FEIMGLIB_API void fftblur_3d(Image& trg, Image& src, float d);

FEIMGLIB_API void fftw_blur_2d(Image& trg, Image& src, float d);
FEIMGLIB_API void fftw_blur_3d(Image& trg, Image& src, float d);
