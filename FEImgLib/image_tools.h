#pragma once
#include <FECore/FECoreClass.h>
#include "feimglib_api.h"

class Image;

// iterative stencil-based blur
enum class BlurMethod
{
	BLUR_AVERAGE,
	BLUR_FFT,
};

FEIMGLIB_API void blur_image_2d(Image& trg, Image& src, float d, BlurMethod blurMethod = BlurMethod::BLUR_AVERAGE);
FEIMGLIB_API void blur_image_3d(Image& trg, Image& src, float d, BlurMethod blurMethod = BlurMethod::BLUR_AVERAGE);


// FFT-based blurs (rely on MKL library. Incompatible with pardiso/dss)
FEIMGLIB_API void fftblur_2d(Image& trg, Image& src, float d);
FEIMGLIB_API void fftblur_3d(Image& trg, Image& src, float d);

// FFT-based blurs (rely on FFTW library)
FEIMGLIB_API void fftw_blur_2d(Image& trg, Image& src, float d[2]);
FEIMGLIB_API void fftw_blur_3d(Image& trg, Image& src, float d[3]);

// Extension functions for image filter edges
FEIMGLIB_API int constant_extension(int N, int n);
FEIMGLIB_API int symmetric_extension(int N, int n);
