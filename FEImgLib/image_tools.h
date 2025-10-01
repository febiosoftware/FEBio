#pragma once
#include <FECore/FECoreClass.h>
#include "feimglib_api.h"

class Image;

FEIMGLIB_API int constant_extension(int N, int n);
FEIMGLIB_API int symmetric_extension(int N, int n);

// iterative stencil-based blur
FEIMGLIB_API void blur_image_2d(Image& trg, Image& src, float d);
FEIMGLIB_API void blur_image(Image& trg, Image& src, float d);

// FFT-based blurs (rely on MKL library. Incompatible with pardiso/dss)
FEIMGLIB_API void fftblur_2d(Image& trg, Image& src, float d);
FEIMGLIB_API void fftblur_3d(Image& trg, Image& src, float d);

// FFT-based blurs (rely on FFTW library)
FEIMGLIB_API void fftw_blur_2d(Image& trg, Image& src, float d[2]);
FEIMGLIB_API void fftw_blur_3d(Image& trg, Image& src, double d[3]);
