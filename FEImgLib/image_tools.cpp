#include "image_tools.h"
#include "Image.h"
#include <math.h>

#ifdef HAVE_MKL
#include <mkl.h>
#endif

//-----------------------------------------------------------------------------
void blur_image_2d(Image& trg, Image& src, float d)
{
	if (d <= 0) { trg = src; return; }

	int n = (int)d;
	float w = d - (float)n;

	int nx = src.width();
	int ny = src.height();
	int nz = src.depth();

	trg = src;
	Image tmp(src);
	float f[4];
	for (int l = 0; l < n; ++l)
	{
		for (int k = 0; k < nz; ++k)
			for (int j = 0; j < ny; ++j)
				for (int i = 0; i < nx; ++i)
				{
					if (i > 0) f[0] = tmp.value(i - 1, j, k); else f[0] = tmp.value(i, j, k);
					if (i < nx - 1) f[1] = tmp.value(i + 1, j, k); else f[1] = tmp.value(i, j, k);
					if (j > 0) f[2] = tmp.value(i, j - 1, k); else f[2] = tmp.value(i, j, k);
					if (j < ny - 1) f[3] = tmp.value(i, j + 1, k); else f[3] = tmp.value(i, j, k);
					trg.value(i, j, k) = 0.25f * (f[0] + f[1] + f[2] + f[3]);
				}
		tmp = trg;
	}

	if (w > 0.0)
	{
		for (int k = 0; k < nz; ++k)
			for (int j = 0; j < ny; ++j)
				for (int i = 0; i < nx; ++i)
				{
					if (i > 0) f[0] = tmp.value(i - 1, j, k); else f[0] = tmp.value(i, j, k);
					if (i < nx - 1) f[1] = tmp.value(i + 1, j, k); else f[1] = tmp.value(i, j, k);
					if (j > 0) f[2] = tmp.value(i, j - 1, k); else f[2] = tmp.value(i, j, k);
					if (j < ny - 1) f[3] = tmp.value(i, j + 1, k); else f[3] = tmp.value(i, j, k);
					float f1 = 0.25f * (f[0] + f[1] + f[2] + f[3]);
					float f2 = trg.value(i, j, k);
					trg.value(i, j, k) = f1 * w + f2 * (1.f - w);
				}
	}
}

//-----------------------------------------------------------------------------
void blur_image(Image& trg, Image& src, float d)
{
	if (d <= 0) { trg = src; return; }

	int n = (int)d;
	float w = d - (float)n;

	int nx = src.width();
	int ny = src.height();
	int nz = src.depth();

	trg = src;
	Image tmp(src);
	float f[6];
	for (int l = 0; l < n; ++l)
	{
		for (int k = 0; k < nz; ++k)
			for (int j = 0; j < ny; ++j)
				for (int i = 0; i < nx; ++i)
				{
					if (i > 0) f[0] = tmp.value(i - 1, j, k); else f[0] = tmp.value(i, j, k);
					if (i < nx - 1) f[1] = tmp.value(i + 1, j, k); else f[1] = tmp.value(i, j, k);
					if (j > 0) f[2] = tmp.value(i, j - 1, k); else f[2] = tmp.value(i, j, k);
					if (j < ny - 1) f[3] = tmp.value(i, j + 1, k); else f[3] = tmp.value(i, j, k);
					if (k > 0) f[4] = tmp.value(i, j, k - 1); else f[4] = tmp.value(i, j, k);
					if (k < nz - 1) f[5] = tmp.value(i, j, k + 1); else f[5] = tmp.value(i, j, k);
					trg.value(i, j, k) = 0.1666667f * (f[0] + f[1] + f[2] + f[3] + f[4] + f[5]);
				}
		tmp = trg;
	}

	if (w > 0.0)
	{
		for (int k = 0; k < nz; ++k)
			for (int j = 0; j < ny; ++j)
				for (int i = 0; i < nx; ++i)
				{
					if (i > 0) f[0] = tmp.value(i - 1, j, k); else f[0] = tmp.value(i, j, k);
					if (i < nx - 1) f[1] = tmp.value(i + 1, j, k); else f[1] = tmp.value(i, j, k);
					if (j > 0) f[2] = tmp.value(i, j - 1, k); else f[2] = tmp.value(i, j, k);
					if (j < ny - 1) f[3] = tmp.value(i, j + 1, k); else f[3] = tmp.value(i, j, k);
					if (k > 0) f[4] = tmp.value(i, j, k - 1); else f[4] = tmp.value(i, j, k);
					if (k < nz - 1) f[5] = tmp.value(i, j, k + 1); else f[5] = tmp.value(i, j, k);
					float f1 = 0.1666667f * (f[0] + f[1] + f[2] + f[3] + f[4] + f[5]);
					float f2 = trg.value(i, j, k);
					trg.value(i, j, k) = f1 * w + f2 * (1.f - w);
				}
	}
}

#ifdef HAVE_MKL

// in fft.cpp
bool mkl_dft2(int nx, int ny, float* x, MKL_Complex8* c);
bool mkl_idft2(int nx, int ny, MKL_Complex8* c, float* y);

FEIMGLIB_API void fftblur_2d(Image& trg, Image& src, float d)
{
	int nx = src.width();
	int ny = src.height();

	float* x = src.data();
	float* y = trg.data();

	// for zero blur radius, we just copy the image
	if (d <= 0.f)
	{
		for (int i = 0; i < nx * ny; ++i) y[i] = x[i];
		return;
	}

	// since the blurring is done in Fourier space,
	// we need to invert the blur radius
	float sigmax = nx / d;
	float sigmay = ny / d;

	// calculate the DFT
	MKL_Complex8* c = new MKL_Complex8[nx * ny];
	mkl_dft2(nx, ny, x, c);

	// multiply the DFT with blur mask
	for (int j = 0; j <= ny/2; ++j)
		for (int i = 0; i < nx; ++i)
		{
			double wx = (i < nx / 2 ? i : i - nx) / sigmax;
			double wy = j / sigmay;
			float v = (float) exp(-(wx * wx + wy * wy));

			c[j * nx + i].real *= v;
			c[j * nx + i].imag *= v;
		}

	// calculate the inverse DFT
	mkl_idft2(nx, ny, c, y);

	// clean up
	delete[] c;
}

// in fft.cpp
bool mkl_dft3(int nx, int ny, int nz, float* x, MKL_Complex8* c);
bool mkl_idft3(int nx, int ny, int nz, MKL_Complex8* c, float* y);
FEIMGLIB_API void fftblur_3d(Image& trg, Image& src, float d)
{
	int nx = src.width();
	int ny = src.height();
	int nz = src.depth();

	float* x = src.data();
	float* y = trg.data();

	// for zero blur radius, we just copy the image
	if (d <= 0.f)
	{
		for (int i = 0; i < nx * ny*nz; ++i) y[i] = x[i];
		return;
	}

	// since the blurring is done in Fourier space,
	// we need to invert the blur radius
	float sigmax = nx / d;
	float sigmay = ny / d;
	float sigmaz = nz / d;

	// calculate the DFT
	MKL_Complex8* c = new MKL_Complex8[nx * ny * nz];
	mkl_dft3(nx, ny, nz, x, c);

	// multiply the DFT with blur mask
	for (int k = 0; k <= nz / 2; ++k)
		for (int j = 0; j < ny; ++j)
			for (int i = 0; i < nx; ++i)
			{
				double wx = (i < nx / 2 ? i : i - nx) / sigmax;
				double wy = (j < ny / 2 ? j : j - ny) / sigmay;
				double wz = k / sigmaz;
				float v = (float)exp(-(wx * wx + wy * wy + wz * wz));

				c[k*nx*ny + j * nx + i].real *= v;
				c[k*nx*ny + j * nx + i].imag *= v;
			}

	// calculate the inverse DFT
	mkl_idft3(nx, ny, nz, c, y);

	// clean up
	delete[] c;
}
#else // HAVE_MKL
FEIMGLIB_API void fftblur_2d(Image& trg, Image& src, float d) {}
FEIMGLIB_API void fftblur_3d(Image& trg, Image& src, float d) {}
#endif // HAVE_MKL
