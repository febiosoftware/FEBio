#include "image_tools.h"
#include "Image.h"
#include <math.h>

#ifdef HAVE_MKL
#include <mkl.h>
#endif

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif // HAVE_FFTW

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
	// SL: Why not copy i.e. if (d <= 0) { trg = src; return; }?
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
	// SL: Why not copy i.e. if (d <= 0) { trg = src; return; }?
	if (d <= 0.f)
	{
		for (int i = 0; i < nx * ny * nz; ++i) y[i] = x[i];
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
void fftblur_2d(Image& trg, Image& src, float d) {}
void fftblur_3d(Image& trg, Image& src, float d) {}
#endif // HAVE_MKL

#ifdef HAVE_FFTW
void create_gaussian_2d(double* kernel, int rows, int cols, double sigma) {
	double sum = 0.0;
	for (int y = 0; y < rows; ++y) {
		int dy = (y <= rows / 2) ? y : (rows - y);  // wrap-around distance
		for (int x = 0; x < cols; ++x) {
			int dx = (x <= cols / 2) ? x : (cols - x);
			double d2 = (double)(dx * dx + dy * dy);
			double val = exp(-d2 / (2.0 * sigma * sigma));
			kernel[y * cols + x] = val;
			sum += val;
		}
	}
	// Normalize
	for (int i = 0; i < rows * cols; ++i) {
		kernel[i] /= sum;
	}
}

// Create 3D Gaussian kernel (centered)
void create_gaussian_3d(double* kernel, int nz, int ny, int nx, double sigma) {
	double sum = 0.0;
	for (int z = 0; z < nz; ++z) {
		int dz = (z <= nz / 2) ? z : (nz - z);  // wrap-around distance
		for (int y = 0; y < ny; ++y) {
			int dy = (y <= ny / 2) ? y : (ny - y);
			for (int x = 0; x < nx; ++x) {
				int dx = (x <= nx / 2) ? x : (nx - x);
				double d2 = (double)(dx * dx + dy * dy + dz * dz);
				double val = exp(-d2 / (2.0 * sigma * sigma));
				kernel[(z * ny + y) * nx + x] = val;
				sum += val;
			}
		}
	}
	// Normalize
	for (int i = 0; i < nx * ny * nz; ++i) {
		kernel[i] /= sum;
	}
}

void fftw_blur_2d(Image& trg, Image& src, float d)
{
	int width = src.width();
	int height = src.height();

	float* img_data = src.data();
	double sigma = d;

	// Allocate FFTW arrays
	double* img_spatial = (double*) fftw_malloc(sizeof(double) * width * height);
	double* kernel_spatial = (double*) fftw_malloc(sizeof(double) * width * height);
	fftw_complex* img_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * height * (width / 2 + 1));
	fftw_complex* kernel_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * height * (width / 2 + 1));

	// Copy image to double buffer
	for (int i = 0; i < width * height; ++i) {
		img_spatial[i] = (double)img_data[i];
	}

	// Create Gaussian kernel
	create_gaussian_2d(kernel_spatial, height, width, sigma);

	// Create FFT plans
	fftw_plan plan_fwd_img = fftw_plan_dft_r2c_2d(height, width, img_spatial, img_freq, FFTW_ESTIMATE);
	fftw_plan plan_fwd_kernel = fftw_plan_dft_r2c_2d(height, width, kernel_spatial, kernel_freq, FFTW_ESTIMATE);
	fftw_plan plan_inv = fftw_plan_dft_c2r_2d(height, width, img_freq, img_spatial, FFTW_ESTIMATE);

	// Forward FFTs
	fftw_execute(plan_fwd_img);
	fftw_execute(plan_fwd_kernel);

	// Multiply in frequency domain
	int nfreq = height * (width / 2 + 1);
	for (int i = 0; i < nfreq; ++i) {
		double a = img_freq[i][0], b = img_freq[i][1];
		double c = kernel_freq[i][0], d = kernel_freq[i][1];
		img_freq[i][0] = a * c - b * d;
		img_freq[i][1] = a * d + b * c;
	}

	// Inverse FFT
	fftw_execute(plan_inv);

	// Normalize and convert back to float
	float* out_data = trg.data();
	for (int i = 0; i < width * height; ++i) {
		double val = img_spatial[i] / (width * height);
		out_data[i] = (float)(val);
	}

	// Cleanup
	fftw_destroy_plan(plan_fwd_img);
	fftw_destroy_plan(plan_fwd_kernel);
	fftw_destroy_plan(plan_inv);
	fftw_free(img_spatial);
	fftw_free(kernel_spatial);
	fftw_free(img_freq);
	fftw_free(kernel_freq);
}

void fftw_blur_3d(Image& trg, Image& src, float d) 
{
	int nx = src.width();
	int ny = src.height();
	int nz = src.depth();

	double sigma = d;
	float* img_data = src.data();

	// Allocate FFTW buffers
	double* img_spatial = (double*)fftw_malloc(sizeof(double) * nx* ny*nz);
	double* kernel_spatial = (double*)fftw_malloc(sizeof(double) * nx * ny * nz);
	fftw_complex* img_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nz * ny * (nx / 2 + 1));
	fftw_complex* kernel_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nz * ny * (nx / 2 + 1));

	// Copy image to double buffer
	for (int i = 0; i < nx*ny*nz; ++i) {
		img_spatial[i] = (double)img_data[i];
	}

	// Create Gaussian kernel
	create_gaussian_3d(kernel_spatial, nz, ny, nx, sigma);

	// Plans
	fftw_plan plan_fwd_img = fftw_plan_dft_r2c_3d(nz, ny, nx, img_spatial, img_freq, FFTW_ESTIMATE);
	fftw_plan plan_fwd_kernel = fftw_plan_dft_r2c_3d(nz, ny, nx, kernel_spatial, kernel_freq, FFTW_ESTIMATE);
	fftw_plan plan_inv = fftw_plan_dft_c2r_3d(nz, ny, nx, img_freq, img_spatial, FFTW_ESTIMATE);

	// Forward FFTs
	fftw_execute(plan_fwd_img);
	fftw_execute(plan_fwd_kernel);

	// Multiply in frequency domain
	int nfreq = nz * ny * (nx / 2 + 1);
	for (int i = 0; i < nfreq; ++i) {
		double a = img_freq[i][0];
		double b = img_freq[i][1];
		double c = kernel_freq[i][0];
		double d = kernel_freq[i][1];
		img_freq[i][0] = a * c - b * d;
		img_freq[i][1] = a * d + b * c;
	}

	// Inverse FFT
	fftw_execute(plan_inv);

	// Normalize and convert back to float
	float* out_data = trg.data();
	double norm_factor = (double)(nx * ny * nz);
	for (int i = 0; i < nx * ny * nz; ++i) {
		double val = img_spatial[i] / norm_factor;
		out_data[i] = (float)(val);
	}

	// Cleanup
	fftw_destroy_plan(plan_fwd_img);
	fftw_destroy_plan(plan_fwd_kernel);
	fftw_destroy_plan(plan_inv);
	fftw_free(img_spatial);
	fftw_free(kernel_spatial);
	fftw_free(img_freq);
	fftw_free(kernel_freq);
}

#else
void fftw_blur_2d(Image& trg, Image& src, float d) {}
void fftw_blur_3d(Image& trg, Image& src, float d) {}
#endif // HAVE_FFTW
