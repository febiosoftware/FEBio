#include "Image.h"
#include "ImageFilter.h"
#include <FECore/log.h>
#include <vector>
#include "image_tools.h"

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif // HAVE_FFTW

ImageFilter::ImageFilter(FEModel* fem) : FECoreClass(fem)
{

}

//=================================================================================================
// IterativeBlur
//=================================================================================================

BEGIN_FECORE_CLASS(IterativeBlur, ImageFilter)
	ADD_PARAMETER(m_blur, "blur");
	ADD_PARAMETER(m_norm_flag, "normalize values");
END_FECORE_CLASS();

IterativeBlur::IterativeBlur(FEModel* fem) : ImageFilter(fem)
{
	m_blur = 0.0;
	m_norm_flag = false;
}

bool IterativeBlur::Init()
{
	return true;
}

//sequential 1D approach
void IterativeBlur::Update(Image& trg, Image& src)
{
	if (m_blur <= 0) { trg = src; return; }
	feLog("Blurring images, blur factor %lg\n", m_blur);

	int K = (int)m_blur;
	float w = float(m_blur) - (float)K;

	int nx = src.width();
	int ny = src.height();
	int nz = src.depth();

	trg = src;
	Image tmp(src);
	
	for (int k = 0; k < K; ++k)
	{
		for (int m = 0; m < 3; ++m)
		{
#ifdef NDEBUG
			#pragma omp parallel for
#endif
			for (int z = 0; z < nz; ++z)
				for (int y = 0; y < ny; ++y)
					for (int x = 0; x < nx; ++x) {
						int m_pos[] = { x, y, z };
						int m_range[] = { nx, ny, nz };
						trg.value(x, y, z) = Apply(tmp, m_pos, m_range, m);
					}
			tmp = trg;
		}
	}

	if (w > 0.0)
	{
		for (int m = 0; m < 3; ++m)
		{
#ifdef NDEBUG
			#pragma omp parallel for
#endif
			for (int z = 0; z < nz; ++z)
				for (int y = 0; y < ny; ++y)
					for (int x = 0; x < nx; ++x) {
						int m_pos[] = { x, y, z };
						int m_range[] = { nx, ny, nz };
						float f1 = Apply(tmp, m_pos, m_range, m);
						float f2 = trg.value(x, y, z);
						trg.value(x, y, z) = Apply(tmp, m_pos, m_range, m);
						trg.value(x, y, z) = f1 * w + f2 * (1.f - w);
					}
			tmp = trg;
		}
	}
}

float IterativeBlur::Apply(Image& img, int m_pos[3], int m_range[3], int m_dir)
{
	int x = m_pos[0]; int y = m_pos[1]; int z = m_pos[2];
	int nx = m_range[0]; int ny = m_range[1]; int nz = m_range[2];
	float f = 0.0;
	if (x > 0) f += img.value(x - 1, y, z); else f += img.value(x, y, z);
	if (x < nx - 1) f += img.value(x + 1, y, z); else f += img.value(x, y, z);
	if (y > 0) f += img.value(x, y - 1, z); else f += img.value(x, y, z);
	if (y < ny - 1) f += img.value(x, y + 1, z); else f += img.value(x, y, z);
	if (z > 0) f += img.value(x, y, z - 1); else f += img.value(x, y, z);
	if (z < nz - 1) f += img.value(x, y, z + 1); else f += img.value(x, y, z);
	return float(0.1666667f) * f;
}

//=================================================================================================
// BoxBlur
//=================================================================================================

BEGIN_FECORE_CLASS(BoxBlur, ImageFilter)
	ADD_PARAMETER(m_blur, "blur"); // blur radius
	ADD_PARAMETER(m_K, "iterations");
	ADD_PARAMETER(m_norm_flag, "normalize values");
	ADD_PARAMETER(m_res, 3, "voxel_resolution");
END_FECORE_CLASS();

BoxBlur::BoxBlur(FEModel* fem) : ImageFilter(fem)
{
	m_blur = 0.0;
	m_K = 0;
	m_norm_flag = false;
	m_res[0] = 1.0; m_res[1] = 1.0; m_res[2] = 1.0;
	m_ri[0] = 1; m_ri[1] = 1; m_ri[2] = 1;
	m_rp = 0.0;
	m_tp = -1.0;
}

bool BoxBlur::Init()
{
	return true;
}

void BoxBlur::Update(Image& trg, Image& src)
{
	if (m_blur < 1) { trg = src; return; }
	// approximate gaussian blur from https://www.ipol.im/pub/art/2013/87/?utm_source=doi
	double r_eff = (int)floor(0.5 * sqrt(((12.0 * m_blur * m_blur / m_K) + 1.0)));
	double sigma_i[3] = { 0.0, 0.0, 0.0 };
	for (int i = 0; i < 3; ++i)
	{
		sigma_i[i] = m_blur / m_res[i]; // (units of px)
		m_ri[i] = (int)floor(0.5 * sqrt(((12.0 * sigma_i[i] * sigma_i[i] / m_K) + 1.0)));
	}
	feLog("Blurring images\n");
	feLog("blur radius = %lg px, rx = %d px, ry = %d px, rz = %d px,\n", r_eff, m_ri[0], m_ri[1], m_ri[2]);
	feLog("blur std = %lg, stdx = %lg, stdy = %lg, stdz = %lg\n", m_blur, sigma_i[0], sigma_i[1], sigma_i[2]);

	int nx = src.width();
	int ny = src.height();
	int nz = src.depth();

	trg = src;
	Image tmp(src);

	// for each iteration
	for (int k = 0; k < m_K; ++k)
	{
		// perform 1D blur along each dimension
		for (int m = 0; m < 3; ++m)
		{
#ifdef NDEBUG
			#pragma omp parallel for
#endif
			for (int z = 0; z < nz; ++z)
				for (int y = 0; y < ny; ++y)
					for (int x = 0; x < nx; ++x)
					{
						int m_pos[] = { x, y, z };
						int m_range[] = { nx, ny, nz };
						
						trg.value(x, y, z) = Apply(tmp, m_pos, m_range, m);
					}
			tmp = trg;
		}
	}
}

float BoxBlur::Apply(Image& img, int m_pos[3], int m_range[3], int m_dir)
{
	//naive implementation
	int x = m_pos[0]; int y = m_pos[1]; int z = m_pos[2];
	int nx = m_range[0]; int ny = m_range[1]; int nz = m_range[2];
	float f = 0.0;
	switch (m_dir)
	{
		case 0:
		{
			for (int i_r = -m_ri[0]; i_r <= m_ri[0]; ++i_r)
			{
				int xr = symmetric_extension(x + i_r, nx);
				f += img.value(xr, y, z);
			}
			f /= (2.f * m_ri[0] + 1.f);
			break;
		}
		case 1:
		{
			for (int i_r = -m_ri[1]; i_r <= m_ri[1]; ++i_r)
			{
				int yr = symmetric_extension(y + i_r, ny);
				f += img.value(x, yr, z);
			}
			f /= (2.f * m_ri[1] + 1.f);
			break;
		}
		case 2:
		{
			for (int i_r = -m_ri[2]; i_r <= m_ri[2]; ++i_r)
			{
				int zr = symmetric_extension(z + i_r, nz);
				f += img.value(x, y, zr);
			}
			f /= (2.f * m_ri[2] + 1.f);
			break;
		}
	}
	return f;
}

//=================================================================================================
// Fastest Fourier Transform in the West (FFTW) Blur
//=================================================================================================

#ifdef HAVE_FFTW
BEGIN_FECORE_CLASS(FFTWBlur, ImageFilter)
ADD_PARAMETER(m_blur, "blur");
ADD_PARAMETER(m_norm_flag, "normalize values");
ADD_PARAMETER(m_res, 3, "voxel_resolution");
END_FECORE_CLASS();

FFTWBlur::FFTWBlur(FEModel* fem) : ImageFilter(fem)
{
	m_blur = 0.0;
	m_norm_flag = false;
	m_res[0] = 1.0; m_res[1] = 1.0; m_res[2] = 1.0;
	m_sigma[0] = 1.0; m_sigma[1] = 1.0; m_sigma[2] = 1.0;
}

bool FFTWBlur::Init() 
{
	return true;
}

void FFTWBlur::Update(Image& trg, Image& src)
{
	if (m_blur < 1) { trg = src; return; }
	// scale blur in each direction for voxel resolution
	for (int i = 0; i < 3; ++i)
		m_sigma[i] = m_blur / m_res[i]; // (units of px)
	feLog("Blurring images\n");
	feLog("blur std = %lg, stdx = %lg, stdy = %lg, stdz = %lg\n", m_blur, m_sigma[0], m_sigma[1], m_sigma[2]);
	fftw_blur_3d(trg, src, m_sigma);
}

// SL: Does nothing for now.
float FFTWBlur::Apply(Image& img, int m_pos[3], int m_range[3], int m_dir)
{
	return 0.f;
}
#endif // HAVE_FFTW