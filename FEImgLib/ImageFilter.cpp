#include "Image.h"
#include "ImageFilter.h"
#include <vector>
ImageFilter::ImageFilter(void) {}

//=================================================================================================
// IterativeBlur
//=================================================================================================

IterativeBlur::IterativeBlur() {}
//sequential 1D approach
void IterativeBlur::eval1D(Image& trg, Image& src, float d)
{
	if (d <= 0) { trg = src; return; }

	int n = (int)d;
	float w = d - (float)n;

	int nx = src.width();
	int ny = src.height();
	int nz = src.depth();

	trg = src;
	Image tmp(src);
	
	for (int l = 0; l < n; ++l)
	{
		for (int m = 0; m < 2; ++m)
		{
			for (int k = 0; k < nz; ++k)
				for (int j = 0; j < ny; ++j)
					for (int i = 0; i < nx; ++i) {
						int m_pos[] = { i, j, k };
						int m_range[] = { nx, ny, nz };
						trg.value(i, j, k) = apply1D(tmp, m_pos, m_range, m);
					}

			tmp = trg;
		}
	}

	if (w > 0.0)
	{
		for (int m = 0; m < 2; ++m)
		{
			for (int k = 0; k < nz; ++k)
				for (int j = 0; j < ny; ++j)
					for (int i = 0; i < nx; ++i) {
						int m_pos[] = { i, j, k };
						int m_range[] = { nx, ny, nz };
						trg.value(i, j, k) = apply1D(tmp, m_pos, m_range, m);
					}

			tmp = trg;
		}
	}
}

//naive approach
void IterativeBlur::eval3D(Image& trg, Image& src, float d)
{
	if (d <= 0) { trg = src; return; }

	int n = (int)d;
	float w = d - (float)n;

	int nx = src.width();
	int ny = src.height();
	int nz = src.depth();

	trg = src;
	Image tmp(src);
	
	//IterativeBlur* m_filt = new IterativeBlur();
	for (int l = 0; l < n; ++l)
	{
		#pragma omp parallel for collapse(3)
		for (int k = 0; k < nz; ++k)
			for (int j = 0; j < ny; ++j)
				for (int i = 0; i < nx; ++i) {
					int m_pos[] = { i, j, k };
					int m_range[] = { nx, ny, nz };
					trg.value(i, j, k) = apply3D(tmp, m_pos, m_range);
				}
		tmp = trg;
	}

	if (w > 0.0)
	{
		#pragma omp parallel for collapse(3)
		for (int k = 0; k < nz; ++k)
			for (int j = 0; j < ny; ++j)
				for (int i = 0; i < nx; ++i)
				{
					int m_pos[] = { i, j, k };
					int m_range[] = { nx, ny, nz };
					float f1 = apply3D(tmp, m_pos, m_range);
					float f2 = trg.value(i, j, k);
					trg.value(i, j, k) = f1 * w + f2 * (1.f - w);
				}
	}
}

double IterativeBlur::apply1D(Image& img, int m_pos[3], int m_range[3], int m_dir)
{
	// sequential 1D implementation
	int i = m_pos[0]; int j = m_pos[1]; int k = m_pos[2];
	int nx = m_range[0]; int ny = m_range[1]; int nz = m_range[2];
	float f[2] = { 0.0f, 0.0f };
	switch (m_dir)
	{
	case 0:
		if (i > 0) f[0] = img.value(i - 1, j, k); else f[0] = img.value(i, j, k);
		if (i < nx - 1) f[1] = img.value(i + 1, j, k); else f[1] = img.value(i, j, k);
		break;
	case 1:
		if (j > 0) f[0] = img.value(i, j - 1, k); else f[0] = img.value(i, j, k);
		if (j < ny - 1) f[1] = img.value(i, j + 1, k); else f[1] = img.value(i, j, k);
		break;
	case 2:
		if (k > 0) f[0] = img.value(i, j, k - 1); else f[0] = img.value(i, j, k);
		if (k < nz - 1) f[1] = img.value(i, j, k + 1); else f[1] = img.value(i, j, k);
		break;
	}
	return 0.5f * (f[0] + f[1]);
}

double IterativeBlur::apply3D(Image& img, int m_pos[3], int m_range[3])
{
	//naive implementation
	int i = m_pos[0]; int j = m_pos[1]; int k = m_pos[2];
	int nx = m_range[0]; int ny = m_range[1]; int nz = m_range[2];
	float f[6];
	if (i > 0) f[0] = img.value(i - 1, j, k); else f[0] = img.value(i, j, k);
	if (i < nx - 1) f[1] = img.value(i + 1, j, k); else f[1] = img.value(i, j, k);
	if (j > 0) f[2] = img.value(i, j - 1, k); else f[2] = img.value(i, j, k);
	if (j < ny - 1) f[3] = img.value(i, j + 1, k); else f[3] = img.value(i, j, k);
	if (k > 0) f[4] = img.value(i, j, k - 1); else f[4] = img.value(i, j, k);
	if (k < nz - 1) f[5] = img.value(i, j, k + 1); else f[5] = img.value(i, j, k);
	return 0.1666667f * (f[0] + f[1] + f[2] + f[3] + f[4] + f[5]);
}

//=================================================================================================
// BoxBlur
//=================================================================================================

BoxBlur::BoxBlur() {}

void BoxBlur::eval1D(Image& trg, Image& src, float d)
{
	if (d <= 0) { trg = src; return; }

	int n = (int)d;
	float w = d - (float)n;

	int nx = src.width();
	int ny = src.height();
	int nz = src.depth();

	trg = src;
	Image tmp(src);

	//IterativeBlur* m_filt = new IterativeBlur();
	for (int l = 0; l < n; ++l)
	{
#pragma omp parallel for collapse(3)
		for (int k = 0; k < nz; ++k)
			for (int j = 0; j < ny; ++j)
				for (int i = 0; i < nx; ++i) {
					int m_pos[] = { i, j, k };
					int m_range[] = { nx, ny, nz };
					trg.value(i, j, k) = apply3D(tmp, m_pos, m_range);
				}
		tmp = trg;
	}

	if (w > 0.0)
	{
#pragma omp parallel for collapse(3)
		for (int k = 0; k < nz; ++k)
			for (int j = 0; j < ny; ++j)
				for (int i = 0; i < nx; ++i)
				{
					int m_pos[] = { i, j, k };
					int m_range[] = { nx, ny, nz };
					float f1 = apply3D(tmp, m_pos, m_range);
					float f2 = trg.value(i, j, k);
					trg.value(i, j, k) = f1 * w + f2 * (1.f - w);
				}
	}
}

void BoxBlur::eval3D(Image& trg, Image& src, float d)
{
	if (d <= 0) { trg = src; return; }

	int n = (int)d;
	float w = d - (float)n;

	int nx = src.width();
	int ny = src.height();
	int nz = src.depth();

	trg = src;
	Image tmp(src);

	//IterativeBlur* m_filt = new IterativeBlur();
	for (int l = 0; l < n; ++l)
	{
#pragma omp parallel for collapse(3)
		for (int k = 0; k < nz; ++k)
			for (int j = 0; j < ny; ++j)
				for (int i = 0; i < nx; ++i) {
					int m_pos[] = { i, j, k };
					int m_range[] = { nx, ny, nz };
					trg.value(i, j, k) = apply3D(tmp, m_pos, m_range);
				}
		tmp = trg;
	}

	if (w > 0.0)
	{
#pragma omp parallel for collapse(3)
		for (int k = 0; k < nz; ++k)
			for (int j = 0; j < ny; ++j)
				for (int i = 0; i < nx; ++i)
				{
					int m_pos[] = { i, j, k };
					int m_range[] = { nx, ny, nz };
					float f1 = apply3D(tmp, m_pos, m_range);
					float f2 = trg.value(i, j, k);
					trg.value(i, j, k) = f1 * w + f2 * (1.f - w);
				}
	}
}

double BoxBlur::apply1D(Image& img, int m_pos[3], int m_range[3], int m_dir)
{
	//naive implementation
	int i = m_pos[0]; int j = m_pos[1]; int k = m_pos[2];
	int nx = m_range[0]; int ny = m_range[1]; int nz = m_range[2];
	float f = 0.0;
	
	//set flags based on position
	bool flag_i0, flag_in, flag_j0, flag_jn, flag_k0, flag_kn = false;
	int n_flag = 0;
	int n_count = 0;
	int i_0, j_0, k_0 = -1;
	int i_f, j_f, k_f = 1;
	if (i == 0) { flag_i0 = true; n_flag++; i_0 = 0; }
	if (i == nx) { flag_in = true; n_flag++; i_f = 0; }
	if (j == 0) { flag_j0 = true; n_flag++; j_0 = 0; }
	if (j == ny) { flag_jn = true; n_flag++; j_f = 0; }
	if (k == 0) { flag_k0 = true; n_flag++; k_0 = 0; }
	if (k == nz) { flag_kn = true; n_flag++; k_f = 0; }

	switch (n_flag) {
	case 0:
		n_count = 27;
		break;
	case 1:
		n_count = 18;
		break;
	case 2:
		n_count = 12;
		break;
	case 3:
		n_count = 8;
		break;	
	}

	for (int a = i_0; a < i_f; ++a)
		for (int b = j_0; b < j_f; ++b)
			for (int c = k_0; c < k_f; ++c)
			{
				f += img.value(i + a, j + b, c + k);
			}	
	
	return (f / n_count);
}

double BoxBlur::apply3D(Image& img, int m_pos[3], int m_range[3])
{
	// naive implementation
	int i = m_pos[0]; int j = m_pos[1]; int k = m_pos[2];
	int nx = m_range[0]; int ny = m_range[1]; int nz = m_range[2];
	float f = 0.0;

	// set flags based on position
	bool flag_i0, flag_in, flag_j0, flag_jn, flag_k0, flag_kn = false;
	int n_flag = 0;	int n_count = 0; 
	int i_0 = -1; int j_0 = -1; int k_0 = -1;
	int i_f = 1; int j_f = 1; int k_f = 1;
	if (i == 0) { flag_i0 = true; n_flag++; i_0 = 0; }
	if (i == nx - 1) { flag_in = true; n_flag++; i_f = 0; }
	if (j == 0) { flag_j0 = true; n_flag++; j_0 = 0; }
	if (j == ny - 1) { flag_jn = true; n_flag++; j_f = 0; }
	if (k == 0) { flag_k0 = true; n_flag++; k_0 = 0; }
	if (k == nz - 1) { flag_kn = true; n_flag++; k_f = 0; }

	// determine # of neighbors
	switch (n_flag) {
	case 0:
		n_count = 27;
		break;
	case 1:
		n_count = 18;
		break;
	case 2:
		n_count = 12;
		break;
	case 3:
		n_count = 8;
		break;
	}

	// get unweighted average value from neighbors
	for (int a = i_0; a <= i_f; ++a)
		for (int b = j_0; b <= j_f; ++b)
			for (int c = k_0; c <= k_f; ++c)
				f += img.value(i + a, j + b, c + k);

	return (f / n_count);
}