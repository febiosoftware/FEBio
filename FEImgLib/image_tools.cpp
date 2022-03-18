#include "image_tools.h"
#include "Image.h"

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
