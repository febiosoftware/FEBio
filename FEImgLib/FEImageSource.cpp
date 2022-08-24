/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/
#include "FEImageSource.h"

FEImageSource::FEImageSource(FEModel* fem) : FECoreClass(fem)
{

}

//========================================================================

BEGIN_FECORE_CLASS(FERawImage, FEImageSource)
	ADD_PARAMETER(m_file, "file", FE_PARAM_ATTRIBUTE);
	ADD_PARAMETER(m_format, "format", 0, "RAW8\0RAW16U\0");
	ADD_PARAMETER(m_dim, 3, "size");
	ADD_PARAMETER(m_bend, "endianess");
END_FECORE_CLASS();

FERawImage::FERawImage(FEModel* fem) : FEImageSource(fem)
{
	m_dim[0] = m_dim[1] = m_dim[2] = 0;
	m_bend = false;
	m_format = 0;
}

bool FERawImage::GetImage3D(Image& im)
{
	// get the file name
	const char* szfile = m_file.c_str();

	int* n = m_dim;
	im.Create(n[0], n[1], n[2]);

	// Try to load the image file
	Image::ImageFormat fmt;
	switch (m_format)
	{
	case 0: fmt = Image::RAW8; break;
	case 1: fmt = Image::RAW16U; break;
	default:
		assert(false);
		return false;
	}
	return Load(szfile, im, fmt, m_bend);
}

//-----------------------------------------------------------------------------
bool FERawImage::Load(const char* szfile, Image& im, Image::ImageFormat fmt, bool endianess)
{
	// check the format
	int m = -1;
	switch (fmt)
	{
	case Image::RAW8: m = 1; break;
	case Image::RAW16U: m = 2; break;
	default:
		return false;
		break;
	}

	// open the file
	FILE* fp = fopen(szfile, "rb");
	if (fp == 0) return false;

	// read the data
	int nx = im.width();
	int ny = im.height();
	int nz = im.depth();
	int voxels = nx * ny * nz;
	int nsize = voxels * m;
	unsigned char* pb = new unsigned char[nsize];
	fread(pb, nsize, 1, fp);

	float* pf = im.data();

	if (fmt == Image::RAW8)
	{
		for (int k = 0; k < nz; ++k)
			for (int j = 0; j < ny; ++j)
				for (int i = 0; i < nx; ++i)
				{
					float f = (float)pb[(k * ny + j) * nx + i] / 255.f;
					pf[(k * ny + (ny - j - 1)) * nx + i] = f;
				}
	}
	else if (fmt == Image::RAW16U)
	{
		for (int k = 0; k < nz; ++k)
			for (int j = 0; j < ny; ++j)
				for (int i = 0; i < nx; ++i)
				{
					unsigned char* d = pb + 2 * ((k * ny + j) * nx + i);

					unsigned short n = *((unsigned short*)d);

					if (endianess)
					{
						unsigned int n1 = (unsigned int)d[0];
						unsigned int n2 = (unsigned int)d[1];
						n = (n1 << 8) + n2;
					}

					float f = (float)n / 65535.f;
					pf[(k * ny + (ny - j - 1)) * nx + i] = f;
				}
	}

	// all done
	delete[] pb;
	fclose(fp);
	return true;
}
