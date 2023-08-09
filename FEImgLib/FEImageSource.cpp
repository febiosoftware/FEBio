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
#include <FECore/log.h>

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

//=============================================================================
BEGIN_FECORE_CLASS(FENRRDImage, FEImageSource)
	ADD_PARAMETER(m_file, "file");
END_FECORE_CLASS();

FENRRDImage::FENRRDImage(FEModel* fem) : FEImageSource(fem)
{
}

bool FENRRDImage::GetImage3D(Image& im)
{
	if (m_file.empty()) return false;
	const char* szfile = m_file.c_str();
	return Load(szfile, im);
}

bool nrrd_read_line(FILE* fp, char* buf)
{
	char c;
	while (fread(&c, 1, 1, fp))
	{
		if ((c != '\r') && (c != '\n')) *buf++ = c;
		if (c == '\n') {
			*buf = 0; return true;
		}
	}
	return false;
}

bool nrrd_read_key_value(const char* szline, char* key, char* val)
{
	const char* c = szline;
	char* d = key;
	bool skip_space = true;
	while (*c)
	{
		if (*c == ':')
		{
			*d = 0;
			d = val;
			skip_space = true;
		}
		else if (skip_space) {
			if (isspace(*c) == 0) {
				*d++ = *c;
				skip_space = false;
			}
		}
		else *d++ = *c;
		c++;
	}
	*d = 0;
	return true;
}

bool FENRRDImage::Load(const char* szfile, Image& im)
{
	if (szfile == nullptr) return false;

	FILE* fp = fopen(szfile, "rb");
	if (fp == nullptr) return false;

	// read the header
	char buf[1024] = { 0 }, key[256] = { 0 }, val[256] = { 0 };
	nrrd_read_line(fp, buf);
	if (strncmp(buf, "NRRD", 4) != 0) { fclose(fp); return false; }

	NRRD_TYPE imType = NRRD_INVALID;
	int dim = 0;
	int sizes[3] = { 0 };
	NRRD_ENCODING enc = NRRD_RAW;

	while (nrrd_read_line(fp, buf))
	{
		if (buf[0] == 0) break;

		// read key-value
		if (buf[0] != '#')
		{
			if (nrrd_read_key_value(buf, key, val) == false) { fclose(fp); return false; }

			if (strcmp(key, "type") == 0)
			{
				if (strcmp(val, "float") == 0) imType = NRRD_FLOAT;
			}
			else if (strcmp(key, "dimension") == 0) dim = atoi(val);
			else if (strcmp(key, "sizes") == 0) sscanf(val, "%d %d %d", sizes, sizes + 1, sizes + 2);
			else if (strcmp(key, "encoding") == 0)
			{
				if (strcmp(val, "raw") == 0) enc = NRRD_RAW;
				if (strcmp(val, "ascii") == 0) enc = NRRD_ASCII;
			}
		}
	}

	// read the image
	if (imType != NRRD_FLOAT) { feLog("Invalid image format (only float supported)"); fclose(fp); return false; }
	if (dim != 3) { feLog("Invalid image dimensions (must be 3)"); fclose(fp); return false; }
	if (enc != NRRD_RAW) { feLog("Invalid encoding (only raw supported)"); fclose(fp); return false; }

	int nsize = sizes[0] * sizes[1] * sizes[2];
	if (nsize <= 0) { feLog("Invalid image sizes"); fclose(fp); return false; }

	im.Create(sizes[0], sizes[1], sizes[2]);
	float* pf = im.data();
	int nread = fread(pf, sizeof(float), nsize, fp);
	if (nread != nsize) { feLog("Failed reading image data"); fclose(fp); return false; }
	
	fclose(fp);

	return true;
}
