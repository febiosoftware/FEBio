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
#include "Image.h"
#include <memory.h>

//-----------------------------------------------------------------------------
Image::Image(void)
{
	m_pf = 0;
	m_nx = m_ny = m_nz = 0;
}

//-----------------------------------------------------------------------------
Image::~Image(void)
{
	delete [] m_pf;
	m_pf = 0;
}

//-----------------------------------------------------------------------------
void Image::Create(int nx, int ny, int nz)
{
	if (m_pf) delete [] m_pf;
	m_nx = nx;
	m_ny = ny;
	m_nz = nz;
	m_pf = new float[m_nx*m_ny*m_nz];
}

//-----------------------------------------------------------------------------
Image::Image(Image& im)
{
	m_pf = 0;
	Create(im.width(), im.height(), im.depth());
	memcpy(m_pf, im.m_pf, m_nx*m_ny*m_nz*sizeof(float));
}

//-----------------------------------------------------------------------------
Image& Image::operator = (Image& im)
{
	if ((m_nx != im.m_nx)||(m_ny != im.m_ny)||(m_nz != im.m_nz)) Create(im.width(), im.height(), im.depth());
	memcpy(m_pf, im.m_pf, m_nx*m_ny*m_nz*sizeof(float));
	return (*this);
}

//-----------------------------------------------------------------------------
void Image::zero()
{
	int n = m_nx*m_ny*m_nz;
	for (int i=0; i<n; ++i) m_pf[i] = 0.f;
}

//-----------------------------------------------------------------------------
float* Image::data()
{
	return m_pf;
}

//-----------------------------------------------------------------------------
void image_derive_x(Image& s, Image& d)
{
	int nx = s.width();
	int ny = s.height();
	int nz = s.depth();
	for (int k=0; k<nz; ++k)
	{
		for (int j=0; j<ny; ++j)
		{
			d.value(0, j, k) = s.value(1, j, k) - s.value(0, j, k);
			for (int i=1; i<nx-1; ++i) d.value(i, j, k) = (s.value(i+1, j, k) - s.value(i-1, j, k))*0.5f;
			d.value(nx-1, j, k) = s.value(nx-1, j, k) - s.value(nx-2, j, k);
		}
	}
}

//-----------------------------------------------------------------------------
void image_derive_y(Image& s, Image& d)
{
	int nx = s.width();
	int ny = s.height();
	int nz = s.depth();
	for (int k=0; k<nz; ++k)
	{
		for (int i=0; i<nx; ++i)
		{
			d.value(i, 0, k) = s.value(i, 1, k) - s.value(i, 0, k);
			for (int j=1; j<ny-1; ++j) d.value(i, j, k) = (s.value(i, j+1, k) - s.value(i, j-1, k)) *0.5f;
			d.value(i, ny-1, k) = s.value(i, ny-1, k) - s.value(i, ny-2, k);
		}
	}
}

//-----------------------------------------------------------------------------
void image_derive_z(Image& s, Image& d)
{
	int nx = s.width();
	int ny = s.height();

	int nz = s.depth();
	if (nz == 1) { d.zero(); return; }

	for (int j=0; j<ny; ++j)
	{
		for (int i=0; i<nx; ++i)
		{
			d.value(i, j, 0) = s.value(i, j, 1) - s.value(i, j, 0);
			for (int k=1; k<nz-1; ++k) d.value(i, j, k) = (s.value(i, j, k+1) - s.value(i, j, k-1)) *0.5f;
			d.value(i, j, nz-1) = s.value(i, j, nz-1) - s.value(i, j, nz-2);
		}
	}
}
