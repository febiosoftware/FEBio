#include "stdafx.h"
#include "Image.h"
#include <stdio.h>
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
bool Image::Load(const char *szfile)
{
	FILE* fp = fopen(szfile, "rb");
	if (fp == 0) return false;

	int n = m_nx*m_ny*m_nz;
	unsigned char* pb = new unsigned char[n];
	fread(pb, n, 1, fp);
	for (int k=0; k<m_nz; ++k) 
		for (int j=0; j<m_ny; ++j) 
			for (int i=0; i<m_nx; ++i) 
			{
				float f =  (float) pb[(k*m_ny + j)*m_nx + i] / 255.f;
				m_pf[(k*m_ny + (m_ny-j-1))*m_nx + i] = f;
			}

	fclose(fp);
	return true;
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
