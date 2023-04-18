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
#pragma once
#include "feimglib_api.h"

//-----------------------------------------------------------------------------
// This class implements a simple single precision grayscale 3D-image
class FEIMGLIB_API Image
{
public:
	enum ImageFormat {
		RAW8,
		RAW16U
	};

public:
	// constructor
	Image(void);

	// destructor
	~Image(void);

	//! copy constructor
	Image(Image& im);

	//! assignment operator
	Image& operator = (Image& im);

	// allocate storage for image data
	void Create(int nx, int ny, int nz);

	// return size attributes
	int width () { return m_nx; }
	int height() { return m_ny; }
	int depth () { return m_nz; }

	// get a particular data value
	float& value(int x, int y, int z) { return m_pf[(z*m_ny + (m_ny-y-1))*m_nx+x]; }

	// zero image data
	void zero();

	// get the data pointer
	float* data();

protected:
	float*	m_pf;				// image data
	int		m_nx, m_ny, m_nz;	// image dimensions
};

//-----------------------------------------------------------------------------
// helper functions for calculating image derivatives
void image_derive_x(Image& s, Image& d);
void image_derive_y(Image& s, Image& d);
void image_derive_z(Image& s, Image& d);
