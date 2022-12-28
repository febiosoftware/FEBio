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
#include "Image.h"
#include <FECore/FECoreClass.h>
#include "feimglib_api.h"

//---------------------------------------------------------------------------
// Base class for image sources. 
class FEIMGLIB_API FEImageSource : public FECoreClass
{
	FECORE_BASE_CLASS(FEImageSource)

public:
	FEImageSource(FEModel* fem);

	virtual bool GetImage3D(Image& im) = 0;
};

//---------------------------------------------------------------------------
// Class for reading raw images
class FEIMGLIB_API FERawImage : public FEImageSource
{
public:
	FERawImage(FEModel* fem);

	bool GetImage3D(Image& im) override;

private:
	// load raw data from file
	bool Load(const char* szfile, Image& im, Image::ImageFormat fmt, bool endianess = false);

protected:
	std::string		m_file;
	int				m_dim[3];
	int				m_format;
	bool			m_bend;

	DECLARE_FECORE_CLASS();
};
