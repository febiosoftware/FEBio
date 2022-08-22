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
#include "FEImageValuator.h"
#include <FECore/FEMaterialPoint.h>

BEGIN_FECORE_CLASS(FEImageValuator, FEScalarValuator)
	ADD_PARAMETER(m_r0, "range_min");
	ADD_PARAMETER(m_r1, "range_max");
	ADD_PROPERTY(m_imSrc, "image")->SetDefaultType("raw");
	ADD_PROPERTY(m_transform, "transform");
END_FECORE_CLASS();

FEImageValuator::FEImageValuator(FEModel* fem) : FEScalarValuator(fem), m_map(m_im)
{
	m_imSrc = nullptr;
}

bool FEImageValuator::Init()
{
	if (m_imSrc == nullptr) return false;
	if (m_imSrc->GetImage3D(m_im) == false) return false;

	m_map.SetRange(m_r0, m_r1);

	return FEScalarValuator::Init();
}

double FEImageValuator::operator()(const FEMaterialPoint& pt)
{
	return m_transform->value(m_map.value(pt.m_r0));
}

FEScalarValuator* FEImageValuator::copy()
{
	assert(false);
	return nullptr;
}
