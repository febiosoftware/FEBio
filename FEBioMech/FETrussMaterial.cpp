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



#include "stdafx.h"
#include "FETrussMaterial.h"

// define the material parameters
BEGIN_FECORE_CLASS(FETrussMaterial, FEMaterial)
	ADD_PARAMETER(m_rho, FE_RANGE_GREATER(0.0), "density");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FETrussMaterial::FETrussMaterial(FEModel* pfem) : FEMaterial(pfem) 
{
	m_rho = 1.0;
}

//-----------------------------------------------------------------------------
FETrussMaterial::~FETrussMaterial() 
{
}

//-----------------------------------------------------------------------------
//! material density
double FETrussMaterial::Density()
{
	return m_rho;
}

//=============================================================================
// define the material parameters
BEGIN_FECORE_CLASS(FELinearTrussMaterial, FETrussMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FELinearTrussMaterial::FELinearTrussMaterial(FEModel* fem) : FETrussMaterial(fem)
{
	m_E = 0.0;
}

//-----------------------------------------------------------------------------
// Note that this function returns the Kirchhoff stress!
double FELinearTrussMaterial::Stress(FEMaterialPoint &mp)
{
	FETrussMaterialPoint& pt = *mp.ExtractData<FETrussMaterialPoint>();
	return m_E*log(pt.m_l);
}

//-----------------------------------------------------------------------------
double FELinearTrussMaterial::Tangent(FEMaterialPoint &pt)
{
	return m_E;
}
