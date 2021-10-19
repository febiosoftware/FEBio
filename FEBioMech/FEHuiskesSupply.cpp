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
#include "FEHuiskesSupply.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEHuiskesSupply, FESolidSupply)
	ADD_PARAMETER(m_B, "B");
	ADD_PARAMETER(m_k, "k");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEHuiskesSupply::FEHuiskesSupply(FEModel* pfem) : FESolidSupply(pfem)
{
	m_B = m_k = 0;
}

//-----------------------------------------------------------------------------
//! Solid supply
double FEHuiskesSupply::Supply(FEMaterialPoint& mp)
{
	FERemodelingMaterialPoint& rpt = *mp.ExtractData<FERemodelingMaterialPoint>();
	double rhor = rpt.m_rhor;
	double sed = rpt.m_sed;
	double rhorhat = m_B*(sed/rhor - m_k);
	return rhorhat;
}

//-----------------------------------------------------------------------------
//! Tangent of solid supply with respect to strain
mat3ds FEHuiskesSupply::Tangent_Supply_Strain(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FERemodelingMaterialPoint& rpt = *mp.ExtractData<FERemodelingMaterialPoint>();
    mat3ds ruhat = pt.m_s*(m_B/rpt.m_rhor);
	return ruhat;
}

//-----------------------------------------------------------------------------
//! Tangent of solid supply with respect to referential density
double FEHuiskesSupply::Tangent_Supply_Density(FEMaterialPoint &mp)
{
	FERemodelingMaterialPoint& rpt = *mp.ExtractData<FERemodelingMaterialPoint>();
    double rhor = rpt.m_rhor;
    double sed = rpt.m_sed;
    double dsed = rpt.m_dsed;
	return (dsed - sed/rhor)*m_B/rhor;
}

