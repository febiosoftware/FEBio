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
#include "FEPerfectOsmometer.h"
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEPerfectOsmometer, FEElasticMaterial)
	ADD_PARAMETER(m_phiwr, FE_RANGE_CLOSED(0.0, 1.0), "phiw0");
	ADD_PARAMETER(m_iosm , FE_RANGE_GREATER_OR_EQUAL(0.0), "iosm")->setUnits("n/L^3");
	ADD_PARAMETER(m_bosm , FE_RANGE_GREATER_OR_EQUAL(0.0), "bosm")->setUnits("n/L^3");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEPerfectOsmometer::FEPerfectOsmometer(FEModel* pfem) : FEElasticMaterial(pfem)
{ 
	m_Rgas = 0; 
	m_Tabs = 0; 
}

//-----------------------------------------------------------------------------

bool FEPerfectOsmometer::Init()
{
	if (FEElasticMaterial::Init() == false) return false;
	
	m_Rgas = GetGlobalConstant("R");
	m_Tabs = GetGlobalConstant("T");
	
	if (m_Rgas <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section"); return false; }
	if (m_Tabs <= 0) { feLogError("A positive absolute temperature T must be defined in Globals section");	 return false; }

	return true;
}

//-----------------------------------------------------------------------------
mat3ds FEPerfectOsmometer::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// jacobian
	double J = pt.m_J;
	
	// calculate internal concentration in current configuration
	double iosm = m_iosm*m_phiwr/(J-1+m_phiwr);
	
	// calculate osmotic pressure
	double p = m_Rgas*m_Tabs*(iosm - m_bosm);
	
	// calculate T = -p*I
	mat3dd I(1.0);	// identity tensor
	mat3ds s = -p*I;
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEPerfectOsmometer::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// jacobian
	double J = pt.m_J;

	// calculate internal osmolarity in current configuration
	double iosm = m_iosm*m_phiwr/(J-1+m_phiwr);
	
	// calculate osmotic pressure
	double p = m_Rgas*m_Tabs*(iosm - m_bosm);
	
	// calculate derivative of osmotic pressure w.r.t. J
	double dp = -m_Rgas*m_Tabs*iosm/(J-1+m_phiwr);
	
	mat3dd I(1.0);	// Identity
	
	tens4ds I1 = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	
	// calculate tangent osmotic modulus
	tens4ds c = -J*dp*I1 + p*(2.0*I4 - I1);
	return c;
}

