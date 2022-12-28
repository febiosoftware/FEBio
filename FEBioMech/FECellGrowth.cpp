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
#include "FECellGrowth.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FECellGrowth, FEElasticMaterial)
	ADD_PARAMETER(m_phir, FE_RANGE_GREATER(0.0), "phir");
	ADD_PARAMETER(m_cr  , FE_RANGE_GREATER(0.0), "cr")->setUnits(UNIT_CONCENTRATION);
	ADD_PARAMETER(m_ce  , FE_RANGE_GREATER(0.0), "ce")->setUnits(UNIT_CONCENTRATION);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FECellGrowth::FECellGrowth(FEModel* pfem) : FEElasticMaterial(pfem) 
{ 
	m_Rgas = 0; 
	m_Tabs = 0; 

	m_phir = 0;
	m_cr = 0;
	m_ce = 0;
}

//-----------------------------------------------------------------------------
bool FECellGrowth::Init()
{
	if (FEElasticMaterial::Init() == false) return false;

	m_Rgas = GetGlobalConstant("R");
	m_Tabs = GetGlobalConstant("T");
	
	if (m_Rgas <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section"); return false; }
	if (m_Tabs <= 0) { feLogError("A positive absolute temperature T must be defined in Globals section");	 return false; }

	return true;
}

//-----------------------------------------------------------------------------
mat3ds FECellGrowth::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// jacobian
	double J = pt.m_J;
	
	// calculate intracellular osmolarity relative to mixture volume in reference configuration
	double c = m_cr/(J-m_phir);
	
	// calculate osmotic pressure
	double p = m_Rgas*m_Tabs*(c - m_ce);
	
	// calculate T = -p*I
	mat3dd I(1.0);	// identity tensor
	mat3ds s = -p*I;
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FECellGrowth::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// jacobian
	double J = pt.m_J;
	
	// calculate intracellular osmolarity relative to mixture volume in reference configuration
	double c = m_cr/(J-m_phir);
		
	// calculate osmotic pressure
	double p = m_Rgas*m_Tabs*(c - m_ce);
	
	mat3dd I(1.0);	// Identity
	
	tens4ds I1 = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	
	// calculate tangent osmotic modulus
	tens4ds C = I4*(2*p) - I1*(p-m_Rgas*m_Tabs*c*J/(J-m_phir));
	return C;
}

