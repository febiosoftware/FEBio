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
#include "FEEFDVerondaWestmann.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEEFDVerondaWestmann, FEUncoupledMaterial)
	ADD_PARAMETER(m_VW.m_c1, "c1")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_VW.m_c2, "c2");
	ADD_PARAMETER(m_EFD.m_beta, 3, "beta");
	ADD_PARAMETER(m_EFD.m_ksi , 3, "ksi" )->setUnits(UNIT_PRESSURE);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEEFDVerondaWestmann::FEEFDVerondaWestmann(FEModel* pfem) : FEUncoupledMaterial(pfem), m_VW(pfem), m_EFD(pfem) 
{
	m_VW.SetParent(this);
	m_EFD.SetParent(this);
}

//-----------------------------------------------------------------------------
bool FEEFDVerondaWestmann::Init()
{
	if (FEUncoupledMaterial::Init() == false) return false;
	if (m_VW.Init() == false) return false;
	if (m_EFD.Init() == false) return false;
	return true;
}

//-----------------------------------------------------------------------------
void FEEFDVerondaWestmann::Serialize(DumpStream& ar)
{
	FEUncoupledMaterial::Serialize(ar);
	m_VW.Serialize(ar);
	m_EFD.Serialize(ar);
}

//-----------------------------------------------------------------------------
mat3ds FEEFDVerondaWestmann::DevStress(FEMaterialPoint& pt)
{
	return m_VW.DevStress(pt) + m_EFD.DevStress(pt);
}

//-----------------------------------------------------------------------------
tens4ds FEEFDVerondaWestmann::DevTangent(FEMaterialPoint& pt)
{
	return m_VW.DevTangent(pt) + m_EFD.DevTangent(pt);
}

//-----------------------------------------------------------------------------
//! calculate deviatoric strain energy density
double FEEFDVerondaWestmann::DevStrainEnergyDensity(FEMaterialPoint& pt)
{
    return m_VW.DevStrainEnergyDensity(pt) + m_EFD.DevStrainEnergyDensity(pt);
}

