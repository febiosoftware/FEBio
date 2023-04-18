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
#include "FEEFDNeoHookean.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEEFDNeoHookean, FEElasticMaterial)
	ADD_PARAMETER(m_NH.m_E, "E")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_NH.m_v, "v");
	ADD_PARAMETER(m_EFD.m_beta, 3, "beta");
	ADD_PARAMETER(m_EFD.m_ksi , 3, "ksi" )->setUnits(UNIT_PRESSURE);
END_FECORE_CLASS();

//////////////////////////////////////////////////////////////////////
// FEEFDNeoHookean
//////////////////////////////////////////////////////////////////////

FEEFDNeoHookean::FEEFDNeoHookean(FEModel* pfem) : FEElasticMaterial(pfem), m_EFD(pfem), m_NH(pfem) 
{
	int a = 0;
}

bool FEEFDNeoHookean::Init()
{
	if (FEElasticMaterial::Init() == false) return false;
	if (m_NH.Init()  == false) return false;
	if (m_EFD.Init() == false) return false;
	return true;
}

void FEEFDNeoHookean::Serialize(DumpStream& ar)
{
	FEElasticMaterial::Serialize(ar);
	m_NH.Serialize(ar);
	m_EFD.Serialize(ar);
}

mat3ds FEEFDNeoHookean::Stress(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	mat3ds s = m_NH.Stress(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// evaluate stress and add it to matrix contribution
	s += m_EFD.Stress(mp);
	
	return s;
}

tens4ds FEEFDNeoHookean::Tangent(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	tens4ds c = m_NH.Tangent(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// evaluate stress and add it to matrix contribution
	c += m_EFD.Tangent(mp);
	
	return c;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FEEFDNeoHookean::StrainEnergyDensity(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	double sed = m_NH.StrainEnergyDensity(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	sed += m_EFD.StrainEnergyDensity(mp);
    
    return sed;
}

//////////////////////////////////////////////////////////////////////
// FEEFDNeoHookeanOld
//////////////////////////////////////////////////////////////////////

// define the material parameters
BEGIN_FECORE_CLASS(FEEFDNeoHookeanOld, FEElasticMaterial)
	ADD_PARAMETER(m_NH.m_E, "E");
	ADD_PARAMETER(m_NH.m_v, "v");
	ADD_PARAMETER(m_EFD.m_beta, 3, "beta");
	ADD_PARAMETER(m_EFD.m_ksi , 3, "ksi" );
END_FECORE_CLASS();

bool FEEFDNeoHookeanOld::Init()
{
	if (FEElasticMaterial::Init() == false) return false;

	if (m_NH.Init() == false) return false;
	if (m_EFD.Init() == false) return false;
	return true;
}

void FEEFDNeoHookeanOld::Serialize(DumpStream& ar)
{
	FEElasticMaterial::Serialize(ar);
	m_NH.Serialize(ar);
	m_EFD.Serialize(ar);
}

mat3ds FEEFDNeoHookeanOld::Stress(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	mat3ds s = m_NH.Stress(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// evaluate stress and add it to matrix contribution
	s += m_EFD.Stress(mp);
	
	return s;
}

tens4ds FEEFDNeoHookeanOld::Tangent(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	tens4ds c = m_NH.Tangent(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// evaluate stress and add it to matrix contribution
	c += m_EFD.Tangent(mp);
	
	return c;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FEEFDNeoHookeanOld::StrainEnergyDensity(FEMaterialPoint& mp)
{
	// --- M A T R I X   C O N T R I B U T I O N ---
	double sed = m_NH.StrainEnergyDensity(mp);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	sed += m_EFD.StrainEnergyDensity(mp);
    
    return sed;
}
