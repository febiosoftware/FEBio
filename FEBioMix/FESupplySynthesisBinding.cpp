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
#include "FESupplySynthesisBinding.h"

// define the material parameters
BEGIN_FECORE_CLASS(FESupplySynthesisBinding, FESoluteSupply)
	ADD_PARAMETER(m_supp, "supp");
	ADD_PARAMETER(m_kf  , FE_RANGE_GREATER         (0.0), "kf"  );
	ADD_PARAMETER(m_kr  , FE_RANGE_GREATER_OR_EQUAL(0.0), "kr"  );
	ADD_PARAMETER(m_crt , FE_RANGE_GREATER_OR_EQUAL(0.0), "Rtot");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FESupplySynthesisBinding::FESupplySynthesisBinding(FEModel* pfem) : FESoluteSupply(pfem)
{
	m_supp = m_kf = m_kr = m_crt = 0;
}

//-----------------------------------------------------------------------------
//! Solute supply
double FESupplySynthesisBinding::Supply(FEMaterialPoint& mp)
{
	double crhat = m_supp-ReceptorLigandSupply(mp);
	return crhat;
}

//-----------------------------------------------------------------------------
//! Tangent of solute supply with respect to strain
double FESupplySynthesisBinding::Tangent_Supply_Strain(FEMaterialPoint &mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Tangent of solute supply with respect to referential concentration
double FESupplySynthesisBinding::Tangent_Supply_Concentration(FEMaterialPoint &mp)
{
	FESolutesMaterialPoint& pt = *mp.ExtractData<FESolutesMaterialPoint>();
	
	double crc = pt.m_sbmr[0];
	double dcrhatdcr = -m_kf*(m_crt - crc);
	
	return dcrhatdcr;
}

//-----------------------------------------------------------------------------
//! Receptor-ligand complex supply
double FESupplySynthesisBinding::ReceptorLigandSupply(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
	double J = et.m_J;
	double ca = spt.m_ca[0];
	double phi0 = ppt.m_phi0t;
	double cr = (J-phi0)*ca;
	double crc = spt.m_sbmr[0];
	double crchat = m_kf*cr*(m_crt-crc) - m_kr*crc;
	return crchat;
}

//-----------------------------------------------------------------------------
//! Solute supply at steady state
double FESupplySynthesisBinding::SupplySS(FEMaterialPoint& mp)
{
	return m_supp;
}

//-----------------------------------------------------------------------------
//! Receptor-ligand concentration at steady-state
double FESupplySynthesisBinding::ReceptorLigandConcentrationSS(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
	double J = et.m_J;
	double ca = spt.m_ca[0];
	double phi0 = ppt.m_phi0t;
	double cr = (J-phi0)*ca;
	double Kd = m_kr/m_kf;	// dissociation constant
	double crc = m_crt*cr/(Kd+cr);
	return crc;
}

//-----------------------------------------------------------------------------
//! Referential solid supply (moles of solid/referential volume/time)
double FESupplySynthesisBinding::SolidSupply(FEMaterialPoint& mp)
{
	return ReceptorLigandSupply(mp);
}

//-----------------------------------------------------------------------------
//! Referential solid concentration (moles of solid/referential volume)
//! at steady-state
double FESupplySynthesisBinding::SolidConcentrationSS(FEMaterialPoint& mp)
{
	return 0;
}
