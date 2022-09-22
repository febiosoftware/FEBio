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
#include "FEBiphasicSolute.h"
#include <FECore/FEModel.h>
#include <FECore/FECoreKernel.h>
#include <FECore/log.h>

//=============================================================================
//                 B I P H A S I C S O L U T E
//=============================================================================

//-----------------------------------------------------------------------------
// Material parameters for the FEBiphasicSolute material
BEGIN_FECORE_CLASS(FEBiphasicSolute, FEMaterial)
	ADD_PARAMETER(m_phi0 , FE_RANGE_CLOSED(0.0, 1.0)     , "phi0");
	ADD_PARAMETER(m_rhoTw, FE_RANGE_GREATER_OR_EQUAL(0.0), "fluid_density");

	// set material properties
	ADD_PROPERTY(m_pSolid , "solid", FEProperty::Required | FEProperty::TopLevel);
	ADD_PROPERTY(m_pPerm  , "permeability");
	ADD_PROPERTY(m_pOsmC  , "osmotic_coefficient");
	ADD_PROPERTY(m_pSolute, "solute");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! FEBiphasicSolute constructor

FEBiphasicSolute::FEBiphasicSolute(FEModel* pfem) : FEMaterial(pfem)
{
	m_phi0 = 0;
	m_rhoTw = 0;
	m_rhoTu = 0;
	m_Mu = 0;
	m_Rgas = 0;
	m_Tabs = 0; 

	m_pSolid = 0;
	m_pPerm = 0;
	m_pOsmC = 0;
	m_pSolute = 0;
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEBiphasicSolute::CreateMaterialPointData() 
{
	FEBiphasicMaterialPoint* pbp = new FEBiphasicMaterialPoint(m_pSolid->CreateMaterialPointData());
	return new FESolutesMaterialPoint(pbp);
}

//-----------------------------------------------------------------------------
bool FEBiphasicSolute::Init()
{
	// we need to set the solute ID before we call FEMaterial::Init()
	// because it is used in FESolute::Init()
	m_pSolute->SetSoluteLocalID(0);

    if (!m_pSolid->Init()) return false;
    if (!m_pPerm->Init()) return false;
    if (!m_pOsmC->Init()) return false;
    if (!m_pSolute->Init()) return false;
    
	// Call base class which calls the Init member of all properties
	if (FEMaterial::Init() == false) return false;
	
	m_Rgas = GetFEModel()->GetGlobalConstant("R");
	m_Tabs = GetFEModel()->GetGlobalConstant("T");
	
	if (m_Rgas <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section"); return false; }
	if (m_Tabs <= 0) { feLogError("A positive absolute temperature T must be defined in Globals section");	 return false; }

	return true;
}

//-----------------------------------------------------------------------------
// update specialized material points
void FEBiphasicSolute::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp)
{
    m_pSolid->UpdateSpecializedMaterialPoints(mp, tp);
    m_pPerm->UpdateSpecializedMaterialPoints(mp, tp);
    m_pOsmC->UpdateSpecializedMaterialPoints(mp, tp);
    m_pSolute->UpdateSpecializedMaterialPoints(mp, tp);
}

//-----------------------------------------------------------------------------
void FEBiphasicSolute::Serialize(DumpStream& ar)
{
	FEMaterial::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_Rgas & m_Tabs & m_Mu;
}

//-----------------------------------------------------------------------------
//! Porosity in current configuration
double FEBiphasicSolute::Porosity(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& pet = *pt.ExtractData<FEBiphasicMaterialPoint>();
	
	// relative volume
	double J = et.m_J;
	// porosity
//	double phiw = 1 - m_phi0/J;
	double phi0 = pet.m_phi0t;
	double phiw = 1 - phi0/J;
	// check for pore collapse
	// TODO: throw an error if pores collapse
	phiw = (phiw > 0) ? phiw : 0;
	
	return phiw;
}

//-----------------------------------------------------------------------------
//! The stress of a solute-poroelastic material is the sum of the fluid pressure
//! and the elastic stress. Note that this function is declared in the base class
//! so you do not have to reimplement it in a derived class, unless additional
//! pressure terms are required.

mat3ds FEBiphasicSolute::Stress(FEMaterialPoint& mp)
{
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
	// calculate solid material stress
	mat3ds s = m_pSolid->Stress(mp);
	
	// add fluid pressure
	s.xx() -= pt.m_pa;
	s.yy() -= pt.m_pa;
	s.zz() -= pt.m_pa;
	
	return s;
}

//-----------------------------------------------------------------------------
//! The tangent is the elastic tangent. Note
//! that this function is declared in the base class, so you don't have to 
//! reimplement it unless additional tangent components are required.

tens4ds FEBiphasicSolute::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& ept = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
	// call solid tangent routine
	tens4ds C = m_pSolid->Tangent(mp);
	
	// relative volume
	double J = ept.m_J;
	
	// fluid pressure and solute concentration
	double p = ppt.m_pa;
	double c = spt.m_c[0];
	
	// solubility and its derivative w.r.t. strain
	double kappa = m_pSolute->m_pSolub->Solubility(mp);
	double dkdJ = m_pSolute->m_pSolub->Tangent_Solubility_Strain(mp);
	
	// osmotic coefficient and its derivative w.r.t. strain
	double osmc = m_pOsmC->OsmoticCoefficient(mp);
	double dodJ = m_pOsmC->Tangent_OsmoticCoefficient_Strain(mp);
	
	double dp = m_Rgas*m_Tabs*c*J*(dodJ*kappa+osmc*dkdJ);
	
	// adjust tangent for pressures
	double D[6][6] = {0};
	C.extract(D);
	
	D[0][0] -= -p + dp;
	D[1][1] -= -p + dp;
	D[2][2] -= -p + dp;
	
	D[0][1] -= p + dp; D[1][0] -= p + dp;
	D[1][2] -= p + dp; D[2][1] -= p + dp;
	D[0][2] -= p + dp; D[2][0] -= p + dp;
	
	D[3][3] -= -p;
	D[4][4] -= -p;
	D[5][5] -= -p;
	
	return tens4ds(D);
}

//-----------------------------------------------------------------------------
//! Calculate fluid flux

vec3d FEBiphasicSolute::FluidFlux(FEMaterialPoint& pt)
{
	FEBiphasicMaterialPoint& ppt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	
	// fluid volume fraction (porosity) in current configuration
	double phiw = Porosity(pt);
	
	// pressure gradient
	vec3d gradp = ppt.m_gradp;
	
	// concentration
	double c = spt.m_c[0];
	
	// concentration gradient
	vec3d gradc = spt.m_gradc[0];
	
	// hydraulic permeability
	mat3ds kt = m_pPerm->Permeability(pt);
	
	// solute diffusivity in mixture
	mat3ds D = m_pSolute->m_pDiff->Diffusivity(pt);
	
	// solute free diffusivity
	double D0 = m_pSolute->m_pDiff->Free_Diffusivity(pt);
	
	// solubility
	double kappa = m_pSolute->m_pSolub->Solubility(pt);
	
	// identity matrix
	mat3dd I(1);
	
	// effective hydraulic permeability
	mat3ds ke = kt.inverse() + (I-D/D0)*(m_Rgas*m_Tabs*kappa*c/phiw/D0);
	ke = ke.inverse();
	
	// fluid flux w
	vec3d w = -(ke*(gradp + (D*gradc)*(m_Rgas*m_Tabs*kappa/D0)));
	
	return w;
}

//-----------------------------------------------------------------------------
//! Calculate solute molar flux

vec3d FEBiphasicSolute::SoluteFlux(FEMaterialPoint& pt)
{
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	
	// fluid volume fraction (porosity) in current configuration
	double phiw = Porosity(pt);
	
	// concentration
	double c = spt.m_c[0];
	
	// concentration gradient
	vec3d gradc = spt.m_gradc[0];
	
	// solute diffusivity in mixture
	mat3ds D = m_pSolute->m_pDiff->Diffusivity(pt);
	
	// solute free diffusivity
	double D0 = m_pSolute->m_pDiff->Free_Diffusivity(pt);
	
	// solubility
	double kappa = m_pSolute->m_pSolub->Solubility(pt);
	
	// fluid flux w
	vec3d w = FluidFlux(pt);
	
	// solute flux j
	vec3d j = D*(w*(c/D0) - gradc*phiw)*kappa;
	
	return j;
}

//-----------------------------------------------------------------------------
//! actual fluid pressure
double FEBiphasicSolute::Pressure(FEMaterialPoint& pt)
{
	FEBiphasicMaterialPoint& ppt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	
	// effective pressure
	double p = ppt.m_p;
	
	// effective concentration
	double c = spt.m_c[0];
	
	// osmotic coefficient
	double osmc = m_pOsmC->OsmoticCoefficient(pt);
	
	// solubility
	double kappa = m_pSolute->m_pSolub->Solubility(pt);
	
	// actual pressure
	double pa = p + m_Rgas*m_Tabs*osmc*kappa*c;
	
	return pa;
}

//-----------------------------------------------------------------------------
//! actual concentration
double FEBiphasicSolute::Concentration(FEMaterialPoint& pt)
{
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	
	// solubility
	double kappa = m_pSolute->m_pSolub->Solubility(pt);
	
	// actual concentration = solubility * effective concentration
	double ca = kappa*spt.m_c[0];
	
	return ca;
}

//-----------------------------------------------------------------------------
//! referential solute concentration
double FEBiphasicSolute::ReferentialConcentration(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& ept = *pt.ExtractData<FEElasticMaterialPoint>();

	double J = ept.m_J;
	double phiw = Porosity(pt);
	double cr = J*phiw*Concentration(pt);
	
	return cr;
}

//-----------------------------------------------------------------------------
//! partition coefficients and their derivatives
void FEBiphasicSolute::PartitionCoefficientFunctions(FEMaterialPoint& mp, double& kappa,
                                                double& dkdJ, double& dkdc)
{
    // evaluate the solubility and its derivatives w.r.t. J and c
    kappa = m_pSolute->m_pSolub->Solubility(mp);
    dkdJ = m_pSolute->m_pSolub->Tangent_Solubility_Strain(mp);
    dkdc = m_pSolute->m_pSolub->Tangent_Solubility_Concentration(mp,0);
}

double FEBiphasicSolute::GetReferentialFixedChargeDensity(const FEMaterialPoint& mp)
{
	const FEElasticMaterialPoint* ept = (mp.ExtractData<FEElasticMaterialPoint >());
	const FEBiphasicMaterialPoint* bpt = (mp.ExtractData<FEBiphasicMaterialPoint>());
	const FESolutesMaterialPoint* spt = (mp.ExtractData<FESolutesMaterialPoint >());
	double cf = (ept->m_J - bpt->m_phi0t) * spt->m_cF / (1 - bpt->m_phi0);
	return cf;
}
