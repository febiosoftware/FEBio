//
//  FEDiffAlbroIso.cpp
//  FEBio
//
//  Created by Gerard Ateshian on 12/17/13.
//  Copyright (c) 2013 febio.org. All rights reserved.
//

#include "FEDiffAlbroIso.h"
#include "FEBiphasicSolute.h"
#include "FETriphasic.h"
#include "FEMultiphasic.h"
#include "stdafx.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEDiffAlbroIso, FESoluteDiffusivity)
	ADD_PARAMETER(m_diff0 , FE_RANGE_GREATER_OR_EQUAL(0.0), "free_diff");
	ADD_PARAMETER(m_cdinv , FE_RANGE_GREATER_OR_EQUAL(0.0), "cdinv"    );
	ADD_PARAMETER(m_alphad, FE_RANGE_GREATER_OR_EQUAL(0.0), "alphad"   );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEDiffAlbroIso::FEDiffAlbroIso(FEModel* pfem) : FESoluteDiffusivity(pfem)
{
	m_diff0 = 1;
	m_cdinv = m_alphad = 0;
    m_lsol = -1;
}

//-----------------------------------------------------------------------------
//! Initialization.
bool FEDiffAlbroIso::Init()
{
	if (FESoluteDiffusivity::Init() == false) return false;

	// get the grandparent material which must be
    // a biphasic-solute/triphasic/multiphasic material
    FESolute* pSol = dynamic_cast<FESolute*> (GetParent());
    m_lsol = pSol->GetSoluteLocalID();
    
	if (m_lsol == -1) return fecore_error("Invalid value for sol");

	return true;
}

//-----------------------------------------------------------------------------
//! Free diffusivity
double FEDiffAlbroIso::Free_Diffusivity(FEMaterialPoint& mp)
{
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
    // solute concentration
    double ca = spt.m_ca[m_lsol];
    
    // diffusivity coefficient
    double d = m_diff0*exp(-m_cdinv*ca);
    
	return d;
}

//-----------------------------------------------------------------------------
//! Tangent of free diffusivity with respect to concentration
double FEDiffAlbroIso::Tangent_Free_Diffusivity_Concentration(FEMaterialPoint& mp, const int isol)
{
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
    // solute concentration
    double ca = spt.m_ca[m_lsol];
    double c = spt.m_c[m_lsol];
    
    // diffusivity coefficient
    double d = m_diff0*exp(-m_cdinv*ca);
    // derivative of d w.r.t. actual concentration
    double dc = -m_cdinv*d;
    double dkdc = spt.m_dkdc[m_lsol][isol];
    
    // tangent w.r.t. concentration
    if (isol == m_lsol) {
        double k = spt.m_k[m_lsol];
        return dc*(k+dkdc*c);
    } else
        return dc*dkdc*c;
}

//-----------------------------------------------------------------------------
//! Diffusivity tensor.
mat3ds FEDiffAlbroIso::Diffusivity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
	// relative volume
	double J = et.m_J;
	
	// solid volume fraction in reference configuration
	double phi0 = ppt.m_phi0;
    // porosity in current configuration
    double phiw = 1 - phi0/J;
    // solute concentration
    double ca = spt.m_ca[m_lsol];
    
    // diffusivity coefficient
    double d = m_diff0*exp(-m_alphad*(1-phiw)/phiw - m_cdinv*ca);
	
	// diffusivity tensor
    mat3dd dt(d);
	
	return dt;
}

//-----------------------------------------------------------------------------
//! Tangent of diffusivity with respect to strain
tens4ds FEDiffAlbroIso::Tangent_Diffusivity_Strain(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
	// Identity
	mat3dd I(1);
	
	// relative volume
	double J = et.m_J;
	
	// solid volume fraction in reference configuration
	double phi0 = ppt.m_phi0;
    // porosity in current configuration
    double phiw = 1 - phi0/J;
    // solute concentration
    double ca = spt.m_ca[m_lsol];
    double c = spt.m_c[m_lsol];
    double dkdJ = spt.m_dkdJ[m_lsol];
    
    // diffusivity coefficient
    double d = m_diff0*exp(-m_alphad*(1-phiw)/phiw - m_cdinv*ca);
    
    // derivative of (J d) w.r.t. J
    double dJ = d*(1+J*(m_alphad*phi0/(J-phi0)/(J-phi0) - m_cdinv*c*dkdJ));
		
	tens4ds D4 = dyad1s(I)*dJ-dyad4s(I)*(2*d);
	
	return D4;
}

//-----------------------------------------------------------------------------
//! Tangent of diffusivity with respect to concentration
mat3ds FEDiffAlbroIso::Tangent_Diffusivity_Concentration(FEMaterialPoint &mp, const int isol)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
	// relative volume
	double J = et.m_J;
	
	// solid volume fraction in reference configuration
	double phi0 = ppt.m_phi0;
    // porosity in current configuration
    double phiw = 1 - phi0/J;
    // solute concentration
    double ca = spt.m_ca[m_lsol];
    double c = spt.m_c[m_lsol];
    
    // diffusivity coefficient
    double d = m_diff0*exp(-m_alphad*(1-phiw)/phiw - m_cdinv*ca);
    // derivative of d w.r.t. actual concentration
    double dc = -m_cdinv*d;
    double dkdc = spt.m_dkdc[m_lsol][isol];
    
    // tangent w.r.t. concentration
    if (isol == m_lsol) {
        double k = spt.m_k[m_lsol];
        return mat3dd(dc*(k+dkdc*c));
    } else
        return mat3dd(dc*dkdc*c);
}
