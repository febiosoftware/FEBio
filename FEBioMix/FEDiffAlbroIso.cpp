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
#include "FEDiffAlbroIso.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEDiffAlbroIso, FESoluteDiffusivity)
	ADD_PARAMETER(m_diff0 , FE_RANGE_GREATER_OR_EQUAL(0.0), "free_diff")->setUnits(UNIT_DIFFUSIVITY)->setLongName("free diffusivity");
	ADD_PARAMETER(m_cdinv , FE_RANGE_GREATER_OR_EQUAL(0.0), "cdinv"    )->setUnits("L^3/n");
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
    
	if (m_lsol == -1) {
		feLogError("Invalid value for sol"); 
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Free diffusivity
double FEDiffAlbroIso::Free_Diffusivity(FEMaterialPoint& mp)
{
    FESoluteInterface* psm = dynamic_cast<FESoluteInterface*>(GetAncestor());
	
    // solute concentration
    double ca = psm->GetActualSoluteConcentration(mp, m_lsol);
    
    // diffusivity coefficient
    double d = m_diff0*exp(-m_cdinv*ca);
    
	return d;
}

//-----------------------------------------------------------------------------
//! Tangent of free diffusivity with respect to concentration
double FEDiffAlbroIso::Tangent_Free_Diffusivity_Concentration(FEMaterialPoint& mp, const int isol)
{
    FESoluteInterface* psm = dynamic_cast<FESoluteInterface*>(GetAncestor());

    // solute concentration
    double ca = psm->GetActualSoluteConcentration(mp, m_lsol);
    double c = psm->GetEffectiveSoluteConcentration(mp, m_lsol);
    double dkdc = psm->dkdc(mp, m_lsol, isol);
    double k = psm->GetPartitionCoefficient(mp, m_lsol);
    
    // diffusivity coefficient
    double d = m_diff0*exp(-m_cdinv*ca);
    // derivative of d w.r.t. actual concentration
    double dc = -m_cdinv*d;    
    
    // tangent w.r.t. concentration
    if (isol == m_lsol)
        return dc*(k+dkdc*c);
    else
        return dc*dkdc*c;
}

//-----------------------------------------------------------------------------
//! Diffusivity tensor.
mat3ds FEDiffAlbroIso::Diffusivity(FEMaterialPoint& mp)
{
    FESoluteInterface* psm = dynamic_cast<FESoluteInterface*>(GetAncestor());
    FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());

	FEElasticMaterialPoint* et = mp.ExtractData<FEElasticMaterialPoint>();
	
	// relative volume
	double J = et->m_J;
	
	// solid volume fraction in reference configuration
	double phi0 = pbm->GetReferentialSolidVolumeFraction(mp);

    // porosity in current configuration
    double phiw = 1 - phi0/J;
    // solute concentration
    double ca = psm->GetActualSoluteConcentration(mp, m_lsol);
    
    // diffusivity coefficient
    double d = m_diff0*exp(-m_alphad*(1-phiw)/phiw - m_cdinv*ca);
	
	// diffusivity tensor
    mat3dd dt(d);
	
	return dt;
}

//-----------------------------------------------------------------------------
//! Tangent of diffusivity with respect to strain
tens4dmm FEDiffAlbroIso::Tangent_Diffusivity_Strain(FEMaterialPoint &mp)
{
    FESoluteInterface* psm = dynamic_cast<FESoluteInterface*>(GetAncestor());
    FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());

    FEElasticMaterialPoint* et = mp.ExtractData<FEElasticMaterialPoint>();
	
	// Identity
	mat3dd I(1);
	
	// relative volume
	double J = et->m_J;
	
	// solid volume fraction in reference configuration
    double phi0 = pbm->GetReferentialSolidVolumeFraction(mp);

    // porosity in current configuration
    double phiw = 1 - phi0/J;
    // solute concentration
    double ca = psm->GetActualSoluteConcentration(mp, m_lsol);
    double c = psm->GetEffectiveSoluteConcentration(mp, m_lsol);
    double dkdJ = psm->dkdJ(mp, m_lsol);
    
    // diffusivity coefficient
    double d = m_diff0*exp(-m_alphad*(1-phiw)/phiw - m_cdinv*ca);
    
    // derivative of (J d) w.r.t. J
    double dJ = d*(1+J*(m_alphad*phi0/(J-phi0)/(J-phi0) - m_cdinv*c*dkdJ));
		
	tens4dmm D4 = dyad1s(I)*dJ-dyad4s(I)*(2*d);
	
	return D4;
}

//-----------------------------------------------------------------------------
//! Tangent of diffusivity with respect to concentration
mat3ds FEDiffAlbroIso::Tangent_Diffusivity_Concentration(FEMaterialPoint &mp, const int isol)
{
    FESoluteInterface* psm = dynamic_cast<FESoluteInterface*>(GetAncestor());
    FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());

    FEElasticMaterialPoint* et = mp.ExtractData<FEElasticMaterialPoint>();
	
	// relative volume
	double J = et->m_J;
	
	// solid volume fraction in reference configuration
    double phi0 = pbm->GetReferentialSolidVolumeFraction(mp);

    // porosity in current configuration
    double phiw = 1 - phi0/J;
    // solute concentration
    double ca = psm->GetActualSoluteConcentration(mp, m_lsol);
    double c = psm->GetEffectiveSoluteConcentration(mp, m_lsol);
    double dkdc = psm->dkdc(mp, m_lsol, isol);
    double k = psm->GetPartitionCoefficient(mp, m_lsol);
    
    // diffusivity coefficient
    double d = m_diff0*exp(-m_alphad*(1-phiw)/phiw - m_cdinv*ca);
    // derivative of d w.r.t. actual concentration
    double dc = -m_cdinv*d;
    
    // tangent w.r.t. concentration
    if (isol == m_lsol) {
        return mat3dd(dc*(k+dkdc*c));
    } else
        return mat3dd(dc*dkdc*c);
}
