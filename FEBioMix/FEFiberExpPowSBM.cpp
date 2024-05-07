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
#include "FEFiberExpPowSBM.h"
#include "FEMultiphasic.h"
#include <FECore/log.h>
#include <FEBioMech/FEElasticMixture.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEFiberExpPowSBM, FEElasticMaterial)
	ADD_PARAMETER(m_alpha, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha");
	ADD_PARAMETER(m_beta , FE_RANGE_GREATER_OR_EQUAL(2.0), "beta" );
	ADD_PARAMETER(m_ksi0 , FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi0" )->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_rho0 , FE_RANGE_GREATER_OR_EQUAL(0.0), "rho0" )->setUnits(UNIT_DENSITY);
	ADD_PARAMETER(m_g    , FE_RANGE_GREATER_OR_EQUAL(0.0), "gamma");
    ADD_PARAMETER(m_sbm , "sbm")->setEnums("$(sbms)");

    ADD_PROPERTY(m_fiber, "fiber");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// FEFiberExpPow
//-----------------------------------------------------------------------------

bool FEFiberExpPowSBM::Init()
{
	if (FEElasticMaterial::Init() == false) return false;

    // get the parent material which must be a multiphasic material
    FEMultiphasic* pMP = dynamic_cast<FEMultiphasic*> (GetAncestor());
	if (pMP == 0) {
		feLogError("Parent material must be multiphasic");
		return false;
	}
    
    // extract the local id of the SBM whose density controls Young's modulus from the global id
    m_lsbm = pMP->FindLocalSBMID(m_sbm);
	if (m_lsbm == -1) {
		feLogError("Invalid value for sbm");
		return false;
	}
    
    FEElasticMaterial* pem = pMP->GetSolid();
    FEElasticMixture* psm = dynamic_cast<FEElasticMixture*>(pem);
    if (psm == nullptr) m_comp = -1;    // in case material is not a solid mixture
    for (int i=0; i<psm->Materials(); ++i) {
        pem = psm->GetMaterial(i);
        if (pem == this) {
            m_comp = i;
            break;
        }
    }
    
	return true;
}

//-----------------------------------------------------------------------------
//! Create material point data
FEMaterialPointData* FEFiberExpPowSBM::CreateMaterialPointData()
{
    return new FERemodelingMaterialPoint(new FEElasticMaterialPoint);
}

//-----------------------------------------------------------------------------
//! update specialize material point data
void FEFiberExpPowSBM::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp)
{
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    FERemodelingMaterialPoint* rpt = mp.ExtractData<FERemodelingMaterialPoint>();
    rpt->m_rhor = spt.m_sbmr[m_lsbm];
    rpt->m_rhorp = spt.m_sbmrp[m_lsbm];
    rpt->m_sed = StrainEnergyDensity(mp);
}

//-----------------------------------------------------------------------------
mat3ds FEFiberExpPowSBM::Stress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    
    // initialize material constants
    double rhor = spt.m_sbmr[m_lsbm];
    double ksi = FiberModulus(mp, rhor);
    
    // deformation gradient
    mat3d &F = pt.m_F;
    double J = pt.m_J;
    
    // loop over all integration points
    vec3d nt;
    double In_1, Wl;
    const double eps = 0;
    mat3ds C = pt.RightCauchyGreen();
    mat3ds s;
    
    // get the material coordinate system
    mat3d Q = GetLocalCS(mp);

    // get local fiber direction
    vec3d fiber = m_fiber->unitVector(mp);

    // convert to global coordinates
    vec3d n0 = Q*fiber;
    
    // Calculate In = n0*C*n0
    In_1 = n0*(C*n0) - 1.0;
    
    // only take fibers in tension into consideration
    if (In_1 >= eps)
    {
        double alpha = m_alpha(mp);
        double beta = m_beta(mp);
        
        // get the global spatial fiber direction in current configuration
        nt = F*n0;
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        
        // calculate strain energy derivative
        Wl = ksi*pow(In_1, beta-1.0)*exp(alpha*pow(In_1, beta));
        
        // calculate the fiber stress
        s = N*(2.0*Wl/J);
    }
    else
    {
        s.zero();
    }
    
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FEFiberExpPowSBM::Tangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    
    // initialize material constants
    double rhor = spt.m_sbmr[m_lsbm];
    double ksi = FiberModulus(mp, rhor);
    
    // deformation gradient
    mat3d &F = pt.m_F;
    double J = pt.m_J;
    
    // loop over all integration points
    vec3d nt;
    double In_1, Wll;
    const double eps = 0;
    mat3ds C = pt.RightCauchyGreen();
    tens4ds c;
    
    // get the material coordinate system
    mat3d Q = GetLocalCS(mp);

    // get local fiber direction
    vec3d fiber = m_fiber->unitVector(mp);

    // convert to global coordinates
    vec3d n0 = Q*fiber;
    
    // Calculate In = n0*C*n0
    In_1 = n0*(C*n0) - 1.0;
    
    // only take fibers in tension into consideration
    if (In_1 >= eps)
    {
        double alpha = m_alpha(mp);
        double beta = m_beta(mp);
        
        // get the global spatial fiber direction in current configuration
        nt = F*n0;
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        tens4ds NxN = dyad1s(N);
        
        // calculate strain energy 2nd derivative
        double tmp = alpha*pow(In_1, beta);
        Wll = ksi*pow(In_1, beta-2.0)*((tmp+1)*beta-1.0)*exp(tmp);
        
        // calculate the fiber tangent
        c = NxN*(4.0*Wll/J);
    }
    else
    {
        c.zero();
    }
    
    return c;
}

//-----------------------------------------------------------------------------
double FEFiberExpPowSBM::StrainEnergyDensity(FEMaterialPoint& mp)
{
    double sed = 0.0;
    
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    
    // initialize material constants
    double rhor = spt.m_sbmr[m_lsbm];
    double ksi = FiberModulus(mp, rhor);
    
    // loop over all integration points
    double In_1;
    const double eps = 0;
    mat3ds C = pt.RightCauchyGreen();
    
    // get the material coordinate system
    mat3d Q = GetLocalCS(mp);

    // get local fiber direction
    vec3d fiber = m_fiber->unitVector(mp);

    // convert to global coordinates
    vec3d n0 = Q*fiber;
    
    // Calculate In = n0*C*n0
    In_1 = n0*(C*n0) - 1.0;
    
    // only take fibers in tension into consideration
    if (In_1 >= eps)
    {
        double alpha = m_alpha(mp);
        double beta = m_beta(mp);
        
        // calculate strain energy derivative
        if (alpha > 0) {
            sed = ksi/(alpha*beta)*(exp(alpha*pow(In_1, beta))-1);
        }
        else
            sed = ksi/beta*pow(In_1, beta);
    }
    
    return sed;
}

//-----------------------------------------------------------------------------
//! evaluate referential mass density
double FEFiberExpPowSBM::Density(FEMaterialPoint& pt)
{
    FERemodelingMaterialPoint* rpt = pt.ExtractData<FERemodelingMaterialPoint>();
    if (rpt) return rpt->m_rhor;
    else {
        FEElasticMixtureMaterialPoint* emp = pt.ExtractData<FEElasticMixtureMaterialPoint>();
        if (emp) {
            rpt = emp->GetPointData(m_comp)->ExtractData<FERemodelingMaterialPoint>();
            if (rpt) return rpt->m_rhor;
        }
    }
    return 0.0;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FEFiberExpPowSBM::StrainEnergy(FEMaterialPoint& mp)
{
    return StrainEnergyDensity(mp);
}

//-----------------------------------------------------------------------------
//! calculate tangent of strain energy density with mass density
double FEFiberExpPowSBM::Tangent_SE_Density(FEMaterialPoint& mp)
{
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    double rhor = spt.m_sbmr[m_lsbm];
    return StrainEnergy(mp)*m_g(mp)/rhor;
}

//-----------------------------------------------------------------------------
//! calculate tangent of stress with mass density
mat3ds FEFiberExpPowSBM::Tangent_Stress_Density(FEMaterialPoint& mp)
{
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    double rhor = spt.m_sbmr[m_lsbm];
    return Stress(mp)*m_g(mp)/rhor;
}

