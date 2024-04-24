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
#include "FEFiberPowLinearSBM.h"
#include "FEMultiphasic.h"
#include <FECore/log.h>
#include <FEBioMech/FEElasticMixture.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEFiberPowLinearSBM, FEElasticMaterial)
    ADD_PARAMETER(m_E0   , FE_RANGE_GREATER_OR_EQUAL(0.0), "E0"   )->setUnits(UNIT_PRESSURE);;
	ADD_PARAMETER(m_lam0 , FE_RANGE_GREATER         (1.0), "lam0" );
	ADD_PARAMETER(m_beta , FE_RANGE_GREATER_OR_EQUAL(2.0), "beta" );
    ADD_PARAMETER(m_rho0 , FE_RANGE_GREATER_OR_EQUAL(0.0), "rho0" )->setUnits(UNIT_DENSITY);;
	ADD_PARAMETER(m_g    , FE_RANGE_GREATER_OR_EQUAL(0.0), "gamma");
    ADD_PARAMETER(m_sbm , "sbm")->setEnums("$(sbms)");

    ADD_PROPERTY(m_fiber, "fiber");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// FEFiberPowLinear
//-----------------------------------------------------------------------------

bool FEFiberPowLinearSBM::Init()
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
FEMaterialPointData* FEFiberPowLinearSBM::CreateMaterialPointData()
{
    return new FERemodelingMaterialPoint(new FEElasticMaterialPoint);
}

//-----------------------------------------------------------------------------
//! update specialize material point data
void FEFiberPowLinearSBM::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp)
{
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    FERemodelingMaterialPoint* rpt = mp.ExtractData<FERemodelingMaterialPoint>();
    rpt->m_rhor = spt.m_sbmr[m_lsbm];
    rpt->m_rhorp = spt.m_sbmrp[m_lsbm];
    rpt->m_sed = StrainEnergyDensity(mp);
}

//-----------------------------------------------------------------------------
mat3ds FEFiberPowLinearSBM::Stress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    
    double lam0 = m_lam0(mp);
    double beta = m_beta(mp);
    
    // initialize material constants
    double rhor = spt.m_sbmr[m_lsbm];
    double E = FiberModulus(mp, rhor);
    double I0 = lam0*lam0;
    double ksi = E/4/(beta-1)*pow(I0, -3./2.)*pow(I0-1, 2-beta);
    double b = ksi*pow(I0-1, beta-1) + E/2/sqrt(I0);
    
    // deformation gradient
    mat3d &F = pt.m_F;
    double J = pt.m_J;
    
    // loop over all integration points
    vec3d nt;
    double In, sn;
    const double eps = 0;
    mat3ds C = pt.RightCauchyGreen();
    mat3ds s;
    
    // get the material coordinate system
    mat3d Q = GetLocalCS(mp);

    // get local fiber direction
    vec3d fiber = m_fiber->unitVector(mp);

    // convert to global coordinates
    vec3d n0 = Q*fiber;
    
    // Calculate In
    In = n0*(C*n0);
    
    // only take fibers in tension into consideration
    if (In - 1 >= eps)
    {
        // get the global spatial fiber direction in current configuration
        nt = F*n0/sqrt(In);
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        
        // calculate the fiber stress magnitude
        sn = (In < I0) ?
        2*In*ksi*pow(In-1, beta-1) :
        2*b*In - E*sqrt(In);
        
        // calculate the fiber stress
        s = N*(sn/J);
    }
    else
    {
        s.zero();
    }
    
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FEFiberPowLinearSBM::Tangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    
    double lam0 = m_lam0(mp);
    double beta = m_beta(mp);
    
    // initialize material constants
    double rhor = spt.m_sbmr[m_lsbm];
    double E = FiberModulus(mp, rhor);
    double I0 = lam0*lam0;
    double ksi = E/4/(beta-1)*pow(I0, -3./2.)*pow(I0-1, 2-beta);
    
    // deformation gradient
    mat3d &F = pt.m_F;
    double J = pt.m_J;
    
    // loop over all integration points
    vec3d nt;
    double In, cn;
    const double eps = 0;
    mat3ds C = pt.RightCauchyGreen();
    tens4ds c;
    
    // get the material coordinate system
    mat3d Q = GetLocalCS(mp);

    // get local fiber direction
    vec3d fiber = m_fiber->unitVector(mp);

    // convert to global coordinates
    vec3d n0 = Q*fiber;
    
    // Calculate In
    In = n0*(C*n0);
    
    // only take fibers in tension into consideration
    if (In - 1 >= eps)
    {
        // get the global spatial fiber direction in current configuration
        nt = F*n0/sqrt(In);
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        tens4ds NxN = dyad1s(N);
        
        // calculate modulus
        cn = (In < I0) ?
        4*In*In*ksi*(beta-1)*pow(In-1, beta-2) :
        E*sqrt(In);
        
        // calculate the fiber tangent
        c = NxN*(cn/J);
    }
    else
    {
        c.zero();
    }
    
    return c;
}

//-----------------------------------------------------------------------------
double FEFiberPowLinearSBM::StrainEnergyDensity(FEMaterialPoint& mp)
{
    double sed = 0.0;
    
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    
    double lam0 = m_lam0(mp);
    double beta = m_beta(mp);
    
    // initialize material constants
    double rhor = spt.m_sbmr[m_lsbm];
    double E = FiberModulus(mp, rhor);
    double I0 = lam0*lam0;
    double ksi = E/4/(beta-1)*pow(I0, -3./2.)*pow(I0-1, 2-beta);
    double b = ksi*pow(I0-1, beta-1) + E/2/sqrt(I0);
    
    // loop over all integration points
    double In;
    const double eps = 0;
    mat3ds C = pt.RightCauchyGreen();
    
    // get the material coordinate system
    mat3d Q = GetLocalCS(mp);

    // get local fiber direction
    vec3d fiber = m_fiber->unitVector(mp);

    // convert to global coordinates
    vec3d n0 = Q*fiber;
    
    // Calculate In = n0*C*n0
    In = n0*(C*n0);
    
    // only take fibers in tension into consideration
    if (In - 1 >= eps)
    {
        // calculate strain energy density
        sed = (In < I0) ?
        ksi/beta*pow(In-1, beta) :
        b*(In-I0) - E*(sqrt(In)-sqrt(I0)) + ksi/beta*pow(I0-1, beta);
    }
    
    return sed;
}

//-----------------------------------------------------------------------------
//! evaluate referential mass density
double FEFiberPowLinearSBM::Density(FEMaterialPoint& pt)
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
double FEFiberPowLinearSBM::StrainEnergy(FEMaterialPoint& mp)
{
    return StrainEnergyDensity(mp);
}

//-----------------------------------------------------------------------------
//! calculate tangent of strain energy density with mass density
double FEFiberPowLinearSBM::Tangent_SE_Density(FEMaterialPoint& mp)
{
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    double rhor = spt.m_sbmr[m_lsbm];
    return StrainEnergy(mp)*m_g(mp)/rhor;
}

//-----------------------------------------------------------------------------
//! calculate tangent of stress with mass density
mat3ds FEFiberPowLinearSBM::Tangent_Stress_Density(FEMaterialPoint& mp)
{
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    double rhor = spt.m_sbmr[m_lsbm];
    return Stress(mp)*m_g(mp)/rhor;
}

