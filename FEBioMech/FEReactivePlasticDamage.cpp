/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2019 University of Utah, The Trustees of Columbia University in
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
#include "FEReactivePlasticDamage.h"
#include "FEReactivePlasticReturnMap.h"
#include "FEDamageCriterion.h"
#include "FEElasticMaterial.h"
#include "FEDamageCDF.h"
#include "FEUncoupledMaterial.h"
#include "FECore/FECoreKernel.h"
#include <FECore/FEMesh.h>
#include <FECore/log.h>
#include <FECore/matrix.h>

#ifndef max
#define max(a, b) ((a)>(b)?(a):(b))
#endif

//////////////////////// PLASTIC DAMAGE MATERIAL  /////////////////////////////////
// define the material parameters
BEGIN_FECORE_CLASS(FEReactivePlasticDamage, FEElasticMaterial)
    // set material properties
    ADD_PROPERTY(m_pBase,   "elastic");
    ADD_PROPERTY(m_pCrit,   "yield_criterion");
    ADD_PROPERTY(m_pFlow,   "flow_curve");
    ADD_PROPERTY(m_pYDamg,  "plastic_damage"          , FEProperty::Optional);
    ADD_PROPERTY(m_pYDCrit, "plastic_damage_criterion", FEProperty::Optional);
    ADD_PROPERTY(m_pIDamg,  "elastic_damage"          , FEProperty::Optional);
    ADD_PROPERTY(m_pIDCrit, "elastic_damage_criterion", FEProperty::Optional);

    ADD_PARAMETER(m_isochrc, "isochoric");
    ADD_PARAMETER(m_rtol   , FE_RANGE_GREATER_OR_EQUAL(0.0), "rtol");

    ADD_PARAMETER(m_secant_tangent, "secant_tangent");
    ADD_PARAMETER(m_formulation, "formulation")->setEnums("original\0symmetric\0");
    ADD_PARAMETER(m_mandel, "mandel_normal");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEReactivePlasticDamage::FEReactivePlasticDamage(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_isochrc = true;
    m_rtol = 1e-4;
    m_pBase = nullptr;
    m_pCrit = nullptr;
    m_pFlow = nullptr;
    m_pYDamg = nullptr;
    m_pYDCrit = nullptr;
    m_pIDamg = nullptr;
    m_pIDCrit = nullptr;
    m_secant_tangent = true;
    m_formulation = 0;
    m_mandel = false;
}

//-----------------------------------------------------------------------------
//! Initialization.
bool FEReactivePlasticDamage::Init()
{
    if (m_pFlow->Init() == false) return false;
    if (m_mandel && (m_formulation != 1)) {
        feLogError("mandel_normal requires formulation = symmetric");
        return false;
    }

    return FEElasticMaterial::Init();
}

//-----------------------------------------------------------------------------
void FEReactivePlasticDamage::Serialize(DumpStream& ar)
{
    FEElasticMaterial::Serialize(ar);
    ar & m_isochrc & m_rtol & m_formulation & m_mandel;
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPointData* FEReactivePlasticDamage::CreateMaterialPointData()
{
    FEMaterialPointData* ep = m_pBase->CreateMaterialPointData();
    FEMaterialPointData* fp = m_pFlow->CreateMaterialPointData();
    fp->SetNext(ep);
    return new FEReactivePlasticDamageMaterialPoint(fp, this);
}

//-----------------------------------------------------------------------------
//! evaluate elastic deformation gradient
void FEReactivePlasticDamage::ElasticDeformationGradient(FEMaterialPoint& pt)
{
    // initialize flow curve (if not done yet)
    if (m_pFlow->InitFlowCurve(pt)) {
        FEReactivePlasticDamageMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticDamageMaterialPoint>();
        pp.Init();
    }
    int n = (int)m_pFlow->BondFamilies(pt);
    
    // extract total deformation gradient
    FEElasticMaterialPoint& pe = *pt.ExtractData<FEElasticMaterialPoint>();
    // extract inverse of plastic deformation gradient and evaluate elastic deformation gradient
    FEReactivePlasticDamageMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticDamageMaterialPoint>();
    FEPlasticFlowCurveMaterialPoint& fp = *pt.ExtractData<FEPlasticFlowCurveMaterialPoint>();

    for (int i=0; i<n; ++i) {
        // total deformation gradient (for EAS and ANS shells, the shell domain sets
        // pe.m_F to the right stretch tensor of the enhanced strain before calling
        // this material, see FESolidMaterial::PK2Stress and MaterialTangent)
        mat3d Fs = pe.m_F;
        mat3d Fe = Fs*pp.m_Fusi[i];
        
        // store safe copy of total deformation gradient
        mat3d Ftmp = pe.m_F;
        double Jtmp = pe.m_J;
        pe.m_F = Fe; pe.m_J = Fe.det();

        // evaluate yield measure
        pp.m_Kv[i] = m_pCrit->DamageCriterion(pt);
        
        // restore total deformation gradient
        pe.m_F = Ftmp; pe.m_J = Jtmp;
        
        // if there is no yielding, we're done
        double phi = pp.m_Kv[i] - fp.m_Ky[i];
        if (phi <= m_rtol*fp.m_Ky[i]) {
            pp.m_Fvsi[i] = pp.m_Fusi[i];
            continue;
        }
        
        // check if i-th bond family is yielding
        if ((pp.m_Kv[i] > pp.m_Ku[i]) && (pp.m_Ku[i] < fp.m_Ky[i]*(1+m_rtol))) {
            if (pp.m_byldt[i] == false) {
                pp.m_byldt[i] = true;
                pp.m_wy[i] = (1.0-pp.m_di[i])*fp.m_w[i];
            }
        }
        // if not, and if this bond family has not yielded at previous times,
        // reset the mass fraction of yielded bonds to zero (in case m_wy[i] was
        // set to non-zero during a prior iteration at current time)
        else if (pp.m_byld[i] == false) {
            pp.m_byldt[i] = false;
            pp.m_wy[i] = 0;
        }

        // find Fv
        FEReactivePlasticReturnMap rm(m_pBase, m_pCrit, m_isochrc, m_rtol);
        FEReactivePlasticReturnMap::Result r = (m_formulation == 1) ?
            rm.Symmetric(pt, Fs, pp.m_Fusi[i], fp.m_Ky[i], m_mandel) :
            rm.Original(pt, Fs, pp.m_Fusi[i], fp.m_Ky[i]);
        if (!r.conv)
            feLogWarning("Plasticity iterations did not converge for bond family %d (phi/Ky = %lg)!\n",i,r.phi/fp.m_Ky[i]);
        // m_Kv holds the actual yield measure evaluated at the final Fv
        pp.m_Kv[i] = r.Kv;
        pp.m_Fvsi[i] = r.Fvsi;
    }
    
    // evaluate octahedral plastic strain
    OctahedralPlasticStrain(pt);
    ReactiveHeatSupplyDensity(pt);
    
    return;
}

//-----------------------------------------------------------------------------
// update plastic damage material point at each iteration
void FEReactivePlasticDamage::UpdateSpecializedMaterialPoints(FEMaterialPoint& pt, const FETimeInfo& tp)
{
    // initialize flow curve (if not done yet)
    if (m_pFlow->InitFlowCurve(pt)) {
        FEReactivePlasticDamageMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticDamageMaterialPoint>();
        pp.Init();
    }
    
    // extract total deformation gradient
    FEElasticMaterialPoint& pe = *pt.ExtractData<FEElasticMaterialPoint>();
    // extract plastic damage material point
    FEReactivePlasticDamageMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticDamageMaterialPoint>();
    FEPlasticFlowCurveMaterialPoint& fp = *pt.ExtractData<FEPlasticFlowCurveMaterialPoint>();
    
    // EAS and ANS shell domains evaluate stresses at the enhanced strain, but call this
    // function with the compatible deformation gradient; for consistency, evaluate the
    // damage at the right stretch tensor of the enhanced strain (as in PK2Stress)
    mat3d Fsafe = pe.m_F; double Jsafe = pe.m_J;
    FEShellElementNew* sel = dynamic_cast<FEShellElementNew*>(pt.m_elem);
    if (sel && (pt.m_index < (int)sel->m_E.size()) && (sel->m_E[pt.m_index].norm() > 0)) {
        mat3ds C = mat3dd(1) + sel->m_E[pt.m_index]*2;
        double eval[3];
        vec3d evec[3];
        C.eigen2(eval,evec);
        mat3ds U = dyad(evec[0])*sqrt(eval[0]) + dyad(evec[1])*sqrt(eval[1]) + dyad(evec[2])*sqrt(eval[2]);
        pe.m_F = U; pe.m_J = U.det();
    }
    
    ElasticDeformationGradient(pt);
    int n = (int)m_pFlow->BondFamilies(pt);

    // zero the total damage variable
    pp.m_D = 0.0;
    
    // get intact damage criterion
    if (m_pIDCrit) pp.m_Etrial = m_pIDCrit->DamageCriterion(pt);
    double Es = max(pp.m_Etrial, pp.m_Emax);
    
    for (int i=0; i<n; ++i) {
        if (pp.m_byldt[i] == false)
        {
            pp.m_di[i] = m_pIDamg ? m_pIDamg->cdf(pt,Es) : 0;
            pp.m_d[i] = pp.m_di[i]*fp.m_w[i];
            // what if we iterate here, update damage, then the next iteration decides we actually are yielding?
            // no mechanism to undo the extra damage we've added
            // do i need to save the final Eim for each family?
        }
        else
        {
            mat3d Fp = pp.m_Fvsi[i].inverse();
            
            // store safe copy of total deformation gradient
            mat3d Ftmp = pe.m_F;
            pe.m_F = Fp;
            
            // calculate damage
            if (m_pYDCrit) pp.m_Eyt[i] = m_pYDCrit->DamageCriterion(pt);
            double Ey = max(pp.m_Eyt[i], pp.m_Eym[i]);
            pp.m_dy[i] = m_pYDamg ? m_pYDamg->cdf(pt,Ey) : 0;
            
            // restore total deformation gradient
            pe.m_F = Ftmp;
            
            // calculate bond fractions
            pp.m_wy[i] = (1.0-pp.m_dy[i])*(1.0-pp.m_di[i])*fp.m_w[i];
            
            // calculate damage
            pp.m_d[i] = (pp.m_di[i]+pp.m_dy[i]*(1.0-pp.m_di[i]))*fp.m_w[i];
        }
        // sum the damage over all bond families
        pp.m_D += pp.m_d[i];
    }
    // add damage to persistent elastic bonds
    pp.m_di[n] = m_pIDamg ? m_pIDamg->cdf(pt,Es) : 0;
    pp.m_d[n] = pp.m_di[n]*fp.m_w[n];
    pp.m_D += pp.m_d[n];
    
    // restore the deformation gradient
    pe.m_F = Fsafe; pe.m_J = Jsafe;
}


//-----------------------------------------------------------------------------
//! calculate stress at material point
mat3ds FEReactivePlasticDamage::Stress(FEMaterialPoint& pt)
{
    ElasticDeformationGradient(pt);
    int n = (int)m_pFlow->BondFamilies(pt);

    // extract elastic material point
    FEElasticMaterialPoint& pe = *pt.ExtractData<FEElasticMaterialPoint>();
    // extract plastic damage material point
    FEReactivePlasticDamageMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticDamageMaterialPoint>();
    
    mat3ds s = m_pBase->Stress(pt)*pp.IntactBonds();
    
    for (int i=0; i<n; ++i) {
        // get the elastic deformation gradient
        mat3d Fv = pe.m_F*pp.m_Fvsi[i];
        
        // store safe copy of total deformation gradient
        mat3d Fs = pe.m_F; double Js = pe.m_J;
        pe.m_F = Fv; pe.m_J = Fv.det();
        
        // evaluate the damaged plastic stress using the elastic deformation gradient
        s += m_pBase->Stress(pt)*pp.m_wy[i];
        
        // restore the original deformation gradient
        pe.m_F = Fs; pe.m_J = Js;
    }
    
    // return the stress
    return s;
}

//-----------------------------------------------------------------------------
//! calculate tangent stiffness at material point
tens4ds FEReactivePlasticDamage::Tangent(FEMaterialPoint& pt)
{
    ElasticDeformationGradient(pt);
    int n = (int)m_pFlow->BondFamilies(pt);

    // extract elastic material point
    FEElasticMaterialPoint& pe = *pt.ExtractData<FEElasticMaterialPoint>();
    // extract plastic material point
    FEReactivePlasticDamageMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticDamageMaterialPoint>();
    
    tens4ds c = m_pBase->Tangent(pt)*pp.IntactBonds();
    
    for (int i=0; i<n; ++i) {
        // get the elastic deformation gradient
        mat3d Fv = pe.m_F*pp.m_Fvsi[i];
        
        // store safe copy of total deformation gradient
        mat3d Fs = pe.m_F; double Js = pe.m_J;
        pe.m_F = Fv; pe.m_J = Fv.det();
        
        // evaluate the tangent using the elastic deformation gradient
        c += m_pBase->Tangent(pt)*pp.m_wy[i];
        
        // restore the original deformation gradient
        pe.m_F = Fs; pe.m_J = Js;
    }
    
    // return the tangent
    return c;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FEReactivePlasticDamage::StrainEnergyDensity(FEMaterialPoint& pt)
{
    ElasticDeformationGradient(pt);
    int n = (int)m_pFlow->BondFamilies(pt);

    // extract elastic material point
    FEElasticMaterialPoint& pe = *pt.ExtractData<FEElasticMaterialPoint>();
    // extract plastic material point
    FEReactivePlasticDamageMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticDamageMaterialPoint>();
    
    double sed = m_pBase->StrainEnergyDensity(pt)*pp.IntactBonds();
    
    for (int i=0; i<n; ++i) {
        // get the elastic deformation gradient
        mat3d Fv = pe.m_F*pp.m_Fvsi[i];
        double Jvsi = m_isochrc ? 1 : pp.m_Fvsi[i].det();

        // store safe copy of total deformation gradient
        mat3d Fs = pe.m_F; double Js = pe.m_J;
        pe.m_F = Fv; pe.m_J = Fv.det();
        
        // evaluate the tangent using the elastic deformation gradient
        sed += m_pBase->StrainEnergyDensity(pt)*pp.m_wy[i]/Jvsi;
        
        // restore the original deformation gradient
        pe.m_F = Fs; pe.m_J = Js;
    }
    
    // return the sed
    return sed;
}

//-----------------------------------------------------------------------------
//! calculate damage at material point
double FEReactivePlasticDamage::Damage(FEMaterialPoint& pt, int k)
{
    //get plastic damage material point data
    FEReactivePlasticDamageMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticDamageMaterialPoint>();
    
    return pp.m_d[k];
}

//-----------------------------------------------------------------------------
// get the yield surface normal
mat3ds FEReactivePlasticDamage::YieldSurfaceNormal(FEMaterialPoint& mp)
{
    FEReactivePlasticReturnMap rm(m_pBase, m_pCrit, m_isochrc, m_rtol);
    return rm.YieldSurfaceNormal(mp);
}

//-----------------------------------------------------------------------------
//! calculate octahedral plastic strain at material point
void FEReactivePlasticDamage::OctahedralPlasticStrain(FEMaterialPoint& pt)
{
    int n = (int)m_pFlow->BondFamilies(pt);
    
    // extract plastic material point
    FEReactivePlasticDamageMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticDamageMaterialPoint>();
    
    for (int i=0; i<n; ++i) {
        // effective plastic strain sqrt(2/3 devE:devE), where E is the Lagrange strain
        // of the plastic deformation gradient Fvs (Zimmerman et al., JMPS 2021, Eq. 5.10);
        // this is the same measure as the "DC octahedral shear strain" criterion
        mat3d Fvs = pp.m_Fvsi[i].inverse();
        mat3ds Evs = ((Fvs.transpose()*Fvs).sym() - mat3dd(1))*0.5;
        mat3ds devE = Evs.dev();
        pp.m_gp[i] = sqrt(devE.dotdot(devE)*(2./3.));
    }
}

//-----------------------------------------------------------------------------
//! evaluate reactive heat supply at material point
void FEReactivePlasticDamage::ReactiveHeatSupplyDensity(FEMaterialPoint& pt)
{
    double Rhat = 0;
    
    double dt = CurrentTimeIncrement();
    
    // extract elastic material point
    FEElasticMaterialPoint& pe = *pt.ExtractData<FEElasticMaterialPoint>();
    // extract plastic material point
    FEReactivePlasticDamageMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticDamageMaterialPoint>();
    
    if (dt == 0) {
        pp.m_Rhat = 0;
        return;
    }
    
    // store safe copy of total deformation gradient
    mat3d Fs = pe.m_F; double Js = pe.m_J;
    
    int n = (int)m_pFlow->BondFamilies(pt);
    
    for (int i=0; i<n; ++i) {
        // get the elastic deformation gradients
        mat3d Fu = Fs*pp.m_Fusi[i];
        
        // evaluate strain energy density in the absence of yielding
        pe.m_F = Fu; pe.m_J = Fu.det();
        
        // evaluate the tangent using the elastic deformation gradient
        Rhat += m_pBase->StrainEnergyDensity(pt)*pp.m_wy[i];
        
        mat3d Fv = Fs*pp.m_Fvsi[i];
        
        // evaluate strain energy density in the absence of yielding
        pe.m_F = Fv; pe.m_J = Fv.det();
        
        // evaluate the tangent using the elastic deformation gradient
        Rhat -= m_pBase->StrainEnergyDensity(pt)*pp.m_wy[i];
    }
    
    // get rate
    Rhat /= dt;
    
    // restore the original deformation gradient
    pe.m_F = Fs; pe.m_J = Js;
    
    // return the reactive heat supply
    pp.m_Rhat = Rhat;
}

