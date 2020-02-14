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
#include "FEReactivePlasticityHardening.h"
#include "FEDamageCriterion.h"
#include "FEElasticMaterial.h"
#include "FEUncoupledMaterial.h"
#include "FECore/FECoreKernel.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <FECore/matrix.h>

//////////////////////// PLASTICITY MATERIAL  /////////////////////////////////
// define the material parameters
BEGIN_FECORE_CLASS(FEReactivePlasticityHardening, FEElasticMaterial)
    // set material properties
    ADD_PROPERTY(m_pBase, "elastic");
    ADD_PROPERTY(m_pCrit, "criterion");

    ADD_PARAMETER(m_Ymin   , FE_RANGE_GREATER_OR_EQUAL(0.0), "ymin"  );
    ADD_PARAMETER(m_Ymax   , FE_RANGE_GREATER_OR_EQUAL(0.0), "ymax"  );
    ADD_PARAMETER(m_H      , FE_RANGE_GREATER_OR_EQUAL(0.0), "H"     );
    ADD_PARAMETER(m_delta  , FE_RANGE_GREATER_OR_EQUAL(0)  , "delta" );
    ADD_PARAMETER(m_rtol   , FE_RANGE_GREATER_OR_EQUAL(0.0), "rtol"  );
    ADD_PARAMETER(m_isochrc, "isochoric");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEReactivePlasticityHardening::FEReactivePlasticityHardening(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_Ymin = m_Ymax = 0;
    m_H = 0;
    m_delta = 0;
    m_isochrc = true;
    m_rtol = 1e-4;
    m_pBase = 0;
    m_pCrit = 0;
}

//-----------------------------------------------------------------------------
//! Initialization.
bool FEReactivePlasticityHardening::Init()
{
    FEUncoupledMaterial* m_pMat = dynamic_cast<FEUncoupledMaterial*>((FEElasticMaterial*)m_pBase);
    if (m_pMat != nullptr) {
        feLogError("Elastic material should not be of type uncoupled");
        return false;
    }
    
    return FEElasticMaterial::Init();
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPoint* FEReactivePlasticityHardening::CreateMaterialPointData()
{
    return new FEReactivePlasticityMaterialPoint(m_pBase->CreateMaterialPointData(), this);
}

//-----------------------------------------------------------------------------
//! evaluate elastic deformation gradient
void FEReactivePlasticityHardening::ElasticDeformationGradient(FEMaterialPoint& pt)
{
    // extract total deformation gradient
    FEElasticMaterialPoint& pe = *pt.ExtractData<FEElasticMaterialPoint>();
    // extract inverse of plastic deformation gradient and evaluate elastic deformation gradient
    FEReactivePlasticityMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticityMaterialPoint>();
    
    // evaluate trial def grad Fe
    mat3d Fs = pe.m_F;
    mat3d Fe = Fs*pp.m_Fusi[0];

    // store safe copy of total deformation gradient
    mat3d Ftmp = pe.m_F;
    double Jtmp = pe.m_J;
    pe.m_F = Fe; pe.m_J = Fe.det();
    
    // evaluate yield measure
    pp.m_Kv[0] = m_pCrit->DamageCriterion(pt);
    
    // restore total deformation gradient
    pe.m_F = Ftmp; pe.m_J = Jtmp;

    // get cumulative plastic strain
    double g = pp.m_gc[0];
    double Y = m_Ymin + m_H*g + (m_Ymax - m_Ymin)*(1-exp(-m_delta*g));


    // if there is no yielding, we're done
    double phi = pp.m_Kv[0] - Y;
    if (phi <= m_rtol*Y) {
        pp.m_Fvsi[0] = pp.m_Fusi[0];

        // evaluate octahedral plastic strain
        OctahedralPlasticStrain(pt);
        ReactiveHeatSupplyDensity(pt);
        return;
    }
    
    if ((pp.m_Kv[0] > pp.m_Ku[0]) && (pp.m_Ku[0] < Y*(1+m_rtol)))
        pp.m_w[0] = 1;
    
    // find Fv
    bool conv = false;
    int iter = 0;
    double lam = 0;
    mat3d Fv = Fe;
    Ftmp = pe.m_F;  // store safe copy
    Jtmp = pe.m_J;
    pe.m_F = Fv; pe.m_J = Fv.det();
    mat3ds Nv = YieldSurfaceNormal(pe);
    double Nvmag = Nv.norm();
    mat3dd I(1);
    double beta = 1;
    mat3ds ImN = I;
    double phi0=0, phi1=0, phi2=0, lam1=0, lam2=0, a, b, c=0, d;
    while (!conv) {
        ++iter;
        pe.m_F = Fv; pe.m_J = Fv.det();
        pp.m_Kv[0] = m_pCrit->DamageCriterion(pt);
        phi = pp.m_Kv[0] - Y;    // phi = 0 => stay on yield surface
        if (iter == 1) {
            phi0 = phi;
            c = phi0;
        }
        else if (iter == 2) {
            phi1 = phi;
            lam1 = lam;
        }
        else if (iter == 3) {
            phi2 = phi;
            lam2 = lam;
        }
        mat3d Rv = Fv*pe.RightStretchInverse();
        mat3d dFvdlam = -Fe*Nv*(beta/Nvmag);
        if (m_isochrc)
            dFvdlam += Fe*ImN*((ImN.inverse()*Nv/Nvmag).trace()*beta/3.);
        double dlam = -phi/(Rv*Nv*dFvdlam.transpose()).trace();
        lam += dlam;
        if (iter == 3) {
            d = lam1*lam2*(lam1-lam2);
            if (d == 0) {
                lam = (lam1*lam2 == 0) ? 0 : lam2;
            }
            else {
                a = (lam2*(phi1-phi0)-lam1*(phi2-phi0))/d;
                b = ((phi2-phi0)*lam1*lam1-(phi1-phi0)*lam2*lam2)/d;
                d = b*b - 4*a*c;
                if (d >= 0) {
                    if (a != 0) {
                        lam1 = (-b+sqrt(d))/(2*a);
                        lam2 = (-b-sqrt(d))/(2*a);
                        lam = (fabs(lam1) < fabs(lam2)) ? lam1 : lam2;
                    }
                    else if (b != 0) lam = -c/b;
                    else lam = 0;
                }
                else if (a != 0) {
                    lam = -b/(2*a);
                }
                else
                    lam = 0;
            }
            conv = true;
        }
        ImN = I - Nv*(lam/Nvmag);
        if (m_isochrc) beta = pow((pp.m_Fusi[0]*ImN).det(), -1./3.);
        Fv = Fe*ImN*beta;
        if (fabs(dlam) <= m_rtol*fabs(lam)) conv = true;
        if (fabs(lam) <= m_rtol*m_rtol) conv = true;
    }
    pe.m_F = Fv; pe.m_J = Fv.det();
    pp.m_Kv[0] = m_pCrit->DamageCriterion(pt);
    pe.m_F = Ftmp; pe.m_J = Jtmp;
    pp.m_Fvsi[0] = Fs.inverse()*Fv;

    // evaluate octahedral plastic strain
    OctahedralPlasticStrain(pt);
    ReactiveHeatSupplyDensity(pt);

    return;
}

//-----------------------------------------------------------------------------
//! calculate stress at material point
mat3ds FEReactivePlasticityHardening::Stress(FEMaterialPoint& pt)
{
    ElasticDeformationGradient(pt);
    
    // extract elastic material point
    FEElasticMaterialPoint& pe = *pt.ExtractData<FEElasticMaterialPoint>();
    // extract plastic material point
    FEReactivePlasticityMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticityMaterialPoint>();
    
    mat3ds s = m_pBase->Stress(pt)*(1 - pp.YieldedBonds());
    
    if (pp.m_w[0] > 0) {
        // get the elastic deformation gradient
        mat3d Fv = pe.m_F*pp.m_Fvsi[0];
        
        // store safe copy of total deformation gradient
        mat3d Fs = pe.m_F; double Js = pe.m_J;
        pe.m_F = Fv; pe.m_J = Fv.det();
        
        // evaluate the stress using the elastic deformation gradient
        s += m_pBase->Stress(pt)*pp.m_w[0];
        
        // restore the original deformation gradient
        pe.m_F = Fs; pe.m_J = Js;
    }
    
    // return the stress
    return s;
}

//-----------------------------------------------------------------------------
//! calculate tangent stiffness at material point
tens4ds FEReactivePlasticityHardening::Tangent(FEMaterialPoint& pt)
{
    ElasticDeformationGradient(pt);
    
    // extract elastic material point
    FEElasticMaterialPoint& pe = *pt.ExtractData<FEElasticMaterialPoint>();
    // extract plastic material point
    FEReactivePlasticityMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticityMaterialPoint>();
    
    tens4ds c = m_pBase->Tangent(pt)*(1 - pp.YieldedBonds());
    
    if (pp.m_w[0] > 0) {
        // get the elastic deformation gradient
        mat3d Fv = pe.m_F*pp.m_Fvsi[0];
        
        // store safe copy of total deformation gradient
        mat3d Fs = pe.m_F; double Js = pe.m_J;
        pe.m_F = Fv; pe.m_J = Fv.det();
        
        // evaluate the tangent using the elastic deformation gradient
        c += m_pBase->Tangent(pt)*pp.m_w[0];
        
        // restore the original deformation gradient
        pe.m_F = Fs; pe.m_J = Js;
    }
    
    // return the tangent
    return c;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FEReactivePlasticityHardening::StrainEnergyDensity(FEMaterialPoint& pt)
{
    ElasticDeformationGradient(pt);
    
    // extract elastic material point
    FEElasticMaterialPoint& pe = *pt.ExtractData<FEElasticMaterialPoint>();
    // extract plastic material point
    FEReactivePlasticityMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticityMaterialPoint>();
    
    double sed = m_pBase->StrainEnergyDensity(pt)*(1 - pp.YieldedBonds());
    
    if (pp.m_w[0] > 0) {
        // get the elastic deformation gradient
        mat3d Fv = pe.m_F*pp.m_Fvsi[0];
        
        // store safe copy of total deformation gradient
        mat3d Fs = pe.m_F; double Js = pe.m_J;
        pe.m_F = Fv; pe.m_J = Fv.det();
        
        // evaluate the tangent using the elastic deformation gradient
        sed += m_pBase->StrainEnergyDensity(pt)*pp.m_w[0];
        
        // restore the original deformation gradient
        pe.m_F = Fs; pe.m_J = Js;
    }
    
    // return the sed
    return sed;
}

//-----------------------------------------------------------------------------
// get the yield surface normal
mat3ds FEReactivePlasticityHardening::YieldSurfaceNormal(FEElasticMaterialPoint& pe)
{
    mat3ds s = m_pBase->Stress(pe);
    tens4ds c = m_pBase->Tangent(pe);
    mat3ds dPhi = m_pCrit->CriterionStressTangent(pe);
    mat3d M = dPhi*s*2 - mat3dd((dPhi*s).trace()) + c.dot(dPhi);
    mat3ds Ui = pe.RightStretchInverse();
    mat3d R = pe.m_F*Ui;
    mat3ds N = (R.transpose()*M*R*Ui).sym();
    return N;
}

//-----------------------------------------------------------------------------
//! calculate stress at material point
void FEReactivePlasticityHardening::OctahedralPlasticStrain(FEMaterialPoint& pt)
{
    // extract plastic material point
    FEReactivePlasticityMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticityMaterialPoint>();
    
    double ev[3];
    mat3ds Cvsi = (pp.m_Fvsi[0].transpose()*pp.m_Fvsi[0]).sym();
    Cvsi.eigen2(ev);
    for (int j=0; j<3; ++j) ev[j] = 1./sqrt(ev[j]);
    pp.m_gp[0] = sqrt(2.)/3.*sqrt(pow(ev[0] - ev[1],2) + pow(ev[1] - ev[2],2) + pow(ev[2] - ev[0],2));
}

//-----------------------------------------------------------------------------
//! evaluate reactive heat supply at material point
void FEReactivePlasticityHardening::ReactiveHeatSupplyDensity(FEMaterialPoint& pt)
{
    double Rhat = 0;
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    // extract elastic material point
    FEElasticMaterialPoint& pe = *pt.ExtractData<FEElasticMaterialPoint>();
    // extract plastic material point
    FEReactivePlasticityMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticityMaterialPoint>();
    
    if (dt == 0) {
        pp.m_Rhat = 0;
        return;
    }
    
    // store safe copy of total deformation gradient
    mat3d Fs = pe.m_F; double Js = pe.m_J;
    
    if (pp.m_w[0] > 0) {
        // get the elastic deformation gradients
        mat3d Fu = Fs*pp.m_Fusi[0];
        
        // evaluate strain energy density in the absence of yielding
        pe.m_F = Fu; pe.m_J = Fu.det();
        
        // evaluate the tangent using the elastic deformation gradient
        Rhat += m_pBase->StrainEnergyDensity(pt)*pp.m_w[0];
        
        mat3d Fv = Fs*pp.m_Fvsi[0];
        
        // evaluate strain energy density in the absence of yielding
        pe.m_F = Fv; pe.m_J = Fv.det();
        
        // evaluate the tangent using the elastic deformation gradient
        Rhat -= m_pBase->StrainEnergyDensity(pt)*pp.m_w[0];
    }
    
    // get rate
    Rhat /= dt;
    
    // restore the original deformation gradient
    pe.m_F = Fs; pe.m_J = Js;

    // return the reactive heat supply
    pp.m_Rhat = Rhat;
}
