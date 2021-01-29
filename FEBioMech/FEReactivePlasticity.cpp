/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include "FEReactivePlasticity.h"
#include "FEDamageCriterion.h"
#include "FEElasticMaterial.h"
#include "FEUncoupledMaterial.h"
#include "FECore/FECoreKernel.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <FECore/matrix.h>

//////////////////////// PLASTICITY MATERIAL  /////////////////////////////////
// define the material parameters
BEGIN_FECORE_CLASS(FEReactivePlasticity, FEElasticMaterial)
    // set material properties
    ADD_PROPERTY(m_pBase, "elastic");
    ADD_PROPERTY(m_pCrit, "yield_criterion");

    ADD_PARAMETER(m_Ymin   , FE_RANGE_GREATER_OR_EQUAL(0.0), "Y0"  );
    ADD_PARAMETER(m_Ymax   , FE_RANGE_GREATER_OR_EQUAL(0.0), "Ymax");
    ADD_PARAMETER(m_wmin   , FE_RANGE_CLOSED(0.0, 1.0)     , "w0"  );
    ADD_PARAMETER(m_we     , FE_RANGE_CLOSED(0.0, 1.0)     , "we"  );
    ADD_PARAMETER(m_n      , FE_RANGE_GREATER(0)           , "nf"  );
    ADD_PARAMETER(m_bias   , FE_RANGE_LEFT_OPEN(0.0, 1.0)  , "r"   );
    ADD_PARAMETER(m_isochrc, "isochoric");
    ADD_PARAMETER(m_rtol   , FE_RANGE_GREATER_OR_EQUAL(0.0), "rtol");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEReactivePlasticity::FEReactivePlasticity(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_n = 1;
    m_wmin = m_wmax = 1;
    m_we = 0;
    m_Ymin = m_Ymax = 0;
    m_isochrc = true;
    m_rtol = 1e-4;
    m_pBase = 0;
    m_pCrit = 0;
    m_bias = 0.9;
	m_secant_tangent = true;
}

//-----------------------------------------------------------------------------
//! Initialization.
bool FEReactivePlasticity::Init()
{
    m_wmax = 1 - m_we;
    FEUncoupledMaterial* m_pMat = dynamic_cast<FEUncoupledMaterial*>((FEElasticMaterial*)m_pBase);
    if (m_pMat != nullptr) {
        feLogError("Elastic material should not be of type uncoupled");
        return false;
    }
    if (m_wmax < m_wmin) {
        feLogError("wmax must be ≥ wmin");
        return false;
    }

    Ky.resize(m_n);
    w.resize(m_n);
    vector<double> Kp(m_n,0);
    
    if (m_n == 1) {
        Ky[0] = m_Ymin;
        w[0] = m_wmin;
    }
    else {
        // use bias r to reduce intervals in Ky and w as they increase proportionally
        double r = m_bias;
        // r= 1 uses uniform intervals
        if (r == 1) {
            w[0] = m_wmin;
            Kp[0] = m_Ymin;
            Ky[0] = Kp[0];
            double sw = w[0];
            for (int i=1; i<m_n; ++i) {
                w[i] = (m_wmax - m_wmin)/(m_n-1);
                Kp[i] = m_Ymin + (m_Ymax - m_Ymin)*i/(m_n-1);
                Ky[i] = Ky[i-1] + (Kp[i]-Kp[i-1])/(1-sw);
                sw += w[i];
            }
        }
        else {
            double c = (1-r)/(1-pow(r, m_n-1));
            w[0] = m_wmin;
            w[1] = c*(m_wmax-m_wmin);
            Kp[0] = m_Ymin;
            Kp[1] = Kp[0] + c*(m_Ymax - m_Ymin);
            double sw = w[0];
            Ky[0] = Kp[0];
            Ky[1] = Ky[0] + (Kp[1]-Kp[0])/(1-sw);
            sw += w[1];
            for (int i=2; i<m_n; ++i) {
                w[i] = w[i-1]*r;
                Kp[i] = Kp[i-1] + (Kp[i-1]-Kp[i-2])*r;
                Ky[i] = Ky[i-1] + (Kp[i]-Kp[i-1])/(1-sw);
                sw += w[i];
            }
        }
    }
    
    return FEElasticMaterial::Init();
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPoint* FEReactivePlasticity::CreateMaterialPointData()
{
    return new FEReactivePlasticityMaterialPoint(m_pBase->CreateMaterialPointData(), this);
}

//-----------------------------------------------------------------------------
//! evaluate elastic deformation gradient
void FEReactivePlasticity::ElasticDeformationGradient(FEMaterialPoint& pt)
{
    // extract total deformation gradient
    FEElasticMaterialPoint& pe = *pt.ExtractData<FEElasticMaterialPoint>();
    // extract inverse of plastic deformation gradient and evaluate elastic deformation gradient
    FEReactivePlasticityMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticityMaterialPoint>();

    for (int i=0; i<m_n; ++i) {
        mat3d Fs = pe.m_F;
        mat3d R = pe.m_F*pe.RightStretchInverse();
        mat3d Fe = Fs*pp.m_Fusi[i];

        // store safe copy of total deformation gradient
        mat3d Ftmp = pe.m_F;
        double Jtmp = pe.m_J;
        pe.m_F = Fe; pe.m_J = Fe.det();
        mat3ds Ue = pe.RightStretch();
        
        // evaluate yield measure
        pp.m_Kv[i] = m_pCrit->DamageCriterion(pt);
        
        // restore total deformation gradient
        pe.m_F = Ftmp; pe.m_J = Jtmp;
        
        // if there is no yielding, we're done
        double phi = pp.m_Kv[i] - Ky[i];
        if (phi <= m_rtol*Ky[i]) {
            pp.m_Fvsi[i] = pp.m_Fusi[i];
            continue;
        }
        
        if ((pp.m_Kv[i] > pp.m_Ku[i]) && (pp.m_Ku[i] < Ky[i]*(1+m_rtol)))
            pp.m_w[i] = w[i];
        
        // find Fv
        bool conv = false;
        int iter = 0;
        double lam = 0;
        mat3d Fv = Fe;
        Ftmp = pe.m_F;  // store safe copy
        Jtmp = pe.m_J;
        pe.m_F = Fv; pe.m_J = Fv.det();
        mat3ds Uv = pe.RightStretch();
        mat3ds Nv = YieldSurfaceNormal(pe);
        double Nvmag = Nv.norm();
        mat3dd I(1);
        double beta = 1;
        mat3ds ImN = I;
        double phi0=0, phi1=0, phi2=0, lam1=0, lam2=0, a, b, c=0, d;
        while (!conv) {
            ++iter;
            pe.m_F = Fv; pe.m_J = Fv.det();
            pp.m_Kv[i] = m_pCrit->DamageCriterion(pt);
            phi = pp.m_Kv[i] - Ky[i];    // phi = 0 => stay on yield surface
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
            mat3d dUvdlam = -Ue*Nv*(beta/Nvmag);
            if (m_isochrc)
                dUvdlam += Ue*ImN*((ImN.inverse()*Nv/Nvmag).trace()*beta/3.);
            double dlam = -phi/(Nv*dUvdlam.transpose()).trace();
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
            if (m_isochrc) beta = pow((pp.m_Fusi[i]*ImN).det(), -1./3.);
            Uv = (Ue*ImN).sym()*beta;
            Fv = R*Uv;
            if (fabs(dlam) <= m_rtol*fabs(lam)) conv = true;
            if (fabs(lam) <= m_rtol*m_rtol) conv = true;
        }
        pe.m_F = Fv; pe.m_J = Fv.det();
        pp.m_Kv[i] = m_pCrit->DamageCriterion(pt);
        pe.m_F = Ftmp; pe.m_J = Jtmp;
        pp.m_Fvsi[i] = Fs.inverse()*Fv;
    }

    // evaluate octahedral plastic strain
    OctahedralPlasticStrain(pt);
    ReactiveHeatSupplyDensity(pt);

    return;
}

//-----------------------------------------------------------------------------
//! calculate stress at material point
mat3ds FEReactivePlasticity::Stress(FEMaterialPoint& pt)
{
    ElasticDeformationGradient(pt);
    
    // extract elastic material point
    FEElasticMaterialPoint& pe = *pt.ExtractData<FEElasticMaterialPoint>();
    // extract plastic material point
    FEReactivePlasticityMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticityMaterialPoint>();
    
    mat3ds s = m_pBase->Stress(pt)*(1 - pp.YieldedBonds());
    
    for (int i=0; i<m_n; ++i) {
        if (pp.m_w[i] > 0) {
            // get the elastic deformation gradient
            mat3d Fv = pe.m_F*pp.m_Fvsi[i];
            
            // store safe copy of total deformation gradient
            mat3d Fs = pe.m_F; double Js = pe.m_J;
            pe.m_F = Fv; pe.m_J = Fv.det();
            
            // evaluate the stress using the elastic deformation gradient
            s += m_pBase->Stress(pt)*pp.m_w[i];
            
            // restore the original deformation gradient
            pe.m_F = Fs; pe.m_J = Js;
        }
    }
    
    // return the stress
    return s;
}

//-----------------------------------------------------------------------------
//! calculate tangent stiffness at material point
tens4ds FEReactivePlasticity::Tangent(FEMaterialPoint& pt)
{
    ElasticDeformationGradient(pt);
    
    // extract elastic material point
    FEElasticMaterialPoint& pe = *pt.ExtractData<FEElasticMaterialPoint>();
    // extract plastic material point
    FEReactivePlasticityMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticityMaterialPoint>();
    
    tens4ds c = m_pBase->Tangent(pt)*(1 - pp.YieldedBonds());
    
    for (int i=0; i<m_n; ++i) {
        if (pp.m_w[i] > 0) {
            // get the elastic deformation gradient
            mat3d Fv = pe.m_F*pp.m_Fvsi[i];
            
            // store safe copy of total deformation gradient
            mat3d Fs = pe.m_F; double Js = pe.m_J;
            pe.m_F = Fv; pe.m_J = Fv.det();
            
            // evaluate the tangent using the elastic deformation gradient
            c += m_pBase->Tangent(pt)*pp.m_w[i];
            
            // restore the original deformation gradient
            pe.m_F = Fs; pe.m_J = Js;
        }
    }
    
    // return the tangent
    return c;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FEReactivePlasticity::StrainEnergyDensity(FEMaterialPoint& pt)
{
    ElasticDeformationGradient(pt);
    
    // extract elastic material point
    FEElasticMaterialPoint& pe = *pt.ExtractData<FEElasticMaterialPoint>();
    // extract plastic material point
    FEReactivePlasticityMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticityMaterialPoint>();
    
    double sed = m_pBase->StrainEnergyDensity(pt)*(1 - pp.YieldedBonds());
    
    for (int i=0; i<m_n; ++i) {
        if (pp.m_w[i] > 0) {
            // get the elastic deformation gradient
            mat3d Fv = pe.m_F*pp.m_Fvsi[i];
            
            // store safe copy of total deformation gradient
            mat3d Fs = pe.m_F; double Js = pe.m_J;
            pe.m_F = Fv; pe.m_J = Fv.det();

            // evaluate the tangent using the elastic deformation gradient
            sed += m_pBase->StrainEnergyDensity(pt)*pp.m_w[i];
            
            // restore the original deformation gradient
            pe.m_F = Fs; pe.m_J = Js;
        }
    }
    
    // return the sed
    return sed;
}

//-----------------------------------------------------------------------------
// get the yield surface normal
mat3ds FEReactivePlasticity::YieldSurfaceNormal(FEElasticMaterialPoint& pe)
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
void FEReactivePlasticity::OctahedralPlasticStrain(FEMaterialPoint& pt)
{
    // extract plastic material point
    FEReactivePlasticityMaterialPoint& pp = *pt.ExtractData<FEReactivePlasticityMaterialPoint>();
    
    double ev[3];
    for (int i=0; i<m_n; ++i) {
        mat3ds Cvsi = (pp.m_Fvsi[i].transpose()*pp.m_Fvsi[i]).sym();
        Cvsi.eigen2(ev);
        for (int j=0; j<3; ++j) ev[j] = 1./sqrt(ev[j]);
        pp.m_gp[i] = sqrt(2.)/3.*sqrt(pow(ev[0] - ev[1],2) + pow(ev[1] - ev[2],2) + pow(ev[2] - ev[0],2));
    }
}

//-----------------------------------------------------------------------------
//! evaluate reactive heat supply at material point
void FEReactivePlasticity::ReactiveHeatSupplyDensity(FEMaterialPoint& pt)
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
    
    for (int i=0; i<m_n; ++i) {
        if (pp.m_w[i] > 0) {
            // get the elastic deformation gradients
            mat3d Fu = Fs*pp.m_Fusi[i];
            
            // evaluate strain energy density in the absence of yielding
            pe.m_F = Fu; pe.m_J = Fu.det();

            // evaluate the tangent using the elastic deformation gradient
            Rhat += m_pBase->StrainEnergyDensity(pt)*pp.m_w[i];
            
            mat3d Fv = Fs*pp.m_Fvsi[i];
            
            // evaluate strain energy density in the absence of yielding
            pe.m_F = Fv; pe.m_J = Fv.det();

            // evaluate the tangent using the elastic deformation gradient
            Rhat -= m_pBase->StrainEnergyDensity(pt)*pp.m_w[i];
        }
    }
    
    // get rate
    Rhat /= dt;
    
    // restore the original deformation gradient
    pe.m_F = Fs; pe.m_J = Js;

    // return the reactive heat supply
    pp.m_Rhat = Rhat;
}
