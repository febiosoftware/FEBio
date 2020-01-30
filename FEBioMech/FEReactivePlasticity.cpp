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
    ADD_PROPERTY(m_pCrit, "criterion");

    ADD_PARAMETER(m_Ymin   , FE_RANGE_GREATER_OR_EQUAL(0.0), "ymin"  );
    ADD_PARAMETER(m_Ymax   , FE_RANGE_GREATER_OR_EQUAL(0.0), "ymax"  );
    ADD_PARAMETER(m_wmin   , FE_RANGE_GREATER_OR_EQUAL(0.0), "wmin"  );
    ADD_PARAMETER(m_n      , FE_RANGE_GREATER_OR_EQUAL(0)  , "n"     );
    ADD_PARAMETER(m_itmax  , FE_RANGE_GREATER_OR_EQUAL(0)  , "maxiter");
    ADD_PARAMETER(m_isochrc, "isochoric"   );

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEReactivePlasticity::FEReactivePlasticity(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_n = 1;
    m_wmin = 1;
    m_Ymin = m_Ymax = 0;
    m_isochrc = true;
    m_itmax = 10;
    m_pBase = 0;
    m_pCrit = 0;
}

//-----------------------------------------------------------------------------
//! Initialization.
bool FEReactivePlasticity::Init()
{
    FEUncoupledMaterial* m_pMat = dynamic_cast<FEUncoupledMaterial*>((FEElasticMaterial*)m_pBase);
    if (m_pMat != nullptr) {
        feLogError("Elastic material should not be of type uncoupled");
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
        w[0] = m_wmin;
        Kp[0] = m_Ymin;
        Ky[0] = Kp[0];
        double sw = w[0];
        for (int i=1; i<m_n; ++i) {
            w[i] = (1 - m_wmin)/(m_n-1);
            Kp[i] = m_Ymin + (m_Ymax - m_Ymin)*i/(m_n-1);
            Ky[i] = Ky[i-1] + (Kp[i]-Kp[i-1])/(1-sw);
            sw += w[i];
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

    double eps = 1e-4;  // convergence criterion for iterative solution of lambda

    for (int i=0; i<m_n; ++i) {
        mat3d Fs = pe.m_F;
        mat3d Fu = Fs*pp.m_Fusi[i];
        
        // store safe copy of total deformation gradient
        mat3d Ftmp = pe.m_F;
        double Jtmp = pe.m_J;
        pe.m_F = Fu; pe.m_J = Fu.det();
        
        // evaluate yield measure
        pp.m_Kv[i] = m_pCrit->DamageCriterion(pt);
        
        // restore total deformation gradient
        pe.m_F = Ftmp; pe.m_J = Jtmp;
        
        // interpolate to approximate Fe on yield surface
        double alpha = 0;
        
        // if there is no yielding, we're done
        double phi = pp.m_Kv[i] - Ky[i];
        if (phi <= 0) return;
        
        if ((pp.m_Kv[i] > pp.m_Ku[i]) && (pp.m_Ku[i] < Ky[i])) {
            alpha = (Ky[i] - pp.m_Ku[i])/(pp.m_Kv[i] - pp.m_Ku[i]);
            pp.m_w[i] = w[i];
        }
        
        // evaluate elastic deformation gradient at previous time
        mat3d Fup = pp.m_Fp*pp.m_Fusi[i];
        
        // interpolate to approximate Fu on yield surface
        Fu = Fu*alpha + Fup*(1-alpha);  // Fu(v) is on yield surface
        
        // find Fv
        bool conv = false;
        int iter = 0;
        double lam = 0;
        mat3d Fe = Fs*pp.m_Fusi[i];
        mat3d Fv = Fe;
        Ftmp = pe.m_F;  // store safe copy
        Jtmp = pe.m_J;
        pe.m_F = Fv; pe.m_J = Fv.det();
        mat3ds Nv = YieldSurfaceNormal(pe);
        double Nvmag = Nv.norm();
        while (!conv) {
            if (++iter > m_itmax) break;
            pe.m_F = Fv; pe.m_J = Fv.det();
            pp.m_Kv[i] = m_pCrit->DamageCriterion(pt);
            phi = pp.m_Kv[i] - Ky[i];    // phi = 0 => stay on yield surface
            double dlam = phi*Nvmag/(Nv*Fe*Nv).trace();
            lam += dlam;
            Fv = Fe*(mat3dd(1) - Nv*(lam/Nvmag));
            if (fabs(dlam) <= eps*fabs(lam)) conv = true;
            if (fabs(phi) <= eps*eps*Ky[i]) conv = true;
        }
        pe.m_F = Ftmp; pe.m_J = Jtmp;
        if (m_isochrc) Fv = Fv*pow(Fs.det()/Fv.det(),1./3.);
        if (iter > m_itmax) feLogError("Max number of iterations exceeded in reactive plasticity solver.");
        else pp.m_Fvsi[i] = Fs.inverse()*Fv;

        // evaluate octahedral plastic shear strain
        double ev[3];
        mat3ds Cvsi = (pp.m_Fvsi[i].transpose()*pp.m_Fvsi[i]).sym();
        Cvsi.eigen2(ev);
        for (int j=0; j<3; ++j) ev[j] = 1./sqrt(ev[j]);
        pp.m_gp[i] = sqrt(2.)/3.*sqrt(pow(ev[0] - ev[1],2) + pow(ev[1] - ev[2],2) + pow(ev[2] - ev[0],2));
    }
    
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
        // get the elastic deformation gradient
        mat3d Fe = pe.m_F*pp.m_Fvsi[i];
        
        // store safe copy of total deformation gradient
        mat3d Ftmp = pe.m_F;
        pe.m_F = Fe;
        
        // evaluate the stress using the elastic deformation gradient
        s += m_pBase->Stress(pt)*pp.m_w[i];
        
        // restore the original deformation gradient
        pe.m_F = Ftmp;
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
        // get the elastic deformation gradient
        mat3d Fe = pe.m_F*pp.m_Fvsi[i];
        
        // store safe copy of total deformation gradient
        mat3d Ftmp = pe.m_F;
        pe.m_F = Fe;
        
        // evaluate the tangent using the elastic deformation gradient
        c += m_pBase->Tangent(pt)*pp.m_w[i];
        
        // restore the original deformation gradient
        pe.m_F = Ftmp;
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
        // get the elastic deformation gradient
        mat3d Fe = pe.m_F*pp.m_Fvsi[i];
        
        // store safe copy of total deformation gradient
        mat3d Ftmp = pe.m_F;
        pe.m_F = Fe;
        
        // evaluate the tangent using the elastic deformation gradient
        sed += m_pBase->StrainEnergyDensity(pt)*pp.m_w[i];
        
        // restore the original deformation gradient
        pe.m_F = Ftmp;
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

