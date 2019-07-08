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

//////////////////////// PLASTICITY MATERIAL  /////////////////////////////////
// define the material parameters
BEGIN_FECORE_CLASS(FEReactivePlasticity, FEElasticMaterial)
ADD_PARAMETER(m_Ymin   , FE_RANGE_GREATER_OR_EQUAL(0.0), "ymin"  );
ADD_PARAMETER(m_Ymax   , FE_RANGE_GREATER_OR_EQUAL(0.0), "ymax"  );
ADD_PARAMETER(m_wmin   , FE_RANGE_GREATER_OR_EQUAL(0.0), "wmin"  );
ADD_PARAMETER(m_n      , FE_RANGE_GREATER_OR_EQUAL(0)  , "n"     );
ADD_PARAMETER(m_stretch, "stretch"     );
ADD_PARAMETER(m_isochrc, "isochoric"   );

// set material properties
ADD_PROPERTY(m_pBase, "elastic");
ADD_PROPERTY(m_pCrit, "criterion");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEReactivePlasticity::FEReactivePlasticity(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_n = 1;
    m_wmin = 1;
    m_Ymin = m_Ymax = 0;
    m_stretch = true;
    m_isochrc = true;
    
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
    
    for (int i=0; i<m_n; ++i) {
        mat3d Fe = pe.m_F*pp.m_Fi[i];
        
        // store safe copy of total deformation gradient
        mat3d Ftmp = pe.m_F;
        pe.m_F = Fe;
        
        // evaluate yield measure
        pp.m_Kt[i] = m_pCrit->DamageCriterion(pt);
        
        // restore total deformation gradient
        pe.m_F = Ftmp;
        
        // interpolate to approximate Fe on yield surface
        double alpha = 0;
        
        // if there is no yielding, we're done
        if (pp.m_Kt[i] < Ky[i]) return;
        
        if ((pp.m_Kt[i] > pp.m_Kp[i]) && (pp.m_Kp[i] < Ky[i])) {
            alpha = (Ky[i] - pp.m_Kp[i])/(pp.m_Kt[i] - pp.m_Kp[i]);
            pp.m_w[i] = w[i];
        }
        
        // evaluate elastic deformation gradient at previous time
        mat3d Fep = pp.m_Fp*pp.m_Fi[i];
        
        // interpolate to approximate Fe on yield surface
        mat3d Fa = Fe*alpha + Fep*(1-alpha);
        double Ja = Fa.det();
        
        // evaluate inverse of trial plasticity deformation gradient
        pp.m_Ft[i] = pe.m_F.inverse()*Fa;
        if (m_isochrc) pp.m_Ft[i] *= pow(pe.m_J/Ja,1./3.);
        
        // force the trial plasticity deformation gradient to be a pure stretch
        mat3ds Ci = (pp.m_Ft[i]*pp.m_Ft[i].transpose()).sym();
        double lam[3];
        vec3d v[3];
        Ci.eigen(lam,v);
        if (m_stretch) pp.m_Ft[i] = (v[0] & v[0])*sqrt(lam[0]) + (v[1] & v[1])*sqrt(lam[1]) + (v[2] & v[2])*sqrt(lam[2]);
        
        // evaluate octahedral plastic shear strain
        for (int j=0; j<3; ++j) lam[j] = 1./sqrt(lam[j]);
        pp.m_gp[i] = sqrt(2.)/3.*sqrt(pow(lam[0] - lam[1],2) + pow(lam[1] - lam[2],2) + pow(lam[2] - lam[0],2));
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
        mat3d Fe = pe.m_F*pp.m_Ft[i];
        
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
        mat3d Fe = pe.m_F*pp.m_Ft[i];
        
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
        mat3d Fe = pe.m_F*pp.m_Ft[i];
        
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
