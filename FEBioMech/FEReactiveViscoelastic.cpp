/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FEReactiveViscoelastic.h"
#include "FECore/FECoreKernel.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <limits>

///////////////////////////////////////////////////////////////////////////////
//
// FEReactiveViscoelasticMaterial
//
///////////////////////////////////////////////////////////////////////////////

// Material parameters for the FEMultiphasic material
BEGIN_FECORE_CLASS(FEReactiveViscoelasticMaterial, FEElasticMaterial)
    ADD_PARAMETER(m_wmin , FE_RANGE_CLOSED(0.0, 1.0), "wmin");
    ADD_PARAMETER(m_btype, FE_RANGE_CLOSED(1,2), "kinetics");
    ADD_PARAMETER(m_ttype, FE_RANGE_CLOSED(0,2), "trigger");

	// set material properties
	ADD_PROPERTY(m_pBase, "elastic");
	ADD_PROPERTY(m_pBond, "bond");
	ADD_PROPERTY(m_pRelx, "relaxation");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEReactiveViscoelasticMaterial::FEReactiveViscoelasticMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_wmin = 0;
    m_btype = 0;
    m_ttype = 0;

	m_pBase = 0;
	m_pBond = 0;
	m_pRelx = 0;
}

//-----------------------------------------------------------------------------
//! data initialization
bool FEReactiveViscoelasticMaterial::Init()
{
    FEUncoupledMaterial* m_pMat = dynamic_cast<FEUncoupledMaterial*>((FEElasticMaterial*)m_pBase);
	if (m_pMat != nullptr) {
		feLogError("Elastic material should not be of type uncoupled");
		return false;
	}
    
    m_pMat = dynamic_cast<FEUncoupledMaterial*>((FEElasticMaterial*)m_pBond);
	if (m_pMat != nullptr) {
		feLogError("Bond material should not be of type uncoupled");
		return false;
	}
    
    return FEElasticMaterial::Init();
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPoint* FEReactiveViscoelasticMaterial::CreateMaterialPointData()
{
	return new FEReactiveVEMaterialPoint(m_pBase->CreateMaterialPointData(), this);
}

//-----------------------------------------------------------------------------
//! detect new generation
bool FEReactiveViscoelasticMaterial::NewGeneration(FEMaterialPoint& mp)
{
    double d;
    double eps = std::numeric_limits<double>::epsilon();

    // get the elastic material poit data
    FEElasticMaterialPoint& pe = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
    // check if the current deformation gradient is different from that of
    // the last generation, in which case store the current state
    // evaluate the relative deformation gradient
    mat3d F = pe.m_F;
    int lg = (int)pt.m_Fi.size() - 1;
    mat3d Fi = (lg > -1) ? pt.m_Fi[lg] : mat3d(mat3dd(1));
    mat3d Fu = F*Fi;

    switch (m_ttype) {
        case 0:
        {
            // trigger in response to any strain
            // evaluate the Lagrangian strain
            mat3ds E = ((Fu.transpose()*Fu).sym() - mat3dd(1))/2;
            
            d = E.norm();
        }
            break;
        case 1:
        {
            // trigger in response to distortional strain
            // evaluate spatial Hencky (logarithmic) strain
            mat3ds Bu = (Fu*Fu.transpose()).sym();
            double l[3];
            vec3d v[3];
            Bu.eigen2(l,v);
            mat3ds h = (dyad(v[0])*log(l[0]) + dyad(v[1])*log(l[1]) + dyad(v[2])*log(l[2]))/2;
            
            // evaluate distortion magnitude (always positive)
            d = (h.dev()).norm();
        }
            break;
        case 2:
        {
            // trigger in response to dilatational strain
            d = fabs(log(Fu.det()));
        }
            break;
            
        default:
            d = 0;
            break;
    }
    
    if (d > eps) return true;
    
    return false;
}

//-----------------------------------------------------------------------------
//! evaluate bond mass fraction
double FEReactiveViscoelasticMaterial::BreakingBondMassFraction(FEMaterialPoint& mp, const int ig, const mat3ds D)
{
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
    // bond mass fraction
    double w = 0;
    
    // current time
    double time = GetFEModel()->GetTime().currentTime;
    
    switch (m_btype) {
        case 1:
        {
            // time when this generation started breaking
            double v = pt.m_v[ig];
            
            if (time >= v)
                w = pt.m_w[ig]*m_pRelx->Relaxation(mp, time - v, D);
        }
            break;
        case 2:
        {
            double tu, tv;
            if (ig == 0) {
                tv = time - pt.m_v[ig];
                w = m_pRelx->Relaxation(mp, tv, D);
            }
            else
            {
                tu = time - pt.m_v[ig-1];
                tv = time - pt.m_v[ig];
                w = m_pRelx->Relaxation(mp, tv, D) - m_pRelx->Relaxation(mp, tu, D);
            }
        }
            break;
            
        default:
            break;
    }
    
    assert((w >= 0) && (w <= 1));
    
    return w;
}

//-----------------------------------------------------------------------------
//! evaluate bond mass fraction of reforming generation
double FEReactiveViscoelasticMaterial::ReformingBondMassFraction(FEMaterialPoint& mp)
{
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
    mat3ds D = ep.RateOfDeformation();
    
    // keep safe copy of deformation gradient
    mat3d F = ep.m_F;
    double J = ep.m_J;
    
    // get current number of generations
    int ng = (int)pt.m_Fi.size();
    
    double w = 1;
    
    for (int ig=0; ig<ng-1; ++ig)
    {
        // evaluate relative deformation gradient for this generation Fu(v)
        ep.m_F = pt.m_Fi[ig+1].inverse()*pt.m_Fi[ig];
        ep.m_J = pt.m_Ji[ig]/pt.m_Ji[ig+1];
        // evaluate the breaking bond mass fraction for this generation
        w -= BreakingBondMassFraction(mp, ig, D);
    }
    
    // restore safe copy of deformation gradient
    ep.m_F = F;
    ep.m_J = J;
    
    assert((w >= 0) && (w <= 1));
    
    // return the bond mass fraction of the reforming generation
    return w;
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FEReactiveViscoelasticMaterial::Stress(FEMaterialPoint& mp)
{
	double dt = GetFEModel()->GetTime().timeIncrement;
	if (dt == 0) return mat3ds(0, 0, 0, 0, 0, 0);
    
	// get the elastic part
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// get the reactive viscoelastic point data
	FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
    mat3ds D = ep.RateOfDeformation();
    
	// calculate the base material Cauchy stress
	mat3ds s = m_pBase->Stress(mp);
    
    // current number of breaking generations
    int ng = (int)pt.m_Fi.size();
    
    // no bonds have broken
    if (ng == 0) {
        s += m_pBond->Stress(mp);
    }
    // bonds have broken
    else {
        // keep safe copy of deformation gradient
        mat3d F = ep.m_F;
        double J = ep.m_J;
        
        double w;
        mat3ds sb;
        
        // calculate the bond stresses for breaking generations
        for (int ig=0; ig<ng; ++ig) {
            // evaluate relative deformation gradient for this generation
            ep.m_F = F*pt.m_Fi[ig];
            ep.m_J = J*pt.m_Ji[ig];
            // evaluate bond mass fraction for this generation
            w = BreakingBondMassFraction(mp, ig, D);
            // evaluate bond stress
            sb = m_pBond->Stress(mp);
            // add bond stress to total stress
            s += sb*(w*pt.m_Ji[ig]);
        }
        
        // restore safe copy of deformation gradient
        ep.m_F = F;
        ep.m_J = J;
    }
    
	// return the total Cauchy stress
	return s;
}

//-----------------------------------------------------------------------------
//! Material tangent
tens4ds FEReactiveViscoelasticMaterial::Tangent(FEMaterialPoint& mp)
{
    CullGenerations(mp);
    
	// get the elastic part
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// get the reactive viscoelastic point data
	FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
    mat3ds D = ep.RateOfDeformation();
    
	// calculate the base material tangent
	tens4ds c = m_pBase->Tangent(mp);
    
    // current number of breaking generations
    int ng = (int)pt.m_Fi.size();
    
    // no bonds have broken
    if (ng == 0) {
        c += m_pBond->Tangent(mp);
    }
    // bonds have broken
    else {
        // keep safe copy of deformation gradient
        mat3d F = ep.m_F;
        double J = ep.m_J;
        
        double w;
        tens4ds cb;
        
        // calculate the bond tangents for breaking generations
        for (int ig=0; ig<ng; ++ig) {
            // evaluate relative deformation gradient for this generation
            ep.m_F = F*pt.m_Fi[ig];
            ep.m_J = J*pt.m_Ji[ig];
            // evaluate bond mass fraction for this generation
            w = BreakingBondMassFraction(mp, ig, D);
            // evaluate bond tangent
            cb = m_pBond->Tangent(mp);
            // add bond tangent to total tangent
            c += cb*(w*pt.m_Ji[ig]);
        }
        
        // restore safe copy of deformation gradient
        ep.m_F = F;
        ep.m_J = J;
    }
    
	// return the total tangent
	return c;
}

//-----------------------------------------------------------------------------
//! strain energy density function
double FEReactiveViscoelasticMaterial::StrainEnergyDensity(FEMaterialPoint& mp)
{
	double dt = GetFEModel()->GetTime().timeIncrement;
	if (dt == 0) return 0;
    
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
    // get the viscous point data
    mat3ds D = ep.RateOfDeformation();
    
    // calculate the base material Cauchy stress
    double sed = m_pBase->StrainEnergyDensity(mp);
    
    // current number of breaking generations
    int ng = (int)pt.m_Fi.size();
    
    // no bonds have broken
    if (ng == 0) {
        sed += m_pBond->StrainEnergyDensity(mp);
    }
    // bonds have broken
    else {
        // keep safe copy of deformation gradient
        mat3d F = ep.m_F;
        double J = ep.m_J;
        
        double w;
        double sedb;
        
        // calculate the strain energy density for breaking generations
        for (int ig=0; ig<ng; ++ig) {
            // evaluate relative deformation gradient for this generation
            ep.m_F = F*pt.m_Fi[ig];
            ep.m_J = J*pt.m_Ji[ig];
            // evaluate bond mass fraction for this generation
            w = BreakingBondMassFraction(mp, ig, D);
            // evaluate bond stress
            sedb = m_pBond->StrainEnergyDensity(mp);
            // add bond stress to total stress
            sed += sedb*w;
        }
        
        // restore safe copy of deformation gradient
        ep.m_F = F;
        ep.m_J = J;
    }
    
    // return the total Cauchy stress
    return sed;
}

//-----------------------------------------------------------------------------
//! Cull generations that have relaxed below a threshold
void FEReactiveViscoelasticMaterial::CullGenerations(FEMaterialPoint& mp)
{
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
    mat3ds D = ep.RateOfDeformation();
    
    if (pt.m_Fi.empty()) return;

    // culling termination flag
    bool done = false;
    
    // always check oldest generation
    while (!done) {
        double w = BreakingBondMassFraction(mp, 0, D);
        if ((w > m_wmin) || (pt.m_Fi.size() == 1))
            done = true;
        else {
            pt.m_Fi.pop_front();
            pt.m_Ji.pop_front();
            pt.m_v.pop_front();
            pt.m_w.pop_front();
        }
    }
    
    return;
}
