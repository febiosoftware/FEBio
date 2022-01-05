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
#include "FEReactiveViscoelastic.h"
#include "FEElasticMixture.h"
#include "FEFiberMaterialPoint.h"
#include "FEScaledElasticMaterial.h"
#include <FECore/FECoreKernel.h>
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
    ADD_PARAMETER(m_emin , FE_RANGE_GREATER_OR_EQUAL(0.0), "emin");

	// set material properties
	ADD_PROPERTY(m_pBase, "elastic");
	ADD_PROPERTY(m_pBond, "bond", FEProperty::Optional);
    ADD_PROPERTY(m_scale, "scale", FEProperty::Optional);
	ADD_PROPERTY(m_pRelx, "relaxation");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEReactiveViscoelasticMaterial::FEReactiveViscoelasticMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_wmin = 0;
    m_btype = 0;
    m_ttype = 0;
    m_emin = 0;
    
    m_nmax = 0;

	m_pBase = nullptr;
	m_pBond = nullptr;
    m_scale = nullptr;
	m_pRelx = nullptr;
    
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
    
    // check number of elastic mixtures -- only one allowed, otherwise FEBio
    // does not know which FEElasticMixtureMaterialPoint to access
    int nmix = 0;
    if (dynamic_cast<FEElasticMixture*>(GetParent())) nmix++;
    if (dynamic_cast<FEElasticMixture*>(m_pBase)) nmix++;
    if (dynamic_cast<FEElasticMixture*>(m_pBond)) nmix++;
    
    if (nmix > 1) {
        feLogError("Parent, Elastic, and Bond materials of reactive viscoelastic material cannot include more than one elastic mixture");
        return false;
    }
    
    if (m_pBond == nullptr) {
        if (m_scale == nullptr) {
            feLogError("Either a bond material or a scale factor must be provided in a reactive viscoelastic material");
            return false;
        }
        else {
            m_pBond = new FEScaledElasticMaterial(GetFEModel(),m_pBase,m_scale);
            assert(m_pBase);
        }
    }
    
    if (!m_pBase->Init()) return false;
    if (!m_pBond->Init()) return false;
    if (!m_pRelx->Init()) return false;

    return FEElasticMaterial::Init();
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPoint* FEReactiveViscoelasticMaterial::CreateMaterialPointData()
{
	return new FEReactiveVEMaterialPoint(m_pBase->CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
//! detect new generation
bool FEReactiveViscoelasticMaterial::NewGeneration(FEMaterialPoint& mp)
{
    double d;
    double eps = max(m_emin, 10*std::numeric_limits<double>::epsilon());

    // get the elastic material point data
    FEElasticMaterialPoint& pe = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
    // check if the current deformation gradient is different from that of
    // the last generation, in which case store the current state
    // evaluate the relative deformation gradient
    mat3d F = pe.m_F;
    int lg = (int)pt.m_Uv.size() - 1;
    mat3ds Ui = (lg > -1) ? pt.m_Uv[lg].inverse() : mat3dd(1);
    mat3d Fu = F*Ui;

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
    double tv = time - pt.m_v[ig];

    switch (m_btype) {
        case 1:
        {
            if (tv >= 0)
                w = pt.m_f[ig]*m_pRelx->Relaxation(mp, tv, D);
        }
            break;
        case 2:
        {
            if (ig == 0) {
                w = m_pRelx->Relaxation(mp, tv, D);
            }
            else
            {
                double tu = time - pt.m_v[ig-1];
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
    int ng = (int)pt.m_Uv.size();
    
    double w = 1;
    
    for (int ig=0; ig<ng-1; ++ig)
    {
        // evaluate deformation gradient when this generation starts breaking
        ep.m_F = F;
        ep.m_J = J;
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
//! Stress function for strong bonds
mat3ds FEReactiveViscoelasticMaterial::StressStrongBonds(FEMaterialPoint& mp)
{
    return m_pBase->Stress(mp);
}

//-----------------------------------------------------------------------------
//! Stress function for weak bonds
mat3ds FEReactiveViscoelasticMaterial::StressWeakBonds(FEMaterialPoint& mp)
{
    double dt = GetFEModel()->GetTime().timeIncrement;
    if (dt == 0) return mat3ds(0, 0, 0, 0, 0, 0);
    
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
    // get fiber material point data (if it exists)
    FEFiberMaterialPoint* fp = mp.ExtractData<FEFiberMaterialPoint>();
    
    mat3ds D = ep.RateOfDeformation();
    
    // calculate the base material Cauchy stress
    mat3ds s; s.zero();
    
    // current number of breaking generations
    int ng = (int)pt.m_Uv.size();
    
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
            // evaluate bond mass fraction for this generation
            ep.m_F = F;
            ep.m_J = J;
            w = BreakingBondMassFraction(mp, ig, D);
            // evaluate relative deformation gradient for this generation
            if (ig > 0) {
                ep.m_F = F*pt.m_Uv[ig-1].inverse();
                ep.m_J = J/pt.m_Jv[ig-1];
                if (fp) fp->SetPreStretch(pt.m_Uv[ig-1]);
            }
            else {
                ep.m_F = F;
                ep.m_J = J;
                if (fp) fp->ResetPreStretch();
            }
            // evaluate bond stress
            sb = m_pBond->Stress(mp);
            // add bond stress to total stress
            s += (ig > 0) ? sb*w/pt.m_Jv[ig-1] : sb*w;
        }
        
        // restore safe copy of deformation gradient
        ep.m_F = F;
        ep.m_J = J;
    }
    
    return s;
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FEReactiveViscoelasticMaterial::Stress(FEMaterialPoint& mp)
{
    mat3ds s = StressStrongBonds(mp);
    s+= StressWeakBonds(mp);
    
	// return the total Cauchy stress
	return s;
}

//-----------------------------------------------------------------------------
//! Material tangent for strong bonds
tens4ds FEReactiveViscoelasticMaterial::TangentStrongBonds(FEMaterialPoint& mp)
{
    // calculate the base material tangent
    return m_pBase->Tangent(mp);
}

//-----------------------------------------------------------------------------
//! Material tangent for weak bonds
tens4ds FEReactiveViscoelasticMaterial::TangentWeakBonds(FEMaterialPoint& mp)
{
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
    // get fiber material point data (if it exists)
    FEFiberMaterialPoint* fp = mp.ExtractData<FEFiberMaterialPoint>();
    
    mat3ds D = ep.RateOfDeformation();
    
    // calculate the base material tangent
    tens4ds c; c.zero();
    
    // current number of breaking generations
    int ng = (int)pt.m_Uv.size();
    
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
            // evaluate bond mass fraction for this generation
            ep.m_F = F;
            ep.m_J = J;
            w = BreakingBondMassFraction(mp, ig, D);
            // evaluate relative deformation gradient for this generation
            if (ig > 0) {
                ep.m_F = F*pt.m_Uv[ig-1].inverse();
                ep.m_J = J/pt.m_Jv[ig-1];
                if (fp) fp->SetPreStretch(pt.m_Uv[ig-1]);
            }
            else {
                ep.m_F = F;
                ep.m_J = J;
                if (fp) fp->ResetPreStretch();
            }
            // evaluate bond tangent
            cb = m_pBond->Tangent(mp);
            // add bond tangent to total tangent
            c += (ig > 0) ? cb*w/pt.m_Jv[ig-1] : cb*w;
        }
        
        // restore safe copy of deformation gradient
        ep.m_F = F;
        ep.m_J = J;
    }
    
    return c;
}

//-----------------------------------------------------------------------------
//! Material tangent
tens4ds FEReactiveViscoelasticMaterial::Tangent(FEMaterialPoint& mp)
{
	tens4ds c = TangentStrongBonds(mp);
    c += TangentWeakBonds(mp);
    
	// return the total tangent
	return c;
}

//-----------------------------------------------------------------------------
//! strain energy density function in strong bonds
double FEReactiveViscoelasticMaterial::StrongBondSED(FEMaterialPoint& mp)
{
    // calculate the base material Cauchy stress
    return m_pBase->StrainEnergyDensity(mp);
}

//-----------------------------------------------------------------------------
//! strain energy density function in weak bonds
double FEReactiveViscoelasticMaterial::WeakBondSED(FEMaterialPoint& mp)
{
    double dt = GetFEModel()->GetTime().timeIncrement;
    if (dt == 0) return 0;
    
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
    // get fiber material point data (if it exists)
    FEFiberMaterialPoint* fp = mp.ExtractData<FEFiberMaterialPoint>();
    
    // get the viscous point data
    mat3ds D = ep.RateOfDeformation();
    
    double sed = 0;
    
    // current number of breaking generations
    int ng = (int)pt.m_Uv.size();
    
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
            // evaluate bond mass fraction for this generation
            ep.m_F = F;
            ep.m_J = J;
            w = BreakingBondMassFraction(mp, ig, D);
            // evaluate relative deformation gradient for this generation
            if (ig > 0) {
                ep.m_F = F*pt.m_Uv[ig-1].inverse();
                ep.m_J = J/pt.m_Jv[ig-1];
                if (fp) fp->SetPreStretch(pt.m_Uv[ig-1]);
            }
            else {
                ep.m_F = F;
                ep.m_J = J;
                if (fp) fp->ResetPreStretch();
            }
            // evaluate bond strain energy density
            sedb = m_pBond->StrainEnergyDensity(mp);
            // add bond strain energy density to total
            sed += sedb*w;
        }
        
        // restore safe copy of deformation gradient
        ep.m_F = F;
        ep.m_J = J;
    }
    
    return sed;
}

//-----------------------------------------------------------------------------
//! strain energy density function
double FEReactiveViscoelasticMaterial::StrainEnergyDensity(FEMaterialPoint& mp)
{
    double sed = StrongBondSED(mp);
    sed += WeakBondSED(mp);
    
    // return the total strain energy density
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
    
    int ng = (int)pt.m_v.size();
    m_nmax = max(m_nmax, ng);
    
    // don't cull if we have too few generations
    if (ng < 3) return;

    // don't reduce number of generations to less than max value achieved so far
    if (ng < m_nmax) return;

    // always check oldest generation
    double w0 = BreakingBondMassFraction(mp, 0, D);
    if (w0 < m_wmin) {
        double w1 = BreakingBondMassFraction(mp, 1, D);
        pt.m_v[1] = (w0*pt.m_v[0] + w1*pt.m_v[1])/(w0+w1);
        pt.m_Uv[1] = (pt.m_Uv[0]*w0 + pt.m_Uv[1]*w1)/(w0+w1);
        pt.m_Jv[1] = pt.m_Uv[1].det();
        pt.m_f[1] = (w0*pt.m_f[0] + w1*pt.m_f[1])/(w0+w1);
        pt.m_Uv.pop_front();
        pt.m_Jv.pop_front();
        pt.m_v.pop_front();
        pt.m_f.pop_front();
    }

    return;
}

//-----------------------------------------------------------------------------
//! Update specialized material points
void FEReactiveViscoelasticMaterial::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp)
{
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();

    FEElasticMaterialPoint& pe = *pt.ExtractData<FEElasticMaterialPoint>();
    
    mat3ds Uv = pe.RightStretch();
    double Jv = pe.m_J;

    // if new generation not already created for current time, check if it should
    if (pt.m_v.empty() || (pt.m_v.back() < tp.currentTime)) {
        // check if the current deformation gradient is different from that of
        // the last generation, in which case store the current state
        if (NewGeneration(mp)) {
            pt.m_v.push_back(tp.currentTime);
            pt.m_Uv.push_back(Uv);
            pt.m_Jv.push_back(Jv);
            double f = (!pt.m_v.empty()) ? ReformingBondMassFraction(mp) : 1;
            pt.m_f.push_back(f);
            CullGenerations(mp);
        }
    }
    // otherwise, if we already have a generation for the current time, update the stored values
    else if (pt.m_v.back() == tp.currentTime) {
        pt.m_Uv.back() = Uv;
        pt.m_Jv.back() = Jv;
        pt.m_f.back() = ReformingBondMassFraction(mp);
    }
}
