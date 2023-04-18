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
	ADD_PROPERTY(m_pBond, "bond");
	ADD_PROPERTY(m_pRelx, "relaxation");
    ADD_PROPERTY(m_pWCDF, "recruitment", FEProperty::Optional);

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
	m_pRelx = nullptr;
    m_pWCDF = nullptr;
    
    m_pDmg = nullptr;
    m_pFtg = nullptr;
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
    
    if (!m_pBase->Init()) return false;
    if (!m_pBond->Init()) return false;
    if (!m_pRelx->Init()) return false;
    if (m_pWCDF && !m_pWCDF->Init()) return false;
    
    m_pDmg = dynamic_cast<FEDamageMaterial*>(m_pBase);
    m_pFtg = dynamic_cast<FEReactiveFatigue*>(m_pBase);

    return FEElasticMaterial::Init();
}

//-----------------------------------------------------------------------------
//! Create material point data  arrayfor this material
FEMaterialPointData* FEReactiveViscoelasticMaterial::CreateMaterialPointData()
{
    FEReactiveViscoelasticMaterialPoint* pt = new FEReactiveViscoelasticMaterialPoint();
    // create materal point for strong bond (base) material
    FEMaterialPoint* pbase = new FEMaterialPoint(m_pBase->CreateMaterialPointData());
    pt->AddMaterialPoint(pbase);

    // create materal point for weak bond material
    FEReactiveVEMaterialPoint* pbond = new FEReactiveVEMaterialPoint(m_pBond->CreateMaterialPointData());
    pt->AddMaterialPoint(new FEMaterialPoint(pbond));
    
	return pt;
}

//-----------------------------------------------------------------------------
//! get base material point
FEMaterialPoint* FEReactiveViscoelasticMaterial::GetBaseMaterialPoint(FEMaterialPoint& mp)
{
    // get the reactive viscoelastic point data
    FEReactiveViscoelasticMaterialPoint& rvp = *mp.ExtractData<FEReactiveViscoelasticMaterialPoint>();
    
    // get the elastic material point
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // extract the strong bond material point data
    FEMaterialPoint* sb = rvp.GetPointData(0);
    sb->m_elem = mp.m_elem;
    sb->m_index = mp.m_index;
	sb->m_rt = mp.m_rt;
	sb->m_r0 = mp.m_r0;

    // copy the elastic material point data to the strong bond component
    FEElasticMaterialPoint& epi = *sb->ExtractData<FEElasticMaterialPoint>();
    epi.m_F = ep.m_F;
    epi.m_J = ep.m_J;
    epi.m_v = ep.m_v;
    epi.m_a = ep.m_a;
    epi.m_L = ep.m_L;

    return sb;
}

//-----------------------------------------------------------------------------
//! get bond material point
FEMaterialPoint* FEReactiveViscoelasticMaterial::GetBondMaterialPoint(FEMaterialPoint& mp)
{
    // get the reactive viscoelastic point data
    FEReactiveViscoelasticMaterialPoint& rvp = *mp.ExtractData<FEReactiveViscoelasticMaterialPoint>();
    
    // get the elastic material point data
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // extract the weak bond material point data
    FEMaterialPoint* wb = rvp.GetPointData(1);
    wb->m_elem = mp.m_elem;
    wb->m_index = mp.m_index;
	wb->m_rt = mp.m_rt;
	wb->m_r0 = mp.m_r0;

    // copy the elastic material point data to the weak bond component
    FEElasticMaterialPoint& epi = *wb->ExtractData<FEElasticMaterialPoint>();
    epi.m_F = ep.m_F;
    epi.m_J = ep.m_J;
    epi.m_v = ep.m_v;
    epi.m_a = ep.m_a;
    epi.m_L = ep.m_L;
    
    return wb;
}

//-----------------------------------------------------------------------------
//! detect new generation
bool FEReactiveViscoelasticMaterial::NewGeneration(FEMaterialPoint& mp)
{
    double d;
    double eps = max(m_emin, 10*std::numeric_limits<double>::epsilon());

    // check if the reforming bond mass fraction is above the minimum threshold wmin
    if (ReformingBondMassFraction(mp) < m_wmin) return false;
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
    // get the elastic point data
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // check if the current deformation gradient is different from that of
    // the last generation, in which case store the current state
    // evaluate the relative deformation gradient
    mat3d F = ep.m_F;
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
    double time = CurrentTime();
    double dtv = time - pt.m_v[ig];

    switch (m_btype) {
        case 1:
        {
            if (dtv >= 0)
                w = pt.m_f[ig]*m_pRelx->Relaxation(mp, dtv, D);
        }
            break;
        case 2:
        {
            if (ig == 0) {
                w = m_pRelx->Relaxation(mp, dtv, D);
            }
            else
            {
                double dtu = time - pt.m_v[ig-1];
                w = m_pRelx->Relaxation(mp, dtv, D) - m_pRelx->Relaxation(mp, dtu, D);
            }
        }
            break;
            
        default:
            break;
    }
    
    assert(w >= 0);
    
    return w;
}

//-----------------------------------------------------------------------------
//! evaluate bond mass fraction of reforming generation
double FEReactiveViscoelasticMaterial::ReformingBondMassFraction(FEMaterialPoint& mp)
{
    // get the elastic material point data
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
        
    mat3ds D = ep.RateOfDeformation();
    
    // keep safe copy of deformation gradient
    mat3d F = ep.m_F;
    double J = ep.m_J;
    
    // get current number of generations
    int ng = (int)pt.m_Uv.size();
    
    double f = (!pt.m_wv.empty()) ? pt.m_wv.back() : 1;
    
    for (int ig=0; ig<ng-1; ++ig)
    {
        // evaluate deformation gradient when this generation starts breaking
        ep.m_F = pt.m_Uv[ig];
        ep.m_J = pt.m_Jv[ig];
        // evaluate the breaking bond mass fraction for this generation
        f -= BreakingBondMassFraction(mp, ig, D);
    }
    
    // restore safe copy of deformation gradient
    ep.m_F = F;
    ep.m_J = J;
    
    assert(f >= 0);
    
    // return the bond mass fraction of the reforming generation
    return f;
}

//-----------------------------------------------------------------------------
//! Stress function for strong bonds
mat3ds FEReactiveViscoelasticMaterial::StressStrongBonds(FEMaterialPoint& mp)
{
    mat3ds s = m_pBase->Stress(*GetBaseMaterialPoint(mp));

    return s;
}

//-----------------------------------------------------------------------------
//! Stress function for weak bonds
mat3ds FEReactiveViscoelasticMaterial::StressWeakBonds(FEMaterialPoint& mp)
{
    double dt = CurrentTime();
    if (dt == 0) return mat3ds(0, 0, 0, 0, 0, 0);
    
    FEMaterialPoint& wb = *GetBondMaterialPoint(mp);

    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *wb.ExtractData<FEReactiveVEMaterialPoint>();
    
    // get the elastic point data
    FEElasticMaterialPoint& ep = *wb.ExtractData<FEElasticMaterialPoint>();
    
    // get fiber material point data (if it exists)
    FEFiberMaterialPoint* fp = wb.ExtractData<FEFiberMaterialPoint>();
    
    mat3ds D = ep.RateOfDeformation();
    
    // calculate the base material Cauchy stress
    mat3ds s; s.zero();
    
    // current number of breaking generations
    int ng = (int)pt.m_Uv.size();
    
    // no bonds have broken
    if (ng == 0) {
        s += m_pBond->Stress(wb);
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
            ep.m_F = pt.m_Uv[ig];
            ep.m_J = pt.m_Jv[ig];
            w = BreakingBondMassFraction(wb, ig, D);
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
            sb = m_pBond->Stress(wb);
            // add bond stress to total stress
            s += (ig > 0) ? sb*w/pt.m_Jv[ig-1] : sb*w;
        }
        
        // restore safe copy of deformation gradient
        ep.m_F = F;
        ep.m_J = J;
    }
    
    ep.m_s = s;
    
    return s;
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FEReactiveViscoelasticMaterial::Stress(FEMaterialPoint& mp)
{
    mat3ds s = StressStrongBonds(mp);
    s+= StressWeakBonds(mp)*(1-Damage(mp));
    
	// return the total Cauchy stress
	return s;
}

//-----------------------------------------------------------------------------
//! Material tangent for strong bonds
tens4ds FEReactiveViscoelasticMaterial::TangentStrongBonds(FEMaterialPoint& mp)
{
     // calculate the base material tangent
    return m_pBase->Tangent(*GetBaseMaterialPoint(mp));
}

//-----------------------------------------------------------------------------
//! Material tangent for weak bonds
tens4ds FEReactiveViscoelasticMaterial::TangentWeakBonds(FEMaterialPoint& mp)
{
    FEMaterialPoint& wb = *GetBondMaterialPoint(mp);
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *wb.ExtractData<FEReactiveVEMaterialPoint>();
    
    // get the elastic point data
    FEElasticMaterialPoint& ep = *wb.ExtractData<FEElasticMaterialPoint>();

    // get fiber material point data (if it exists)
    FEFiberMaterialPoint* fp = wb.ExtractData<FEFiberMaterialPoint>();
    
    mat3ds D = ep.RateOfDeformation();
    
    // calculate the base material tangent
    tens4ds c; c.zero();
    
    // current number of breaking generations
    int ng = (int)pt.m_Uv.size();
    
    // no bonds have broken
    if (ng == 0) {
        c += m_pBond->Tangent(wb);
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
            ep.m_F = pt.m_Uv[ig];
            ep.m_J = pt.m_Jv[ig];
            w = BreakingBondMassFraction(wb, ig, D);
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
            cb = m_pBond->Tangent(wb);
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
    c += TangentWeakBonds(mp)*(1-Damage(mp));
    
	// return the total tangent
	return c;
}

//-----------------------------------------------------------------------------
//! strain energy density function in strong bonds
double FEReactiveViscoelasticMaterial::StrongBondSED(FEMaterialPoint& mp)
{
    // calculate the base material strain energy density
    return m_pBase->StrainEnergyDensity(*GetBaseMaterialPoint(mp));
}

//-----------------------------------------------------------------------------
//! strain energy density function in weak bonds
double FEReactiveViscoelasticMaterial::WeakBondSED(FEMaterialPoint& mp)
{
    double dt = CurrentTime();
    if (dt == 0) return 0;
    
    FEMaterialPoint& wb = *GetBondMaterialPoint(mp);
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *wb.ExtractData<FEReactiveVEMaterialPoint>();
    
    // get the elastic point data
    FEElasticMaterialPoint& ep = *wb.ExtractData<FEElasticMaterialPoint>();

    // get fiber material point data (if it exists)
    FEFiberMaterialPoint* fp = wb.ExtractData<FEFiberMaterialPoint>();
    
    // get the viscous point data
    mat3ds D = ep.RateOfDeformation();
    
    double sed = 0;
    
    // current number of breaking generations
    int ng = (int)pt.m_Uv.size();
    
    // no bonds have broken
    if (ng == 0) {
        sed += m_pBond->StrainEnergyDensity(wb);
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
            ep.m_F = pt.m_Uv[ig];
            ep.m_J = pt.m_Jv[ig];
            w = BreakingBondMassFraction(wb, ig, D);
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
            sedb = m_pBond->StrainEnergyDensity(wb);
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
    sed += WeakBondSED(mp)*(1-Damage(mp));
    
    // return the total strain energy density
    return sed;
}

//-----------------------------------------------------------------------------
//! Cull generations that have relaxed below a threshold
void FEReactiveViscoelasticMaterial::CullGenerations(FEMaterialPoint& mp)
{
    // get the elastic material point data
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
    mat3ds D = ep.RateOfDeformation();
    
    // keep safe copy of deformation gradient
    mat3d F = ep.m_F;
    double J = ep.m_J;
    
    int ng = (int)pt.m_v.size();
    m_nmax = max(m_nmax, ng);
    
    // don't cull if we have too few generations
    if (ng < 3) return;

    // don't reduce number of generations to less than max value achieved so far
    if (ng < m_nmax) return;

    // always check oldest generation
    ep.m_F = pt.m_Uv[0];
    ep.m_J = pt.m_Jv[0];
    double w0 = BreakingBondMassFraction(mp, 0, D);
    if (w0 < m_wmin) {
        ep.m_F = pt.m_Uv[1];
        ep.m_J = pt.m_Jv[1];
        double w1 = BreakingBondMassFraction(mp, 1, D);
        pt.m_v[1] = (w0*pt.m_v[0] + w1*pt.m_v[1])/(w0+w1);
        pt.m_Uv[1] = (pt.m_Uv[0]*w0 + pt.m_Uv[1]*w1)/(w0+w1);
        pt.m_Jv[1] = pt.m_Uv[1].det();
        pt.m_f[1] = (w0*pt.m_f[0] + w1*pt.m_f[1])/(w0+w1);
        pt.m_wv[1] = (w0*pt.m_wv[0] + w1*pt.m_wv[1])/(w0+w1);
        pt.m_Uv.pop_front();
        pt.m_Jv.pop_front();
        pt.m_v.pop_front();
        pt.m_f.pop_front();
        pt.m_wv.pop_front();
    }
    
    // restore safe copy of deformation gradient
    ep.m_F = F;
    ep.m_J = J;

    return;
}

//-----------------------------------------------------------------------------
//! Update specialized material points
void FEReactiveViscoelasticMaterial::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp)
{
    FEMaterialPoint& sb = *GetBaseMaterialPoint(mp);
    FEMaterialPoint& wb = *GetBondMaterialPoint(mp);
    
    // start by updating specialized material points of base and bond materials
    m_pBase->UpdateSpecializedMaterialPoints(sb, tp);
    m_pBond->UpdateSpecializedMaterialPoints(wb, tp);
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *wb.ExtractData<FEReactiveVEMaterialPoint>();
    
    // get the elastic point data
    FEElasticMaterialPoint& ep = *wb.ExtractData<FEElasticMaterialPoint>();

    mat3ds Uv = ep.RightStretch();
    double Jv = ep.m_J;

    // if new generation not already created for current time, check if it should
    if (pt.m_v.empty() || (pt.m_v.back() < tp.currentTime)) {
        // check if the current deformation gradient is different from that of
        // the last generation, in which case store the current state
        if (NewGeneration(wb)) {
            pt.m_v.push_back(tp.currentTime);
            pt.m_Uv.push_back(Uv);
            pt.m_Jv.push_back(Jv);
            if (m_pWCDF) {
                pt.m_Et = ScalarStrain(mp);
                if (pt.m_Et > pt.m_Em)
                    pt.m_wv.push_back(m_pWCDF->cdf(mp,pt.m_Et));
                else
                    pt.m_wv.push_back(m_pWCDF->cdf(mp,pt.m_Em));
            }
            else pt.m_wv.push_back(1);
            double f = (!pt.m_v.empty()) ? ReformingBondMassFraction(wb) : 1;
            pt.m_f.push_back(f);
            CullGenerations(wb);
        }
    }
    // otherwise, if we already have a generation for the current time, update the stored values
    else if (pt.m_v.back() == tp.currentTime) {
        pt.m_Uv.back() = Uv;
        pt.m_Jv.back() = Jv;
        if (m_pWCDF) {
            pt.m_Et = ScalarStrain(mp);
            if (pt.m_Et > pt.m_Em)
                pt.m_wv.back() = m_pWCDF->cdf(mp,pt.m_Et);
            else
                pt.m_wv.back() = m_pWCDF->cdf(mp,pt.m_Em);
        }
        pt.m_f.back() = ReformingBondMassFraction(wb);
    }
}

//-----------------------------------------------------------------------------
//! evaluate bond mass fraction of reforming generation
int FEReactiveViscoelasticMaterial::RVEGenerations(FEMaterialPoint& mp)
{
    FEMaterialPoint& wb = *GetBondMaterialPoint(mp);
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *wb.ExtractData<FEReactiveVEMaterialPoint>();
    
    // return the bond mass fraction of the reforming generation
    return (int)pt.m_v.size();
}

//-----------------------------------------------------------------------------
//! evaluate trigger strain
double FEReactiveViscoelasticMaterial::ScalarStrain(FEMaterialPoint& mp)
{
    double d;
    // get the elastic point data
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    switch (m_ttype) {
        case 0:
        {
            // evaluate the Lagrangian strain
            mat3ds E = ep.Strain();
            
            d = E.norm();
        }
            break;
        case 1:
        {
            // distortional strain
            // evaluate spatial Hencky (logarithmic) strain
            mat3ds h = ep.LeftHencky();
            
            // evaluate distortion magnitude (always positive)
            d = (h.dev()).norm();
        }
            break;
        case 2:
        {
            // dilatational strain
            d = fabs(log(ep.m_J));
        }
            break;
            
        default:
            d = 0;
            break;
    }
    
    return d;
}

//-----------------------------------------------------------------------------
double FEReactiveViscoelasticMaterial::Damage(FEMaterialPoint& mp)
{
    double D = 0;
    if (m_pDmg) D = m_pDmg->Damage(*GetBaseMaterialPoint(mp));
    else if (m_pFtg) D = m_pFtg->Damage(*GetBaseMaterialPoint(mp));
    return D;
}
