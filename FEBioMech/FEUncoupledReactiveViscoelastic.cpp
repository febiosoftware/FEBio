//
//  FEUncoupledReactiveViscoelastic.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 12/12/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEUncoupledReactiveViscoelastic.h"
#include "FECore/FECoreKernel.h"
#include <limits>

///////////////////////////////////////////////////////////////////////////////
//
// FEUncoupledReactiveViscoelasticMaterial
//
///////////////////////////////////////////////////////////////////////////////

// Material parameters for the FEUncoupledReactiveViscoelastic material
BEGIN_PARAMETER_LIST(FEUncoupledReactiveViscoelasticMaterial, FEUncoupledMaterial)
ADD_PARAMETER(m_wmin, FE_PARAM_DOUBLE, "wmin");
ADD_PARAMETER(m_btype, FE_PARAM_INT, "kinetics");
ADD_PARAMETER(m_ttype, FE_PARAM_INT, "trigger");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor
FEUncoupledReactiveViscoelasticMaterial::FEUncoupledReactiveViscoelasticMaterial(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
    m_pBase = 0;
    m_pBond = 0;
    m_pRelx = 0;
    
    m_wmin = 0;
    m_btype = 0;
    m_ttype = 0;
}

//-----------------------------------------------------------------------------
void FEUncoupledReactiveViscoelasticMaterial::SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp)
{
    FEElasticMaterial::SetLocalCoordinateSystem(el, n, mp);
    FEElasticMaterial* pme = GetBaseMaterial();
    pme->SetLocalCoordinateSystem(el, n, mp);
    FEElasticMaterial* pmb = GetBondMaterial();
    pmb->SetLocalCoordinateSystem(el, n, mp);
}

//-----------------------------------------------------------------------------
//! data initialization
void FEUncoupledReactiveViscoelasticMaterial::Init()
{
    FEUncoupledMaterial::Init();
    if (m_pBase == 0) throw MaterialError("This material needs an elastic base.");
    if (m_pBond == 0) throw MaterialError("This material needs an elastic bond.");
    if (m_pRelx == 0) throw MaterialError("This material needs a bond relaxation.");
    if (!INRANGE(m_wmin, 0.0, 1.0)) throw MaterialError("wmin must be in the range 0 <= wmin <= 1");
    if ((m_btype != 1) && (m_btype != 2)) throw MaterialError("kinetics must be 1 or 2");
    if (!INRANGE(m_ttype, 0, 2)) throw MaterialError("trigger must be 0, 1 or 2");
    
    m_pBase->Init();
    m_pBond->Init();
    m_pRelx->Init();
    
    // set the mixture's bulk modulus based on the base and bond materials
    m_K = m_pBase->m_K + m_pBond->m_K;
}

//-----------------------------------------------------------------------------
//! This material only has one property
int FEUncoupledReactiveViscoelasticMaterial::Properties()
{
    return 3;
}

//-----------------------------------------------------------------------------
FECoreBase* FEUncoupledReactiveViscoelasticMaterial::GetProperty(int i)
{
    if (i == 0) return m_pBase;
    else if (i == 1) return m_pBond;
    else if (i == 2) return m_pRelx;
    assert(false);
    return 0;
}

//-----------------------------------------------------------------------------
//! find a material property index ( returns <0 for error)
int FEUncoupledReactiveViscoelasticMaterial::FindPropertyIndex(const char* szname)
{
    if (strcmp(szname, "elastic") == 0) return 0;
    else if (strcmp(szname, "bond") == 0) return 1;
    else if (strcmp(szname, "relaxation") == 0) return 2;
    else return -1;
}

//-----------------------------------------------------------------------------
//! set a material property (returns false on error)
bool FEUncoupledReactiveViscoelasticMaterial::SetProperty(int i, FECoreBase* pm)
{
    switch(i)
    {
        case 0:
        {
            FEUncoupledMaterial* pme = dynamic_cast<FEUncoupledMaterial*>(pm);
            if (pme) { m_pBase = pme; return true; }
        }
            break;
        case 1:
        {
            FEUncoupledMaterial* pmb = dynamic_cast<FEUncoupledMaterial*>(pm);
            if (pmb) { m_pBond = pmb; return true; }
        }
            break;
        case 2:
        {
            FEBondRelaxation* pmr = dynamic_cast<FEBondRelaxation*>(pm);
            if (pmr) { m_pRelx = pmr; return true; }
        }
            break;
    }
    return false;
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPoint* FEUncoupledReactiveViscoelasticMaterial::CreateMaterialPointData()
{
    return new FEReactiveVEMaterialPoint(m_pBase->CreateMaterialPointData(), this);
}

//-----------------------------------------------------------------------------
//! detect new generation
bool FEUncoupledReactiveViscoelasticMaterial::NewGeneration(FEMaterialPoint& mp)
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
double FEUncoupledReactiveViscoelasticMaterial::BreakingBondMassFraction(FEMaterialPoint& mp, const int ig)
{
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
    // bond mass fraction
    double w = 0;
    
    // current time
    double time = FEMaterialPoint::time;
    
    switch (m_btype) {
        case 1:
        {
            // time when this generation started breaking
            double v = pt.m_v[ig];
            
            if (time >= v)
                w = pt.m_w[ig]*m_pRelx->Relaxation(mp, time - v);
        }
            break;
        case 2:
        {
            double tu, tv;
            if (ig == 0) {
                tv = time - pt.m_v[ig];
                w = m_pRelx->Relaxation(mp, tv);
            }
            else
            {
                tu = time - pt.m_v[ig-1];
                tv = time - pt.m_v[ig];
                w = m_pRelx->Relaxation(mp, tv) - m_pRelx->Relaxation(mp, tu);
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
double FEUncoupledReactiveViscoelasticMaterial::ReformingBondMassFraction(FEMaterialPoint& mp)
{
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
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
        w -= BreakingBondMassFraction(mp, ig);
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
mat3ds FEUncoupledReactiveViscoelasticMaterial::DevStress(FEMaterialPoint& mp)
{
    if (mp.dt == 0) return mat3ds(0,0,0,0,0,0);
    
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
    // calculate the base material Cauchy stress
    mat3ds s = m_pBase->DevStress(mp);
    
    // current number of breaking generations
    int ng = (int)pt.m_Fi.size();
    
    // no bonds have broken
    if (ng == 0) {
        s += m_pBond->DevStress(mp);
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
            w = BreakingBondMassFraction(mp, ig);
            // evaluate bond stress
            sb = m_pBond->DevStress(mp);
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
tens4ds FEUncoupledReactiveViscoelasticMaterial::DevTangent(FEMaterialPoint& mp)
{
    CullGenerations(mp);
    
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
    // calculate the base material tangent
    tens4ds c = m_pBase->DevTangent(mp);
    
    // current number of breaking generations
    int ng = (int)pt.m_Fi.size();
    
    // no bonds have broken
    if (ng == 0) {
        c += m_pBond->DevTangent(mp);
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
            w = BreakingBondMassFraction(mp, ig);
            // evaluate bond tangent
            cb = m_pBond->DevTangent(mp);
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
double FEUncoupledReactiveViscoelasticMaterial::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
    if (mp.dt == 0) return 0;
    
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
    // calculate the base material Cauchy stress
    double sed = m_pBase->DevStrainEnergyDensity(mp);
    
    // current number of breaking generations
    int ng = (int)pt.m_Fi.size();
    
    // no bonds have broken
    if (ng == 0) {
        sed += m_pBond->DevStrainEnergyDensity(mp);
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
            w = BreakingBondMassFraction(mp, ig);
            // evaluate bond stress
            sedb = m_pBond->DevStrainEnergyDensity(mp);
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
void FEUncoupledReactiveViscoelasticMaterial::CullGenerations(FEMaterialPoint& mp)
{
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
    if (pt.m_Fi.empty()) return;
    
    // culling termination flag
    bool done = false;
    
    // always check oldest generation
    while (!done) {
        double w = BreakingBondMassFraction(mp, 0);
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

//-----------------------------------------------------------------------------
//! Get a material parameter
FEParam* FEUncoupledReactiveViscoelasticMaterial::GetParameter(const ParamString& s)
{
    // see if this is a composite parameter
    if (s.count() == 1) return FEMaterial::GetParameter(s);
    
    // else find the component's parameter
    if      (s == "elastic"       ) return m_pBase->GetParameter(s.next());
    else if (s == "bond"          ) return m_pBond->GetParameter(s.next());
    else if (s == "relaxation"    ) return m_pRelx->GetParameter(s.next());
    else return 0;
}

//-----------------------------------------------------------------------------
//! Save data to dump file

void FEUncoupledReactiveViscoelasticMaterial::Serialize(DumpFile& ar)
{
    // serialize material parameters
    FEMaterial::Serialize(ar);
    
    // serialize sub-materials
    if (ar.IsSaving())
    {
        ar << m_pBase->GetTypeStr();
        m_pBase->Serialize(ar);
        
        ar << m_pBond->GetTypeStr();
        m_pBond->Serialize(ar);
        
        ar << m_pRelx->GetTypeStr();
        m_pRelx->Serialize(ar);
    }
    else
    {
        char sz[256] = {0};
        
        ar >> sz;
        m_pBase = dynamic_cast<FEUncoupledMaterial*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
        assert(m_pBase);
        m_pBase->Serialize(ar);
        m_pBase->Init();
        
        ar >> sz;
        m_pBond = dynamic_cast<FEUncoupledMaterial*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
        assert(m_pBond);
        m_pBond->Serialize(ar);
        m_pBond->Init();
        
        ar >> sz;
        m_pRelx = dynamic_cast<FEBondRelaxation*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
        assert(m_pRelx);
        m_pRelx->Serialize(ar);
        m_pRelx->Init();
    }
}
