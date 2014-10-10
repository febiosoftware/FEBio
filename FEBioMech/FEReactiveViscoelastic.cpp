//
//  FEReactiveViscoelastic.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/25/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEReactiveViscoelastic.h"
#include "FECore/FECoreKernel.h"

///////////////////////////////////////////////////////////////////////////////
//
// FEReactiveVEMaterialPoint
//
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//! Create a shallow copy of the material point data
FEMaterialPoint* FEReactiveVEMaterialPoint::Copy()
{
	FEReactiveVEMaterialPoint* pt = new FEReactiveVEMaterialPoint(*this);
	if (m_pt) pt->m_pt = m_pt->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
//! Initializes material point data.
void FEReactiveVEMaterialPoint::Init(bool bflag)
{
	FEElasticMaterialPoint& pt = *m_pt->ExtractData<FEElasticMaterialPoint>();
	if (bflag)
	{
		// initialize data to zero
		m_Fi.clear(); m_Fi.push_back(mat3dd(1.0));
        m_Ji.clear(); m_Ji.push_back(1.0);
        m_tgen.clear(); m_tgen.push_back(FEMaterialPoint::time);
	}
	else
	{
        // check if the current deformation gradient is different from that of
        // the last generation, in which case store the current state
        double eps = 1e-9;
        int lg = (int)m_Fi.size() - 1;
        mat3d D = pt.m_F*m_Fi[lg] - mat3dd(1);
        if (D.norm() > eps) {
            m_Fi.push_back(pt.m_F.inverse());
            m_Ji.push_back(1./pt.m_J);
            m_tgen.push_back(FEMaterialPoint::time);
        }
	}
    
	// don't forget to initialize the nested data
	if (m_pt) m_pt->Init(bflag);
}

//-----------------------------------------------------------------------------
//! Serialize data to the archive
void FEReactiveVEMaterialPoint::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (m_pt) m_pt->ShallowCopy(dmp, bsave);
    
	if (bsave)
	{
        int n = (int)m_Fi.size();
		dmp << n;
		for (int i=0; i<n; ++i) dmp << m_Fi[i] << m_Ji[i] << m_tgen[i];
	}
	else
	{
        int n;
		dmp >> n;
		for (int i=0; i<n; ++i) dmp >> m_Fi[i] >> m_Ji[i] >> m_tgen[i];
	}
}

//-----------------------------------------------------------------------------
//! Serialize data to the archive
void FEReactiveVEMaterialPoint::Serialize(DumpFile& ar)
{
	if (m_pt) m_pt->Serialize(ar);
    
	if (ar.IsSaving())
	{
        int n = (int)m_Fi.size();
		ar << n;
		for (int i=0; i<n; ++i) ar << m_Fi[i] << m_Ji[i] << m_tgen[i];
	}
	else
	{
        int n;
		ar >> n;
		for (int i=0; i<n; ++i) ar >> m_Fi[i] >> m_Ji[i] >> m_tgen[i];
	}
}

///////////////////////////////////////////////////////////////////////////////
//
// FEReactiveViscoelasticMaterial
//
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//! constructor
FEReactiveViscoelasticMaterial::FEReactiveViscoelasticMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_pBase = 0;
    m_pBond = 0;
    m_pRelx = 0;
}

//-----------------------------------------------------------------------------
void FEReactiveViscoelasticMaterial::SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp)
{
	FEElasticMaterial::SetLocalCoordinateSystem(el, n, mp);
	FEElasticMaterial* pme = GetBaseMaterial();
	pme->SetLocalCoordinateSystem(el, n, mp);
	FEElasticMaterial* pmb = GetBondMaterial();
	pmb->SetLocalCoordinateSystem(el, n, mp);
}

//-----------------------------------------------------------------------------
//! data initialization
void FEReactiveViscoelasticMaterial::Init()
{
	FEElasticMaterial::Init();
	if (m_pBase == 0) throw MaterialError("This material needs an elastic base.");
	if (m_pBond == 0) throw MaterialError("This material needs an elastic bond.");
	if (m_pRelx == 0) throw MaterialError("This material needs a bond relaxation.");
	m_pBase->Init();
	m_pBond->Init();
    m_pRelx->Init();
}

//-----------------------------------------------------------------------------
//! This material only has one property
int FEReactiveViscoelasticMaterial::Properties()
{
	return 3;
}

//-----------------------------------------------------------------------------
FECoreBase* FEReactiveViscoelasticMaterial::GetProperty(int i)
{
	if (i == 0) return m_pBase;
	else if (i == 1) return m_pBond;
	else if (i == 2) return m_pRelx;
	assert(false);
	return 0;
}

//-----------------------------------------------------------------------------
//! find a material property index ( returns <0 for error)
int FEReactiveViscoelasticMaterial::FindPropertyIndex(const char* szname)
{
	if (strcmp(szname, "elastic") == 0) return 0;
    else if (strcmp(szname, "bond") == 0) return 1;
    else if (strcmp(szname, "relaxation") == 0) return 2;
    else return -1;
}

//-----------------------------------------------------------------------------
//! set a material property (returns false on error)
bool FEReactiveViscoelasticMaterial::SetProperty(int i, FECoreBase* pm)
{
	switch(i)
	{
        case 0:
		{
			FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(pm);
			if (pme) { m_pBase = pme; return true; }
		}
            break;
        case 1:
		{
			FEElasticMaterial* pmb = dynamic_cast<FEElasticMaterial*>(pm);
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
FEMaterialPoint* FEReactiveViscoelasticMaterial::CreateMaterialPointData()
{
	return new FEReactiveVEMaterialPoint(m_pBase->CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FEReactiveViscoelasticMaterial::Stress(FEMaterialPoint& mp)
{
    if (mp.dt == 0) return mat3ds(0,0,0,0,0,0);
    
    int ig;
    
	// get the elastic part
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// get the reactive viscoelastic point data
	FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
	// calculate the base material Cauchy stress
	mat3ds s = m_pBase->Stress(mp);
    
    // keep safe copy of deformation gradient
    mat3d F = ep.m_F;
    double J = ep.m_J;

    // current time
    double time = FEMaterialPoint::time;
    
    // current number of generations
    int ng = (int)pt.m_Fi.size();
    
    double w;
    mat3ds sb;
    double ta, tb;
    
    // calculate the bond stresses for past generations
    for (ig=0; ig<ng; ++ig) {
        // evaluate relative deformation gradient for this generation
        ep.m_F = F*pt.m_Fi[ig];
        ep.m_J = J*pt.m_Ji[ig];
        // evaluate bond mass fraction for this generation
        if (ig == 0) {
            ta = time - pt.m_tgen[ig];
            w = m_pRelx->Relaxation(mp, ta);
        }
        else if (ig == ng-1) {
            ta = time - pt.m_tgen[ig];
            w = 1 - m_pRelx->Relaxation(mp, ta);
        }
        else
        {
            ta = time - pt.m_tgen[ig];
            tb = time - pt.m_tgen[ig+1];
            w = m_pRelx->Relaxation(mp, tb) - m_pRelx->Relaxation(mp, ta);
        }
        // evaluate bond stress
        sb = m_pBond->Stress(mp);
        // add bond stress to total stress
        s += sb*w;
    }
    
    // restore safe copy of deformation gradient
    ep.m_F = F;
    ep.m_J = J;
    
	// return the total Cauchy stress
	return s;
}

//-----------------------------------------------------------------------------
//! Material tangent
tens4ds FEReactiveViscoelasticMaterial::Tangent(FEMaterialPoint& mp)
{
    int ig;
    
	// get the elastic part
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// get the reactive viscoelastic point data
	FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
	// calculate the base material tangent
	tens4ds c = m_pBase->Tangent(mp);
    
    // keep safe copy of deformation gradient
    mat3d F = ep.m_F;
    double J = ep.m_J;
    
    // current time
    double time = FEMaterialPoint::time;
    
    // current number of generations
    int ng = (int)pt.m_Fi.size();
    
    double w;
    tens4ds cb;
    double ta, tb;
    
    // calculate the bond tangents for past generations
    for (ig=0; ig<ng; ++ig) {
        // evaluate relative deformation gradient for this generation
        ep.m_F = F*pt.m_Fi[ig];
        ep.m_J = J*pt.m_Ji[ig];
        // evaluate bond mass fraction for this generation
        if (ig == 0) {
            ta = time - pt.m_tgen[ig];
            w = m_pRelx->Relaxation(mp, ta);
        }
        else if (ig == ng-1) {
            ta = time - pt.m_tgen[ig];
            w = 1 - m_pRelx->Relaxation(mp, ta);
        }
        else
        {
            ta = time - pt.m_tgen[ig];
            tb = time - pt.m_tgen[ig+1];
            w = m_pRelx->Relaxation(mp, tb) - m_pRelx->Relaxation(mp, ta);
        }
        // evaluate bond tangent
        cb = m_pBond->Tangent(mp);
        // add bond tangent to total tangent
        c += cb*w;
    }
    
    // restore safe copy of deformation gradient
    ep.m_F = F;
    ep.m_J = J;
    
	// return the total tangent
	return c;
}

//-----------------------------------------------------------------------------
//! strain energy density function
double FEReactiveViscoelasticMaterial::StrainEnergyDensity(FEMaterialPoint& mp)
{
    if (mp.dt == 0) return 0;
    
    int ig;
    
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the reactive viscoelastic point data
    FEReactiveVEMaterialPoint& pt = *mp.ExtractData<FEReactiveVEMaterialPoint>();
    
    // calculate the base material Cauchy stress
    double sed = m_pBase->StrainEnergyDensity(mp);
    
    // keep safe copy of deformation gradient
    mat3d F = ep.m_F;
    double J = ep.m_J;
    
    // current time
    double time = FEMaterialPoint::time;
    
    // current number of generations
    int ng = (int)pt.m_Fi.size();
    
    double w;
    double sedb;
    double ta, tb;
    
    // calculate the bond stresses for past generations
    for (ig=0; ig<ng; ++ig) {
        // evaluate relative deformation gradient for this generation
        ep.m_F = F*pt.m_Fi[ig];
        ep.m_J = J*pt.m_Ji[ig];
        // evaluate bond mass fraction for this generation
        if (ig == 0) {
            ta = time - pt.m_tgen[ig];
            w = m_pRelx->Relaxation(mp, ta);
        }
        else if (ig == ng-1) {
            ta = time - pt.m_tgen[ig];
            w = 1 - m_pRelx->Relaxation(mp, ta);
        }
        else
        {
            ta = time - pt.m_tgen[ig];
            tb = time - pt.m_tgen[ig+1];
            w = m_pRelx->Relaxation(mp, tb) - m_pRelx->Relaxation(mp, ta);
        }
        // evaluate bond stress
        sedb = m_pBond->StrainEnergyDensity(mp);
        // add bond stress to total stress
        sed += sedb*w;
    }
    
    // restore safe copy of deformation gradient
    ep.m_F = F;
    ep.m_J = J;
    
    // return the total Cauchy stress
    return sed;
}

//-----------------------------------------------------------------------------
//! Get a material parameter
FEParam* FEReactiveViscoelasticMaterial::GetParameter(const ParamString& s)
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

void FEReactiveViscoelasticMaterial::Serialize(DumpFile& ar)
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
		m_pBase = dynamic_cast<FEElasticMaterial*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
		assert(m_pBase);
		m_pBase->Serialize(ar);
		m_pBase->Init();
        
		ar >> sz;
		m_pBond = dynamic_cast<FEElasticMaterial*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
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
