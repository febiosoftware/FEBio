//
//  FEDamageMaterial.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 9/18/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEDamageMaterial.h"
#include "FEDamageCriterion.h"
#include "FEDamageCDF.h"
#include "FEUncoupledMaterial.h"
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageMaterial::FEDamageMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
}

//-----------------------------------------------------------------------------
//! Initialization.
void FEDamageMaterial::Init()
{
    FEUncoupledMaterial* m_pMat = dynamic_cast<FEUncoupledMaterial*>(m_pBase);
    if (m_pMat != NULL)
        throw MaterialError("Elastic material should not be of type uncoupled");
    
    // set parent materials
    m_pBase->SetParent(this);
    m_pDamg->SetParent(this);
    m_pCrit->SetParent(this);
}

//-----------------------------------------------------------------------------
//! calculate stress at material point
mat3ds FEDamageMaterial::Stress(FEMaterialPoint& pt)
{
    // get the damage material point data
	FEDamageMaterialPoint& pd = *pt.ExtractData<FEDamageMaterialPoint>();
    
    // get the damage
    double d = pd.m_D;

    // evaluate the stress
    mat3ds s = m_pBase->Stress(pt);
    
    // return damaged stress
    return s*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate tangent stiffness at material point
tens4ds FEDamageMaterial::Tangent(FEMaterialPoint& pt)
{
    // get the damage material point data
	FEDamageMaterialPoint& pd = *pt.ExtractData<FEDamageMaterialPoint>();
    
    // evaluate the trial value of the damage criterion
    // this must be done before evaluating the damage
    pd.m_Etrial = m_pCrit->DamageCriterion(pt);
    
    // evaluate the damage
    double d = Damage(pt);
    
    // evaluate the tangent
    tens4ds c = m_pBase->Tangent(pt);
    
    // return damaged tangent
    return c*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FEDamageMaterial::StrainEnergyDensity(FEMaterialPoint& pt)
{
    // get the damage material point data
	FEDamageMaterialPoint& pd = *pt.ExtractData<FEDamageMaterialPoint>();
    
    // get the damage
    double d = pd.m_D;
    
    // evaluate the strain energy density
    double sed = m_pBase->StrainEnergyDensity(pt);

    // return damaged sed
    return sed*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate damage at material point
double FEDamageMaterial::Damage(FEMaterialPoint& pt)
{
    // get the damage material point data
    FEDamageMaterialPoint& pd = *pt.ExtractData<FEDamageMaterialPoint>();

    // evaluate the trial value of the damage criterion
    // this must be done before evaluating the damage
    pd.m_Etrial = m_pCrit->DamageCriterion(pt);
    
    // evaluate and set the damage
    double d = m_pDamg->Damage(pt);
    pd.m_D = d;
    
    return d;
}

//-----------------------------------------------------------------------------
int FEDamageMaterial::Properties()
{
	return 3;
}

//-----------------------------------------------------------------------------
//! get a specific material property
FECoreBase* FEDamageMaterial::GetProperty(int i)
{
	switch(i)
	{
        case 0: return m_pBase; break;
        case 1: return m_pDamg; break;
        case 2: return m_pCrit; break;
	}
	return 0;
}

//-----------------------------------------------------------------------------
//! find a material property index ( returns <0 for error)
int FEDamageMaterial::FindPropertyIndex(const char* szname)
{
	if      (strcmp(szname, "elastic"   ) == 0) return 0;
	else if (strcmp(szname, "damage"    ) == 0) return 1;
	else if (strcmp(szname, "criterion" ) == 0) return 2;
	else return -1;
}

//-----------------------------------------------------------------------------
//! set a material property (returns false on error)
bool FEDamageMaterial::SetProperty(int i, FECoreBase* pm)
{
	switch(i)
	{
        case 0:
		{
			m_pBase = dynamic_cast<FEElasticMaterial*>(pm);
			if ((m_pBase == 0) || (m_pBase->IsRigid())) return false;
			return true;
		}
            break;
        case 1:
		{
			m_pDamg = dynamic_cast<FEDamageCDF*>(pm);
			if (m_pDamg == 0) return false;
			return true;
		}
            break;
        case 2:
		{
			m_pCrit = dynamic_cast<FEDamageCriterion*>(pm);
			if (m_pCrit == 0) return false;
			return true;
		}
            break;
	}
	return false;
}

//-----------------------------------------------------------------------------
FEParam* FEDamageMaterial::GetParameter(const ParamString& s)
{
	if (s.count() == 1) return FEMaterial::GetParameter(s);
    
	if      (s == "elastic"     ) return m_pBase->GetParameter(s.next());
	else if (s == "damage"      ) return m_pDamg->GetParameter(s.next());
	else if (s == "criterion"   ) return m_pCrit->GetParameter(s.next());
	return 0;
}

//-----------------------------------------------------------------------------
void FEDamageMaterial::SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp)
{
	FEElasticMaterial::SetLocalCoordinateSystem(el, n, mp);
	m_pBase->SetLocalCoordinateSystem(el, n, mp);
}

//-----------------------------------------------------------------------------
void FEDamageMaterial::Serialize(DumpFile& ar)
{
	// serialize material parameters
	FEElasticMaterial::Serialize(ar);

	// serialize sub-materials
	if (ar.IsSaving())
	{
		ar << m_pBase->GetTypeStr();
		m_pBase->Serialize(ar);

		ar << m_pDamg->GetTypeStr();
		m_pDamg->Serialize(ar);

		ar << m_pCrit->GetTypeStr();
		m_pCrit->Serialize(ar);
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
		m_pDamg = dynamic_cast<FEDamageCDF*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
		assert(m_pDamg);
		m_pDamg->Serialize(ar);
		m_pDamg->Init();

		ar >> sz;
		m_pCrit = dynamic_cast<FEDamageCriterion*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
		assert(m_pCrit);
		m_pCrit->Serialize(ar);
		m_pCrit->Init();
	}
}
