//
//  FEDamageMaterialUC.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 9/19/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEDamageMaterialUC.h"
#include "FEDamageCriterionUC.h"
#include "FEDamageCDF.h"
#include "FEUncoupledMaterial.h"
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageMaterialUC::FEDamageMaterialUC(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
	// set material properties
	AddProperty(&m_pBase, "elastic"  );
	AddProperty(&m_pDamg, "damage"   );
	AddProperty(&m_pCrit, "criterion");
}

//-----------------------------------------------------------------------------
//! Initialization.
void FEDamageMaterialUC::Init()
{
    FEUncoupledMaterial::Init();
    
    // set bulk modulus to that of base elastic material
    m_K = m_pBase->m_K;
}

//-----------------------------------------------------------------------------
//! calculate stress at material point
mat3ds FEDamageMaterialUC::DevStress(FEMaterialPoint& pt)
{
    // get the damage material point data
	FEDamageMaterialPoint& pd = *pt.ExtractData<FEDamageMaterialPoint>();
    
    // evaluate the trial value of the damage criterion
    // this must be done before evaluating the damage
    pd.m_Etrial = m_pCrit->DamageCriterion(pt);
    
    // evaluate the damage
    double d = m_pDamg->Damage(pt);
    
    // evaluate the stress
    mat3ds s = m_pBase->DevStress(pt);
    
    // return the damaged stress
    return s*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate tangent stiffness at material point
tens4ds FEDamageMaterialUC::DevTangent(FEMaterialPoint& pt)
{
    // get the damage material point data
	FEDamageMaterialPoint& pd = *pt.ExtractData<FEDamageMaterialPoint>();
    
    // evaluate the trial value of the damage criterion
    // this must be done before evaluating the damage
    pd.m_Etrial = m_pCrit->DamageCriterion(pt);
    
    // evaluate the damage
    double d = m_pDamg->Damage(pt);
    
    // evaluate the tangent
    tens4ds c = m_pBase->DevTangent(pt);
    
    // return the damaged tangent
    return c*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FEDamageMaterialUC::DevStrainEnergyDensity(FEMaterialPoint& pt)
{
    // get the damage material point data
	FEDamageMaterialPoint& pd = *pt.ExtractData<FEDamageMaterialPoint>();
    
    // evaluate the trial value of the damage criterion
    // this must be done before evaluating the damage
    pd.m_Etrial = m_pCrit->DamageCriterion(pt);

    // evaluate the damage
    double d = m_pDamg->Damage(pt);
    
    // evaluate the strain energy density
    double sed = m_pBase->DevStrainEnergyDensity(pt);
    
    // return the damaged sed
    return sed*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate damage at material point
double FEDamageMaterialUC::Damage(FEMaterialPoint& pt)
{
    // get the damage material point data
    FEDamageMaterialPoint& pd = *pt.ExtractData<FEDamageMaterialPoint>();
    
    // evaluate the trial value of the damage criterion
    // this must be done before evaluating the damage
    pd.m_Etrial = m_pCrit->DamageCriterion(pt);
    
    // evaluate the damage
    double d = m_pDamg->Damage(pt);
    
    return d;
}

//-----------------------------------------------------------------------------
void FEDamageMaterialUC::SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp)
{
	FEUncoupledMaterial::SetLocalCoordinateSystem(el, n, mp);
	m_pBase->SetLocalCoordinateSystem(el, n, mp);
}
