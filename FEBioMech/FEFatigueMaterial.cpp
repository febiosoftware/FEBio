//
//  FEFatigueMaterial.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 9/29/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#include "stdafx.h"
#include "FEFatigueMaterial.h"
#include "FEDamageCriterion.h"
#include "FEDamageCDF.h"
#include "FEUncoupledMaterial.h"
#include "FECore/FECoreKernel.h"
#include <FECore/DumpStream.h>

////////////////////// FATIGUE MATERIAL POINT /////////////////////////////////
//-----------------------------------------------------------------------------
FEMaterialPoint* FEFatigueMaterialPoint::Copy()
{
    FEFatigueMaterialPoint* pt = new FEFatigueMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FEFatigueMaterialPoint::Init()
{
    FEMaterialPoint::Init();
    
    // intialize data
    m_D = 0;
    
    // initialize intact bond fraction to 1
    m_wi = m_wip = m_wit = 1.0;
    m_awi = m_awip = m_awit = 0;
    
    // initialize fatigued bond fraction to 0
    m_wf = m_wfp = m_wft = 0;
    m_awf = m_awfp = m_awft = 0;

    // initialize damage criterion
    m_Xd = m_Xdp = m_Xdt = 0;
    m_aXd = m_aXdp = m_aXdt = 0;
    
    // initialize fatigue criterion
    m_Xf = m_Xfp = m_Xft = 0;
    m_aXf = m_aXfp = m_aXft = 0;
}

//-----------------------------------------------------------------------------
void FEFatigueMaterialPoint::Update(const FETimeInfo& timeInfo)
{
    FEMaterialPoint::Update(timeInfo);
    
    // update damage response
    if (m_Xdt > m_Xdp) {
        // intact bonds
        m_wip = m_wit;
        m_awip = m_awit;
        // fatigued bonds
        m_wfp = m_wft;
        m_awfp = m_awft;
        // damage criterion
        m_Xdp = m_Xdt;
        m_aXdp = m_aXdt;
    }
    
    // update fatigue response
    m_wfp = m_wft;
    m_awfp = m_awft;
    m_Xfp = m_Xft;
    m_aXfp = m_aXft;
}

//-----------------------------------------------------------------------------
void FEFatigueMaterialPoint::Serialize(DumpStream& ar)
{
    if (ar.IsSaving())
    {
        ar << m_D;
        ar << m_wi << m_wip << m_wit << m_awi << m_awip << m_awit;
        ar << m_wf << m_wfp << m_wft << m_awf << m_awfp << m_awft;
        ar << m_Xd << m_Xdp << m_Xdt << m_aXd << m_aXdp << m_aXdt;
        ar << m_Xf << m_Xfp << m_Xft << m_aXf << m_aXfp << m_aXft;
    }
    else
    {
        ar >> m_D;
        ar >> m_wi >> m_wip >> m_wit >> m_awi >> m_awip >> m_awit;
        ar >> m_wf >> m_wfp >> m_wft >> m_awf >> m_awfp >> m_awft;
        ar >> m_Xd >> m_Xdp >> m_Xdt >> m_aXd >> m_aXdp >> m_aXdt;
        ar >> m_Xf >> m_Xfp >> m_Xft >> m_aXf >> m_aXfp >> m_aXft;
    }
    FEMaterialPoint::Serialize(ar);
}

//////////////////////////// FATIGUE MATERIAL /////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEFatigueMaterial, FEMaterial)
	ADD_PARAMETER(m_k0   , FE_RANGE_GREATER_OR_EQUAL(0.0), "k0"  );
	ADD_PARAMETER(m_beta , FE_RANGE_GREATER_OR_EQUAL(0.0), "beta");

	// set material properties
	ADD_PROPERTY(m_pBase, "elastic");
	ADD_PROPERTY(m_pIdmg, "intact_damage");
	ADD_PROPERTY(m_pFdmg, "fatigue_damage");
	ADD_PROPERTY(m_pDcrt, "damage_criterion");
	ADD_PROPERTY(m_pFcrt, "fatigue_criterion");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEFatigueMaterial::FEFatigueMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_pBase = 0;
	m_pIdmg = 0;
	m_pFdmg = 0;
	m_pDcrt = 0;
	m_pFcrt = 0;
}

//-----------------------------------------------------------------------------
//! Initialization.
bool FEFatigueMaterial::Init()
{
    FEUncoupledMaterial* m_pMat = dynamic_cast<FEUncoupledMaterial*>((FEElasticMaterial*)m_pBase);
    if (m_pMat != nullptr)
        return fecore_error("Elastic material should not be of type uncoupled");
    
    return FEElasticMaterial::Init();
}

//-----------------------------------------------------------------------------
//! calculate stress at material point
mat3ds FEFatigueMaterial::Stress(FEMaterialPoint& pt)
{
    // evaluate the damage
    double d = Damage(pt);
    
    // evaluate the stress
    mat3ds s = m_pBase->Stress(pt);
    
    // return damaged stress
    return s*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate tangent stiffness at material point
tens4ds FEFatigueMaterial::Tangent(FEMaterialPoint& pt)
{
    // evaluate the damage
    double d = Damage(pt);
    
    // evaluate the tangent
    tens4ds c = m_pBase->Tangent(pt);
    
    // return damaged tangent
    return c*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FEFatigueMaterial::StrainEnergyDensity(FEMaterialPoint& pt)
{
    // evaluate the damage
    double d = Damage(pt);
    
    // evaluate the strain energy density
    double sed = m_pBase->StrainEnergyDensity(pt);
    
    // return damaged sed
    return sed*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate damage at material point
double FEFatigueMaterial::Damage(FEMaterialPoint& pt)
{
    // get the damage material point data
    FEFatigueMaterialPoint& pd = *pt.ExtractData<FEFatigueMaterialPoint>();
    
    pd.m_D = 1 - pd.m_wi - pd.m_wf;
    
    return pd.m_D;
}

//-----------------------------------------------------------------------------
// update fatigue material point at each iteration
void FEFatigueMaterial::UpdateSpecializedMaterialPoints(FEMaterialPoint& pt, const FETimeInfo& tp)
{
    // get the fatigue material point data
    FEFatigueMaterialPoint& pd = *pt.ExtractData<FEFatigueMaterialPoint>();
    
    // evaluate intermediate time point variables
    pd.m_wi = (1-tp.alphaf)*pd.m_wip + tp.alphaf*pd.m_wit;
    pd.m_awi = (1-tp.alpham)*pd.m_awip + tp.alpham*pd.m_awit;
    pd.m_wf = (1-tp.alphaf)*pd.m_wfp + tp.alphaf*pd.m_wft;
    pd.m_awf = (1-tp.alpham)*pd.m_awfp + tp.alpham*pd.m_awft;

    // get damage and fatigue criteria at intermediate time point
    pd.m_Xd = m_pDcrt->DamageCriterion(pt);
    pd.m_Xf = m_pFcrt->DamageCriterion(pt);
    pd.m_Xdt = pd.m_Xdp + (pd.m_Xd - pd.m_Xdp)/tp.alphaf;
    pd.m_Xft = pd.m_Xfp + (pd.m_Xf - pd.m_Xfp)/tp.alphaf;
    pd.m_aXdt = (pd.m_Xd - pd.m_Xdp)/(tp.alphaf*tp.gamma*tp.timeIncrement) - (1.0/tp.gamma - 1.0)*pd.m_aXdp;
    pd.m_aXft = (pd.m_Xf - pd.m_Xfp)/(tp.alphaf*tp.gamma*tp.timeIncrement) - (1.0/tp.gamma - 1.0)*pd.m_aXfp;
    pd.m_aXd = (1-tp.alpham)*pd.m_aXdp + tp.alpham*pd.m_aXdt;
    pd.m_aXf = (1-tp.alpham)*pd.m_aXfp + tp.alpham*pd.m_aXft;

    // evaluate mass supplies from damage
    double fi = 0, ff = 0;
    if (pd.m_Xd > pd.m_Xdp) {
        fi = -m_pIdmg->pdf(pd.m_Xd)/(1-m_pIdmg->cdf(pd.m_Xd))*pd.m_aXd;
        ff = -m_pFdmg->pdf(pd.m_Xd)/(1-m_pFdmg->cdf(pd.m_Xd))*pd.m_aXd;
    }
    
    // evaluate mass supplies from fatigue
    double f = m_k0*pow(fabs(pd.m_Xf), m_beta);
    // deduct from intact bond fraction
    fi -= f;
    // add to fatigued bond fraction
    ff += f;
    
    // update bond mass fractions based on mass supplies
    double dwi = (fi*pd.m_wi - pd.m_awi)/(tp.alpham/tp.gamma/tp.timeIncrement - tp.alphaf*fi);
    double dwf = (ff*pd.m_wf - pd.m_awf)/(tp.alpham/tp.gamma/tp.timeIncrement - tp.alphaf*ff);
    pd.m_wit += dwi;
    pd.m_wft += dwf;
    pd.m_wi = (1-tp.alphaf)*pd.m_wip + tp.alphaf*pd.m_wit;
    pd.m_wf = (1-tp.alphaf)*pd.m_wfp + tp.alphaf*pd.m_wft;

    // update time derivatives
    pd.m_awit = (pd.m_wit - pd.m_wip)/(tp.gamma*tp.timeIncrement) + (1.0-1.0/tp.gamma)*pd.m_awip;
    pd.m_awft = (pd.m_wft - pd.m_wfp)/(tp.gamma*tp.timeIncrement) + (1.0-1.0/tp.gamma)*pd.m_awfp;
    pd.m_awi = (1-tp.alpham)*pd.m_awip + tp.alpham*pd.m_awit;
    pd.m_awf = (1-tp.alpham)*pd.m_awfp + tp.alpham*pd.m_awft;
}
