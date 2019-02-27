#include "stdafx.h"
#include "FERemodelingElasticMaterial.h"
#include "FECore/FECoreKernel.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
FEMaterialPoint* FERemodelingMaterialPoint::Copy()
{
	FERemodelingMaterialPoint* pt = new FERemodelingMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
void FERemodelingMaterialPoint::Init()
{
	// intialize data to zero
	m_sed = m_dsed = 0; 
	m_rhor = m_rhorp = 0;
        
	// don't forget to initialize the base class
    FEMaterialPoint::Init();
}

//-----------------------------------------------------------------------------
void FERemodelingMaterialPoint::Update(const FETimeInfo& timeInfo)
{
	m_rhorp = m_rhor;
        
	// don't forget to initialize the base class
    FEMaterialPoint::Update(timeInfo);
}

//-----------------------------------------------------------------------------
void FERemodelingMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPoint::Serialize(ar);
	ar & m_sed & m_dsed;
	ar & m_rhor & m_rhorp;
}

//=============================================================================
// FERemodelingElasticMaterial
//=============================================================================

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FERemodelingElasticMaterial, FEElasticMaterial)
	ADD_PARAMETER(m_rhormin, "min_density");
	ADD_PARAMETER(m_rhormax, "max_density");

	ADD_PROPERTY(m_pBase, "solid");
	ADD_PROPERTY(m_pSupp, "supply");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERemodelingElasticMaterial::FERemodelingElasticMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_pBase = 0;
	m_pSupp = 0;
}

//-----------------------------------------------------------------------------
//! Strain energy density function
double FERemodelingElasticMaterial::StrainEnergyDensity(FEMaterialPoint& mp)
{
	return (dynamic_cast<FERemodelingInterface*>((FEElasticMaterial*)m_pBase))->StrainEnergy(mp);
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FERemodelingElasticMaterial::Stress(FEMaterialPoint& mp)
{
	double dt = GetFEModel()->GetTime().timeIncrement;

    FERemodelingMaterialPoint& rpt = *(mp.ExtractData<FERemodelingMaterialPoint>());

	// calculate the strain energy density at this material point
	rpt.m_sed = StrainEnergyDensity(mp);

	// calculate the sed derivative with respect to mass density at this material point
    rpt.m_dsed = Tangent_SE_Density(mp);
                
	double rhorhat = m_pSupp->Supply(mp);
	rpt.m_rhor = rhorhat*dt + rpt.m_rhorp;
	if (rpt.m_rhor > m_rhormax) rpt.m_rhor = m_rhormax;
	if (rpt.m_rhor < m_rhormin) rpt.m_rhor = m_rhormin;

	return m_pBase->Stress(mp);
}

//-----------------------------------------------------------------------------
//! Tangent of stress with strain
tens4ds FERemodelingElasticMaterial::Tangent(FEMaterialPoint& mp)
{
	return m_pBase->Tangent(mp);
}

//-----------------------------------------------------------------------------
//! Tangent of strain energy density with mass density
double FERemodelingElasticMaterial::Tangent_SE_Density(FEMaterialPoint& pt)
{
	return (dynamic_cast<FERemodelingInterface*>((FEElasticMaterial*)m_pBase))->Tangent_SE_Density(pt);
}

//-----------------------------------------------------------------------------
//! Tangent of stress with mass density
mat3ds FERemodelingElasticMaterial::Tangent_Stress_Density(FEMaterialPoint& pt)
{
	return (dynamic_cast<FERemodelingInterface*>((FEElasticMaterial*)m_pBase))->Tangent_Stress_Density(pt);
}
