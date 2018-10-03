#include "FEPermHolmesMow.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEPermHolmesMow, FEHydraulicPermeability)
	ADD_PARAMETER(m_perm , FE_RANGE_GREATER_OR_EQUAL(0.0), "perm" );
	ADD_PARAMETER(m_M    , FE_RANGE_GREATER_OR_EQUAL(0.0), "M"    );
	ADD_PARAMETER(m_alpha, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEPermHolmesMow::FEPermHolmesMow(FEModel* pfem) : FEHydraulicPermeability(pfem)
{
	m_perm = 1;
	m_M = m_alpha = 0;
}

//-----------------------------------------------------------------------------
//! Permeability tensor.
mat3ds FEPermHolmesMow::Permeability(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
	// relative volume
	double J = et.m_J;
	// referential solid volume fraction
	double phi0 = pt.m_phi0;
	
	// --- strain-dependent isotropic permeability ---
	
	return mat3dd(m_perm*pow((J-phi0)/(1.0-phi0),m_alpha)*exp(m_M*(J*J-1.0)/2.0));
}

//-----------------------------------------------------------------------------
//! Tangent of permeability
tens4ds FEPermHolmesMow::Tangent_Permeability_Strain(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
	// relative volume
	double J = et.m_J;
	// referential solid volume fraction
	double phi0 = pt.m_phi0;
	
	mat3dd I(1);	// Identity
	
	double k0 = m_perm*pow((J-phi0)/(1.0-phi0),m_alpha)*exp(m_M*(J*J-1.0)/2.0);
	double K0prime = (J*J*m_M+(J*(m_alpha+1)-phi0)/(J-phi0))*k0;
	mat3ds k0hat = I*K0prime;
	
	return dyad1s(I,k0hat)/2.0-dyad4s(I)*2*k0;
}
