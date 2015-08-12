#pragma once
#include "FERemodelingElasticMaterial.h"

//-----------------------------------------------------------------------------
// This class implements a material that has a constant solute supply

class FEHuiskesSupply :	public FESolidSupply
{
public:
	//! constructor
	FEHuiskesSupply(FEModel* pfem);
	
	//! solid supply
	double Supply(FEMaterialPoint& pt);
	
	//! tangent of solute supply with respect to strain
	mat3ds Tangent_Supply_Strain(FEMaterialPoint& mp);
	
	//! tangent of solute supply with respect to referential density
	double Tangent_Supply_Density(FEMaterialPoint& mp);
	
public:
	double	m_B;			//!< mass supply coefficient
	double	m_k;			//!< specific strain energy at homeostasis
	
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
