#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
class FEIsotropicElastic : public FEElasticMaterial
{
public:
	FEIsotropicElastic(FEModel* pfem) : FEElasticMaterial(pfem) {}

public:
	double	m_E;	//!< Young's modulus
	double	m_v;	//!< Poisson's ratio

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt);
    
    //! calculate the 2nd Piola-Kirchhoff stress at material point
    mat3ds PK2Stress(FEMaterialPoint& pt, const mat3ds E);
    
    //! calculate material tangent stiffness at material point
    tens4ds MaterialTangent(FEMaterialPoint& pt, const mat3ds E);
    
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
