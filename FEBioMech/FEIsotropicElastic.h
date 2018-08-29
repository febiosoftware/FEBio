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
	virtual mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) override;

	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
    //! calculate the 2nd Piola-Kirchhoff stress at material point
    mat3ds PK2Stress(FEMaterialPoint& pt, const mat3ds E) override;
    
    //! calculate material tangent stiffness at material point
    tens4ds MaterialTangent(FEMaterialPoint& pt, const mat3ds E) override;
    
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
