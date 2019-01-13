#pragma once
#include "FEElasticMaterial.h"
#include <FECore/FEModelParam.h>

//-----------------------------------------------------------------------------
//! Neo Hookean material

//! Implementation of a neo-Hookean hyperelastic material.
class FEBIOMECH_API FENeoHookean : public FEElasticMaterial
{
public:
	FENeoHookean(FEModel* pfem) : FEElasticMaterial(pfem) {}

public:
	FEParamDouble		m_E;	//!< Young's modulus
	FEParamDouble		m_v;	//!< Poisson's ratio

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
	DECLARE_FECORE_CLASS();
};
