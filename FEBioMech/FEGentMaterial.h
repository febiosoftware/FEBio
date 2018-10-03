#pragma once
#include <FEBioMech/FEElasticMaterial.h>
#include <FEBioMech/FEUncoupledMaterial.h>

//-----------------------------------------------------------------------------
// uncoupled Gent material
class FEGentMaterial : public FEUncoupledMaterial
{
public:
	//! constructor
	FEGentMaterial(FEModel* pfem);

	//! deviatoric Cauchy stress
	mat3ds DevStress(FEMaterialPoint& mp) override;

	//! Deviatoric spatial Tangent
	tens4ds DevTangent(FEMaterialPoint& mp) override;

private: // material parameters
	double	m_G;	//!< shear modulus
	double	m_Jm;	//!< Jm = Im - 3, where Im is max first invariant

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// compressible Gent material
class FECompressibleGentMaterial : public FEElasticMaterial
{
public:
	// constructor
	FECompressibleGentMaterial(FEModel* pfem);

	// Cauchy stress
	mat3ds Stress(FEMaterialPoint& mp) override;

	// spatial elasticity tangent
	tens4ds Tangent(FEMaterialPoint& mp) override;

public: // material parameter
	double	m_G;	//!< shear modulus
	double	m_K;	//!< bulk modulus
	double	m_Jm;	//!< Jm = Im - 3, where Im is max first invariant

	DECLARE_FECORE_CLASS();
};
