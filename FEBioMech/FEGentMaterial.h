#pragma once
#include <FEBioMech/FEElasticMaterial.h>

//-----------------------------------------------------------------------------
class FEGentMaterial : public FEElasticMaterial
{
public:
	// constructor
	FEGentMaterial(FEModel* pfem);

	// Cauchy stress
	mat3ds Stress(FEMaterialPoint& mp);

	// spatial elasticity tangent
	tens4ds Tangent(FEMaterialPoint& mp);

public: // material parameter
	double	m_G;	//!< shear modulus
	double	m_K;	//!< bulk modulus
	double	m_Jm;	//!< Jm = Im - 3, where Im is max first invariant

	DECLARE_PARAMETER_LIST();
};
