#pragma once
#include "FEBiphasic.h"

//-----------------------------------------------------------------------------
// This class implements a poroelastic material that has a strain-dependent
// permeability which is isotropic in the reference state, but exhibits
// strain-induced anisotropy, according to the constitutive relation
// of Ateshian and Weiss (JBME 2010)

class FEPermRefIso :	public FEHydraulicPermeability
{
public:
	//! constructor
	FEPermRefIso(FEModel* pfem);
		
	//! permeability
	mat3ds Permeability(FEMaterialPoint& pt) override;
		
	//! Tangent of permeability
	tens4ds Tangent_Permeability_Strain(FEMaterialPoint& mp) override;
		
	//! data initialization and checking
	bool Validate() override;
		
public:
	double	m_perm0;		//!< permeability for I term
	double	m_perm1;		//!< permeability for b term
	double	m_perm2;		//!< permeability for b^2 term
	double	m_phi0;			//!< solid volume fraction in reference state
	double	m_M;			//!< nonlinear exponential coefficient
	double	m_alpha;		//!< nonlinear power exponent
		
	// declare parameter list
	DECLARE_FECORE_CLASS();
};
