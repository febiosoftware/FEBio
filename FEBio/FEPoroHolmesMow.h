#pragma once
#include "FEPoroElastic.h"

//-----------------------------------------------------------------------------
// This class implements a poroelastic material that has a strain-dependent
// permeability according to the constitutive relation of Holmes & Mow (JB 1990)

class FEPoroHolmesMow :	public FEPoroElastic
{
public:
	//! constructor
	FEPoroHolmesMow();
		
	//! calculate fluid flux
	virtual vec3d Flux(FEMaterialPoint& pt);
		
	//! permeability
	virtual void Permeability(double k[3][3], FEMaterialPoint& pt);

	//! Tangent of permeability
	virtual tens4ds Tangent_Permeability(FEMaterialPoint& mp);

	//! data initialization and checking
	void Init();
	
public:
	double	m_perm;			//!< permeability
	double	m_permv[3];		//!< permeability for diagonal tensor
	double	m_phi0;			//!< solid volume fraction in reference state
	double	m_M;			//!< nonlinear exponential coefficient
	double	m_alpha;		//!< nonlinear power exponent
	
	// declare as registered
	DECLARE_REGISTERED(FEPoroHolmesMow);
		
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
