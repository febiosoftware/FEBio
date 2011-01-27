#pragma once
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
// This class implements a poroelastic material that has a strain-dependent
// permeability which is orthotropic in the reference state, but exhibits
// further strain-induced anisotropy, according to the constitutive relation
// of Ateshian and Weiss (JBME 2010)

class FEPermRefOrtho :	public FEHydraulicPermeability
	{
	public:
		//! constructor
		FEPermRefOrtho();
		
		//! permeability
		mat3ds Permeability(FEMaterialPoint& pt);
		
		//! Tangent of permeability
		tens4ds Tangent_Permeability_Strain(FEMaterialPoint& mp);
		
		//! data initialization and checking
		void Init();
		
	public:
		double	m_perm0;		//!< permeability for I term
		double	m_perm1[3];		//!< permeability for b term
		double	m_perm2[3];		//!< permeability for b^2 term
		double	m_phi0;			//!< solid volume fraction in reference state
		double	m_M0;			//!< nonlinear exponential coefficient
		double	m_alpha0;		//!< nonlinear power exponent
		double	m_M[3];			//!< nonlinear exponential coefficient
		double	m_alpha[3];		//!< nonlinear power exponent
		
		// declare as registered
		DECLARE_REGISTERED(FEPermRefOrtho);
		
		// declare parameter list
		DECLARE_PARAMETER_LIST();
	};
