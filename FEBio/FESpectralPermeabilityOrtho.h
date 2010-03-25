#pragma once
#include "FEPoroElastic.h"

//-----------------------------------------------------------------------------
// This class implements a poroelastic material that has a strain-dependent
// orthotropic permeability, according to the constitutive relation
// of Ateshian and Weiss (JBME 2010)

class FESpectralPermeabilityOrtho :	public FEPoroElastic
	{
	public:
		//! constructor
		FESpectralPermeabilityOrtho();
		
		//! calculate fluid flux
		virtual vec3d Flux(FEMaterialPoint& pt);
		
		//! permeability
		virtual void Permeability(double k[3][3], FEMaterialPoint& pt);
		
		//! Tangent of permeability
		virtual tens4ds Tangent_Permeability(FEMaterialPoint& mp);
		
		//! data initialization and checking
		void Init();
		
	public:
		double	m_perm[3];		//!< permeability along preferred directions
		double	m_phi0;			//!< fluid volume fraction in reference state
		double	m_M[3];			//!< nonlinear exponential coefficient
		double	m_alpha[3];		//!< nonlinear power exponent
		
		// declare as registered
		DECLARE_REGISTERED(FESpectralPermeabilityOrtho);
		
		// declare parameter list
		DECLARE_PARAMETER_LIST();
	};
