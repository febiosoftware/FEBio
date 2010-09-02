#pragma once
#include "FEPoroElastic.h"

//-----------------------------------------------------------------------------
// This class implements a poroelastic material that has a strain-dependent
// permeability which is isotropic in the reference state, but exhibits
// strain-induced anisotropy, according to the constitutive relation
// of Ateshian and Weiss (JBME 2010)

class FERefIsoPerm :	public FEPoroElastic
	{
	public:
		//! constructor
		FERefIsoPerm();
		
		//! calculate fluid flux
		virtual vec3d Flux(FEMaterialPoint& pt);
		
		//! permeability
		virtual void Permeability(double k[3][3], FEMaterialPoint& pt);
		
		//! Tangent of permeability
		virtual tens4ds Tangent_Permeability(FEMaterialPoint& mp);
		
		//! data initialization and checking
		void Init();
		
	public:
		double	m_perm0;		//!< permeability for I term
		double	m_perm1;		//!< permeability for b term
		double	m_perm2;		//!< permeability for b^2 term
		double	m_phi0;			//!< solid volume fraction in reference state
		double	m_M;			//!< nonlinear exponential coefficient
		double	m_alpha;		//!< nonlinear power exponent
		
		// declare as registered
		DECLARE_REGISTERED(FERefIsoPerm);
		
		// declare parameter list
		DECLARE_PARAMETER_LIST();
	};
