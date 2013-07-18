#pragma once
#include "FEPoroElastic.h"

//-----------------------------------------------------------------------------
// This class implements a poroelastic material that has a strain-dependent
// permeability which is orthotropic in the reference state, but exhibits
// further strain-induced anisotropy, according to the constitutive relation
// of Ateshian and Weiss (JBME 2010)

class FERefTransIsoPerm :	public FEPoroElastic
	{
	public:
		//! constructor
		FERefTransIsoPerm();
		
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
		double	m_perm1T;		//!< transverse permeability for b term
		double	m_perm1A;		//!< axial permeability for b term
		double	m_perm2T;		//!< transverse permeability for b^2 term
		double	m_perm2A;		//!< axial permeability for b^2 term
		double	m_phi0;			//!< solid volume fraction in reference state
		double	m_M0;			//!< nonlinear exponential coefficient for I term
		double	m_MT;			//!< nonlinear exponential coefficient for transverse direction
		double	m_MA;			//!< nonlinear exponential coefficient for axial direction
		double	m_alpha0;		//!< nonlinear power exponent for I term
		double	m_alphaT;		//!< nonlinear power exponent for transverse direction
		double	m_alphaA;		//!< nonlinear power exponent for axial direction
		
		// declare as registered
		DECLARE_REGISTERED(FERefTransIsoPerm);
		
		// declare parameter list
		DECLARE_PARAMETER_LIST();
	};
