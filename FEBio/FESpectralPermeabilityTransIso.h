#pragma once
#include "FEPoroElastic.h"

//-----------------------------------------------------------------------------
// This class implements a poroelastic material that has a strain-dependent
// transversely isotropic permeability, according to the constitutive relation
// of Ateshian and Weiss (JBME 2010)

class FESpectralPermeabilityTransIso :	public FEPoroElastic
	{
	public:
		//! constructor
		FESpectralPermeabilityTransIso();
		
		//! calculate fluid flux
		virtual vec3d Flux(FEMaterialPoint& pt);
		
		//! permeability
		virtual void Permeability(double k[3][3], FEMaterialPoint& pt);
		
		//! Tangent of permeability
		virtual tens4ds Tangent_Permeability(FEMaterialPoint& mp);
		
		//! data initialization and checking
		void Init();
		
	public:
		double	m_permA;		//!< axial permeability
		double	m_permT;		//!< transverse permeability
		double	m_phi0;			//!< fluid volume fraction in reference state
		double	m_MA;			//!< nonlinear exponential coefficient along axial direction
		double	m_MT;			//!< nonlinear exponential coefficient in transverse plane
		double	m_alphaA;		//!< nonlinear power exponent along axial direction
		double	m_alphaT;		//!< nonlinear power exponent in transverse plane
		
		// declare as registered
		DECLARE_REGISTERED(FESpectralPermeabilityTransIso);
		
		// declare parameter list
		DECLARE_PARAMETER_LIST();
	};
