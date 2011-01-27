#pragma once
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
// This class implements a material that has a strain-dependent
// diffusivity which is isotropic in the reference state, but exhibits
// strain-induced anisotropy, according to the constitutive relation
// of Ateshian and Weiss (JBME 2010)

class FEDiffRefIso : public FESoluteDiffusivity
	{
	public:
		//! constructor
		FEDiffRefIso();
		
		//! free diffusivity
		double Free_Diffusivity(FEMaterialPoint& pt);
		
		//! diffusivity
		mat3ds Diffusivity(FEMaterialPoint& pt);
		
		//! Tangent of diffusivity with respect to strain
		tens4ds Tangent_Diffusivity_Strain(FEMaterialPoint& mp);
		
		//! Tangent of diffusivity with respect to concentration
		mat3ds Tangent_Diffusivity_Concentration(FEMaterialPoint& mp);
		
		//! data initialization and checking
		void Init();
		
	public:
		double	m_free_diff;	//!< free diffusivity
		double	m_diff0;		//!< diffusivity for I term
		double	m_diff1;		//!< diffusivity for b term
		double	m_diff2;		//!< diffusivity for b^2 term
		double	m_phi0;			//!< solid volume fraction in reference state
		double	m_M;			//!< nonlinear exponential coefficient
		double	m_alpha;		//!< nonlinear power exponent
		
		// declare as registered
		DECLARE_REGISTERED(FEDiffRefIso);
		
		// declare parameter list
		DECLARE_PARAMETER_LIST();
	};
