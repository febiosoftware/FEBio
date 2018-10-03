#pragma once
#include "FEBiphasicSolute.h"

//-----------------------------------------------------------------------------
// This class implements a material that has a strain-dependent
// diffusivity which is isotropic in the reference state, but exhibits
// strain-induced anisotropy, according to the constitutive relation
// of Ateshian and Weiss (JBME 2010)

class FEDiffRefIso : public FESoluteDiffusivity
{
public:
	//! constructor
	FEDiffRefIso(FEModel* pfem);
	
	//! free diffusivity
	double Free_Diffusivity(FEMaterialPoint& pt) override;

	//! Tangent of free diffusivity with respect to concentration
	double Tangent_Free_Diffusivity_Concentration(FEMaterialPoint& pt, const int isol) override;
		
	
	//! diffusivity
	mat3ds Diffusivity(FEMaterialPoint& pt) override;
	
	//! Tangent of diffusivity with respect to strain
	tens4ds Tangent_Diffusivity_Strain(FEMaterialPoint& mp) override;
	
	//! Tangent of diffusivity with respect to concentration
	mat3ds Tangent_Diffusivity_Concentration(FEMaterialPoint& mp, const int isol=0) override;
	
public:
	double	m_free_diff;	//!< free diffusivity
	double	m_diff0;		//!< diffusivity for I term
	double	m_diff1;		//!< diffusivity for b term
	double	m_diff2;		//!< diffusivity for b^2 term
	double	m_M;			//!< nonlinear exponential coefficient
	double	m_alpha;		//!< nonlinear power exponent
	
	// declare parameter list
	DECLARE_FECORE_CLASS();
};
