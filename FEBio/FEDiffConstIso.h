#pragma once
#include "FEBiphasic.h"

//-----------------------------------------------------------------------------
// This class implements a material that has a constant diffusivity

class FEDiffConstIso :	public FESoluteDiffusivity
	{
	public:
		//! constructor
		FEDiffConstIso();
		
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
		double	m_diff;			//!< diffusivity
		
		// declare as registered
		DECLARE_REGISTERED(FEDiffConstIso);
		
		// declare parameter list
		DECLARE_PARAMETER_LIST();
	};
