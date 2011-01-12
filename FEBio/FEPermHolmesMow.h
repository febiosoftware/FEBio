#pragma once
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
// This class implements a poroelastic material that has a strain-dependent
// permeability according to the constitutive relation of Holmes & Mow (JB 1990)

class FEPermHolmesMow :	public FEHydraulicPermeability
	{
	public:
		//! constructor
		FEPermHolmesMow();
		
		//! permeability
		mat3ds Permeability(FEMaterialPoint& pt);
		
		//! Tangent of permeability
		tens4ds Tangent_Permeability_Strain(FEMaterialPoint& mp);
		
		//! data initialization and checking
		void Init();
		
	public:
		double	m_perm;			//!< permeability
		double	m_M;			//!< nonlinear exponential coefficient
		double	m_alpha;		//!< nonlinear power exponent
		
		// declare as registered
		DECLARE_REGISTERED(FEPermHolmesMow);
		
		// declare parameter list
		DECLARE_PARAMETER_LIST();
	};
