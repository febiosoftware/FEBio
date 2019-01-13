#pragma once
#include "FEBiphasic.h"

//-----------------------------------------------------------------------------
// This class implements a poroelastic material that has a strain-dependent
// permeability according to the constitutive relation of Holmes & Mow (JB 1990)

class FECORE_API FEPermHolmesMow :	public FEHydraulicPermeability
{
public:
	//! constructor
	FEPermHolmesMow(FEModel* pfem);
		
	//! permeability
	mat3ds Permeability(FEMaterialPoint& pt) override;
		
	//! Tangent of permeability
	tens4ds Tangent_Permeability_Strain(FEMaterialPoint& mp) override;
		
public:
	double	m_perm;			//!< permeability
	double	m_M;			//!< nonlinear exponential coefficient
	double	m_alpha;		//!< nonlinear power exponent
		
	// declare parameter list
	DECLARE_FECORE_CLASS();
};
