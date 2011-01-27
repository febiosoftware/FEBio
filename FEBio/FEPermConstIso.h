#pragma once
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
// This class implements a poroelastic material that has a constant permeability

class FEPermConstIso :	public FEHydraulicPermeability
	{
	public:
		//! constructor
		FEPermConstIso();
		
		//! permeability
		mat3ds Permeability(FEMaterialPoint& pt);
		
		//! Tangent of permeability
		tens4ds Tangent_Permeability_Strain(FEMaterialPoint& mp);
		
		//! data initialization and checking
		void Init();
		
	public:
		double	m_perm;			//!< permeability
		
		// declare as registered
		DECLARE_REGISTERED(FEPermConstIso);
		
		// declare parameter list
		DECLARE_PARAMETER_LIST();
	};
