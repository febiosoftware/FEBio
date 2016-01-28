#pragma once
#include "FEBiphasic.h"

//-----------------------------------------------------------------------------
// This class implements a poroelastic material that has a constant permeability

class FEPermConstIso :	public FEHydraulicPermeability
{
public:
	//! constructor
	FEPermConstIso(FEModel* pfem);
		
	//! permeability
	mat3ds Permeability(FEMaterialPoint& pt);
		
	//! Tangent of permeability
	tens4ds Tangent_Permeability_Strain(FEMaterialPoint& mp);
		
public:
	double	m_perm;			//!< permeability
		
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
