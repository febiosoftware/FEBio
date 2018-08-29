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
	mat3ds Permeability(FEMaterialPoint& pt) override;
		
	//! Tangent of permeability
	tens4ds Tangent_Permeability_Strain(FEMaterialPoint& mp) override;
		
public:
	double	m_perm;			//!< permeability
		
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
