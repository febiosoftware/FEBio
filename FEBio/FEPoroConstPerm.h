#pragma once
#include "FEPoroElastic.h"

//-----------------------------------------------------------------------------
// This class implements a poroelastic material that has a constant permeability

class FEPoroConstPerm :	public FEPoroElastic
{
public:
	//! constructor
	FEPoroConstPerm();

	//! calculate fluid flux
	virtual vec3d Flux(FEMaterialPoint& pt);

	//! permeability
	virtual void Permeability(double k[3][3], FEMaterialPoint& pt);

	//! Tangent of permeability
	virtual tens4ds Tangent_Permeability(FEMaterialPoint& mp);

public:
	double	m_perm;			//!< permeability
	double	m_permv[3];		//!< permeability for diagonal tensor

	// declare as registered
	DECLARE_REGISTERED(FEPoroConstPerm);

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
