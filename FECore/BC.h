#pragma once
#include "FEBoundaryCondition.h"

//-----------------------------------------------------------------------------
//! concentrated nodal force boundary condition

class FENodalForce : public FEBoundaryCondition
{
public:
	double	s;		// scale factor
	int		node;	// node number
	int		bc;		// force direction
	int		lc;		// load curve
};

//-----------------------------------------------------------------------------
//! prescribed boundary condition data

class FEPrescribedBC : public FEBoundaryCondition
{
public:
	double	s;		// scale factor
	int		node;	// node number
	int		bc;		// dof
	int		lc;		// load curve
};

//-----------------------------------------------------------------------------
//! rigid node

class FERigidNode : public FEBoundaryCondition
{
public:
	int	nid;	// node number
	int	rid;	// rigid body number
};
