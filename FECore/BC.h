#pragma once
#include "FEBoundaryCondition.h"
using namespace FECore;

//-----------------------------------------------------------------------------
//! concentrated nodal force boundary condition

class FENodalForce : public FEBoundaryCondition
{
public:
	FENodalForce() : FEBoundaryCondition(FEBC_ID){}

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
	FEPrescribedBC() : FEBoundaryCondition(FEBC_ID){}

public:
	double	s;		// scale factor
	int		node;	// node number
	int		bc;		// dof
	int		lc;		// load curve
	double	r;		// initial value // GAA
	bool	br;		// flag for relative bc
};

//-----------------------------------------------------------------------------
//! rigid node

class FERigidNode : public FEBoundaryCondition
{
public:
	FERigidNode() : FEBoundaryCondition(FEBC_ID){}

public:
	int	nid;	// node number
	int	rid;	// rigid body number
};

//-----------------------------------------------------------------------------
//! rigid body force

class FERigidBodyForce : public FEBoundaryCondition
{
public:
	FERigidBodyForce() : FEBoundaryCondition(FEBC_ID){}

public:
	int		id;	// rigid body id
	int		bc;	// force direction
	int		lc;	// load curve number
	double	sf;	// scale factor
};

//-----------------------------------------------------------------------------
//! rigid body displacement

class FERigidBodyDisplacement : public FEBoundaryCondition
{
public:
	FERigidBodyDisplacement() : FEBoundaryCondition(FEBC_ID){}

public:
	int		id;	// rigid body id
	int		bc;	// displacement direction
	int		lc;	// load curve number
	double	sf;	// scale factor
};
