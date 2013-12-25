#pragma once
#include "FEBoundaryCondition.h"
using namespace FECore;

//-----------------------------------------------------------------------------
//! concentrated nodal force boundary condition

class FENodalForce : public FEBoundaryCondition
{
public:
	FENodalForce(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem){}

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
	FEPrescribedBC(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem){}

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
	FERigidNode(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem){}

public:
	int	nid;	// node number
	int	rid;	// rigid body number
};

//-----------------------------------------------------------------------------
//! rigid body force

class FERigidBodyForce : public FEBoundaryCondition
{
public:
	FERigidBodyForce(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem){}

public:
	int		id;	// rigid body id
	int		bc;	// force direction
	int		lc;	// load curve number
	double	sf;	// scale factor
};

//-----------------------------------------------------------------------------
//! fixed rigid body constraint
class FERigidBodyFixedBC : public FEBoundaryCondition
{
public:
	FERigidBodyFixedBC(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem){}

public:
	int		id;	//!< rigid body ID
	int		bc;	//!< constrained dof
};

//-----------------------------------------------------------------------------
//! rigid body displacement

class FERigidBodyDisplacement : public FEBoundaryCondition
{
public:
	FERigidBodyDisplacement(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem) { ref= 0.0; }

public:
	int		id;		//!< rigid body id
	int		bc;		//!< displacement direction
	int		lc;		//!< load curve number
	double	sf;		//!< scale factor
	double	ref;	//!< reference value
};
