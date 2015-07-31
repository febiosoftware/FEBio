#pragma once
#include "FEBoundaryCondition.h"
#include "FEGlobalVector.h"
#include "FETypes.h"

using namespace FECore;

class FESolver;

//-----------------------------------------------------------------------------
//! concentrated nodal force boundary condition

class FENodalForce : public FEBoundaryCondition
{
public:
	FENodalForce(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem){}

	void Serialize(DumpFile& ar);

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

	void Serialize(DumpFile& ar);

public:
	double	s;		// scale factor
	int		node;	// node number
	int		bc;		// dof
	int		lc;		// load curve
	double	r;		// initial value
	bool	br;		// flag for relative bc
};

//-----------------------------------------------------------------------------
//! rigid node

class FERigidNode : public FEBoundaryCondition
{
public:
	FERigidNode(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem){}

	void Serialize(DumpFile& ar);

public:
	int	nid;	// node number
	int	rid;	// rigid body number
};

//-----------------------------------------------------------------------------
//! fixed rigid body constraint
class FERigidBodyFixedBC : public FEBoundaryCondition
{
public:
	FERigidBodyFixedBC(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem){}

	void Serialize(DumpFile& ar);

public:
	int		id;	//!< rigid body ID
	int		bc;	//!< constrained dof
};

//-----------------------------------------------------------------------------
//! rigid body displacement

class FERigidBodyDisplacement : public FEBoundaryCondition
{
public:
	FERigidBodyDisplacement(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem) { ref= 0.0; brel = false; }

	double Value();

	void Serialize(DumpFile& ar);

public:
	int		id;		//!< rigid body id
	int		bc;		//!< displacement direction
	int		lc;		//!< load curve number
	double	sf;		//!< scale factor
	double	ref;	//!< reference value for relative displacement
	bool	brel;	//!< relative displacement flag
};

//-----------------------------------------------------------------------------
//! rigid body initial velocity
class FERigidBodyVelocity : public FEBoundaryCondition
{
public:
	FERigidBodyVelocity(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem){}

public:
	int		id;	//!< rigid body ID
	vec3d	v;	//!< value
};

//-----------------------------------------------------------------------------
//! rigid body initial angular velocity
class FERigidBodyAngularVelocity : public FEBoundaryCondition
{
public:
	FERigidBodyAngularVelocity(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem){}

public:
	int		id;	//!< rigid body ID
	vec3d	w;	//!< value
};
