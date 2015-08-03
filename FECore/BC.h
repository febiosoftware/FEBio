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

	bool Init();

	void Serialize(DumpFile& ar);

public:
	double	s;		// scale factor
	int		node;	// node number
	int		bc;		// force direction
	int		lc;		// load curve
};

//-----------------------------------------------------------------------------
class FEFixedBC : public FEBoundaryCondition
{
public:
	FEFixedBC(FEModel* pfem);
	FEFixedBC(FEModel* pfem, int node, int dof);

	void Serialize(DumpFile& ar);

	void Activate();

public:
	int		m_node;
	int		m_dof;
};

//-----------------------------------------------------------------------------
//! prescribed boundary condition data

class FEPrescribedBC : public FEBoundaryCondition
{
public:
	FEPrescribedBC(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem){}

	void Serialize(DumpFile& ar);

	bool Init();

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

	void Activate();
	void Deactivate();

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

	bool Init();

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

	bool Init();

	double Value();

	void Serialize(DumpFile& ar);

//	void Activate();

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

	bool Init();

	void Activate();

public:
	int		m_rid;	//!< rigid body ID
	vec3d	m_vel;	//!< initial velocity
};

//-----------------------------------------------------------------------------
//! rigid body initial angular velocity
class FERigidBodyAngularVelocity : public FEBoundaryCondition
{
public:
	FERigidBodyAngularVelocity(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem){}

	bool Init();

	void Activate();

public:
	int		m_rid;	//!< rigid body ID
	vec3d	m_w;	//!< value
};
