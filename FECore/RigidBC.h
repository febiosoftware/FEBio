#pragma once
#include "FEBoundaryCondition.h"

class FENodeSet;

//-----------------------------------------------------------------------------
//! rigid node set
class FECORE_API FERigidNodeSet : public FEBoundaryCondition
{
public:
	enum SHELL_BC {
		HINGED_SHELL,
		CLAMPED_SHELL
	};

public:
	FERigidNodeSet(FEModel* pfem);
	FERigidNodeSet(const FERigidNodeSet& rs);
	void operator = (const FERigidNodeSet& rs);

	bool Init() override;

	void Serialize(DumpStream& ar) override;

	void Activate() override;
	void Deactivate() override;

	int GetRigidID() const { return m_rid; }
	void SetRigidID(int rid) { m_rid = rid; }

	void SetNodeSet(FENodeSet& ns);

	void AddNode(int nid);

	void SetShellBC(SHELL_BC bc);

public:
	vector<int>		m_node;	// node number
	int				m_rid;	// rigid body number

private: // parameters
	int	m_nshellBC;		//!< flag defining how shells are attached (0=hinged, 1=clamped)

	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! fixed rigid body constraint
class FECORE_API FERigidBodyFixedBC : public FEBoundaryCondition
{
public:
	FERigidBodyFixedBC(FEModel* pfem);

	bool Init();

	void Serialize(DumpStream& ar);

	void Activate();

	void Deactivate();

public:
	int		id;	//!< rigid body ID
	int		bc;	//!< constrained dof

private:
	bool	m_binit;
};

//-----------------------------------------------------------------------------
//! rigid body displacement

class FECORE_API FERigidBodyDisplacement : public FEBoundaryCondition
{
public:
	FERigidBodyDisplacement(FEModel* pfem);

	bool Init();

	double Value();

	void Serialize(DumpStream& ar);

	void Activate();

	void Deactivate();

public:
	int		id;		//!< rigid body id
	int		bc;		//!< displacement direction
	int		lc;		//!< load curve number
	double	sf;		//!< scale factor
	double	ref;	//!< reference value for relative displacement
	bool	brel;	//!< relative displacement flag

private:
	bool	m_binit;	//!init flag
};

//-----------------------------------------------------------------------------
//! rigid body initial velocity
class FECORE_API FERigidBodyVelocity : public FEBoundaryCondition
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
class FECORE_API FERigidBodyAngularVelocity : public FEBoundaryCondition
{
public:
	FERigidBodyAngularVelocity(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem){}

	bool Init();

	void Activate();

public:
	int		m_rid;	//!< rigid body ID
	vec3d	m_w;	//!< value
};
