#pragma once
#include "FEBoundaryCondition.h"
#include "FEGlobalVector.h"
#include "FETypes.h"

using namespace FECore;

class FESolver;

//-----------------------------------------------------------------------------
//! Nodal load boundary condition
class FENodalLoad : public FEBoundaryCondition
{
public:
	FENodalLoad(FEModel* pfem);

	bool Init();

	void Serialize(DumpFile& ar);

	double Value();

public:
	double	m_s;		// scale factor
	int		m_node;	// node number
	int		m_bc;		// dof
	int		m_lc;		// load curve
};

//-----------------------------------------------------------------------------
//! This class represents a fixed degree of freedom
//! This boundary conditions sets the BC attribute of the nodes in the nodeset
//! to DOF_FIXED when activated.
class FEFixedBC : public FEBoundaryCondition
{
public:
	//! constructors
	FEFixedBC(FEModel* pfem);
	FEFixedBC(FEModel* pfem, int node, int dof);

	//! add a node to the node set
	void AddNode(int node);

	//! set the degree of freedom that will be fixed
	void SetDOF(int dof);

public:
	//! serialization
	void Serialize(DumpFile& ar);

	//! activation
	void Activate();

	//! deactivations
	void Deactivate();

public:
	vector<int>		m_node;		//!< node set
	int				m_dof;		//!< fixed degree of freedom
};

//-----------------------------------------------------------------------------
//! prescribed boundary condition data
//! \todo Should I make a derived class for the relative prescribed BC's?
class FEPrescribedBC : public FEBoundaryCondition
{
	struct ITEM
	{
		int		nid;	// nodal ID
		double	scale;	// nodal scale factor
		double	ref;	// reference value (for relative BC's)
	};

public:
	FEPrescribedBC(FEModel* pfem);
	FEPrescribedBC(FEModel* pfem, const FEPrescribedBC& bc);

	void AddNode(int node, double scale = 1.0);
	int NodeID(int i) { return m_item[i].nid; }

	size_t Items() const { return m_item.size(); }

	void Serialize(DumpFile& ar);

	void Activate();

	void Deactivate();

	bool Init();

	double NodeValue(int n) const;

	void Update();

	void PrepStep(std::vector<double>& ui, bool brel = true);

public:
	FEPrescribedBC& SetScale(double s) { m_scale = s; return *this; }
	FEPrescribedBC& SetDOF(int dof) { m_dof = dof; return *this; }
	FEPrescribedBC& SetRelativeFlag(bool br) { m_br = br; return *this; }
	FEPrescribedBC& SetLoadCurveIndex(int lc) { m_lc = lc; return *this; }

	double GetScaleFactor() const { return m_scale; }
	int GetDOF() const { return m_dof; }

private:
	double	m_scale;	//!< overall scale factor
	int		m_dof;		//!< dof
	int		m_lc;		//!< load curve
	bool	m_br;		//!< flag for relative bc

	vector<ITEM>	m_item;		//!< item list
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
	FERigidBodyFixedBC(FEModel* pfem);

	bool Init();

	void Serialize(DumpFile& ar);

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

class FERigidBodyDisplacement : public FEBoundaryCondition
{
public:
	FERigidBodyDisplacement(FEModel* pfem);

	bool Init();

	double Value();

	void Serialize(DumpFile& ar);

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
