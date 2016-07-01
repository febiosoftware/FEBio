#pragma once
#include "FEBoundaryCondition.h"
#include "FEGlobalVector.h"
#include "FEDataArray.h"
#include "FETypes.h"

class FESolver;
class FENodeSet;

//-----------------------------------------------------------------------------
//! Nodal load boundary condition
class FENodalLoad : public FEBoundaryCondition
{
public:
	//! constructor
	FENodalLoad(FEModel* pfem);

	//! initialization
	bool Init();

	//! serialiation
	void Serialize(DumpStream& ar);

	//! Add a node to the node set
	void AddNode(int nid, double scale = 1.0);

	//! add a node set
	void AddNodes(const FENodeSet& ns, double scale = 1.0);

	//! number of nodes
	int Nodes() const { return (int) m_item.size(); }

	//! Node ID
	int NodeID(int n) const { return m_item[n]; }

	//! get nodal value
	double NodeValue(int n) const;

	//! get/set load 
	void SetLoad(double s, int lc = -1);
	double GetLoad() const { return m_scale; }

	//! get/set degree of freedom
	void SetDOF(int ndof) { m_dof = ndof; }
	int GetDOF() const { return m_dof; }

private:
	int		m_dof;		// degree of freedom index

	double			m_scale;	// applied load scale factor
	vector<int>		m_item;		// item list
	FEDataArray		m_data;		// nodal data

	DECLARE_PARAMETER_LIST();
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

	//! add a node set
	void AddNodes(const FENodeSet& ns);

	//! set the degree of freedom that will be fixed
	void SetDOF(int dof);

public:
	//! serialization
	void Serialize(DumpStream& ar);

	//! activation
	void Activate();

	//! deactivations
	void Deactivate();

public:
	vector<int>		m_node;		//!< node set
	int				m_dof;		//!< fixed degree of freedom
};

//-----------------------------------------------------------------------------
// base class for prescribed boundary conditions
class FEPrescribedBC : public FEBoundaryCondition
{
public:
	FEPrescribedBC(FEModel* pfem);

public:
	// implement these functions

	// assign a node set to the prescribed BC
	virtual void AddNodes(const FENodeSet& set) = 0;

	// This function is called when the solver needs to know the 
	// prescribed dof values. The brel flag indicates wheter the total 
	// value is needed or the value with respect to the current nodal dof value
	virtual void PrepStep(std::vector<double>& ui, bool brel = true) = 0;

	// This is called during nodal update and should be used to enforce the 
	// nodal degrees of freedoms
	virtual void Update() = 0;

	// copy data from another class
	virtual void CopyFrom(FEPrescribedBC* pbc) = 0;

	// Also implement the following functions.
	// These are already declared in base classes.
//  bool Init();
//  void Activate();
//  void Deactivate();
//  void Serialize(DumpStream& ar);
};

//-----------------------------------------------------------------------------
//! prescribed boundary condition data
//! \todo Should I make a derived class for the relative prescribed BC's?
class FEPrescribedDOF : public FEPrescribedBC
{
	struct ITEM
	{
		int		nid;	// nodal ID
		double	ref;	// reference value (for relative BC's)
	};

public:
	FEPrescribedDOF(FEModel* pfem);
	FEPrescribedDOF(FEModel* pfem, const FEPrescribedDOF& bc);

	void AddNode(int node, double scale = 1.0);
	void AddNodes(const FENodeSet& s, double scale);
	void AddNodes(const FENodeSet& s) { AddNodes(s, 1.0); }

	int NodeID(int i) { return m_item[i].nid; }

	size_t Items() const { return m_item.size(); }

public:
	void Serialize(DumpStream& ar);

	void Activate();

	void Deactivate();

	bool Init();

	void Update();

	void PrepStep(std::vector<double>& ui, bool brel = true);

	void CopyFrom(FEPrescribedBC* pbc);

public:
	FEPrescribedDOF& SetScale(double s, int lc = -1);
	FEPrescribedDOF& SetDOF(int dof) { m_dof = dof; return *this; }
	FEPrescribedDOF& SetRelativeFlag(bool br) { m_br = br; return *this; }

	void SetNodeScale(int n, double s) { m_data.set(n, s); }

	double GetScaleFactor() const { return m_scale; }
	int GetDOF() const { return m_dof; }

	double NodeValue(int n) const;

private:
	int			m_dof;		//!< dof
	double		m_scale;	//!< overall scale factor
	bool		m_br;		//!< flag for relative bc
	FEDataArray	m_data;		//!< nodal displacement values

	vector<ITEM>	m_item;		//!< item list

	DECLARE_PARAMETER_LIST();
};
