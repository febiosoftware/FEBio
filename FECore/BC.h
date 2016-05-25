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
//! prescribed boundary condition data
//! \todo Should I make a derived class for the relative prescribed BC's?
class FEPrescribedBC : public FEBoundaryCondition
{
	struct ITEM
	{
		int		nid;	// nodal ID
		double	ref;	// reference value (for relative BC's)
	};

public:
	FEPrescribedBC(FEModel* pfem);
	FEPrescribedBC(FEModel* pfem, const FEPrescribedBC& bc);

	void AddNode(int node, double scale = 1.0);
	void AddNodes(const FENodeSet& s, double scale = 1.0);

	int NodeID(int i) { return m_item[i].nid; }

	size_t Items() const { return m_item.size(); }

	void Serialize(DumpStream& ar);

	void Activate();

	void Deactivate();

	bool Init();

	double NodeValue(int n) const;

	void Update();

	void PrepStep(std::vector<double>& ui, bool brel = true);

public:
	FEPrescribedBC& SetScale(double s, int lc = -1);
	FEPrescribedBC& SetDOF(int dof) { m_dof = dof; return *this; }
	FEPrescribedBC& SetRelativeFlag(bool br) { m_br = br; return *this; }

	void SetNodeScale(int n, double s) { m_data.set(n, s); }

	double GetScaleFactor() const { return m_scale; }
	int GetDOF() const { return m_dof; }

private:
	int			m_dof;		//!< dof
	double		m_scale;	//!< overall scale factor
	bool		m_br;		//!< flag for relative bc
	FEDataArray	m_data;		//!< nodal displacement values

	vector<ITEM>	m_item;		//!< item list

	DECLARE_PARAMETER_LIST();
};
