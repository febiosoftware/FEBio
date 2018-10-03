#pragma once
#include "FEBoundaryCondition.h"
#include "FENodeDataMap.h"

//-----------------------------------------------------------------------------
class FENodeSet;

//-----------------------------------------------------------------------------
//! Nodal load boundary condition
class FECORE_API FENodalLoad : public FEBoundaryCondition
{
public:
	//! constructor
	FENodalLoad(FEModel* pfem);

	//! initialization
	bool Init() override;

	//! serialiation
	void Serialize(DumpStream& ar) override;

	//! Add a node to the node set
	void AddNode(int nid, double scale = 1.0);

	//! add a node set
	void AddNodes(const FENodeSet& ns, double scale = 1.0);

	//! number of nodes
	int Nodes() const { return (int)m_item.size(); }

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
	FENodeDataMap	m_data;		// nodal data

	DECLARE_FECORE_CLASS();
};

