#include "stdafx.h"
#include "FENodalLoad.h"
#include "FENodeSet.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FENodalLoad, FEBoundaryCondition)
	ADD_PARAMETER(m_scale, "scale");
	ADD_PARAMETER(m_data , "value");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FENodalLoad::FENodalLoad(FEModel* pfem) : FEBoundaryCondition(pfem), m_data(FE_DOUBLE)
{
	m_scale = 1.0;
	m_dof = -1;
}

//-----------------------------------------------------------------------------
void FENodalLoad::Serialize(DumpStream& ar)
{
	FEBoundaryCondition::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_dof & m_item;
}

//-----------------------------------------------------------------------------
bool FENodalLoad::Init()
{
	return true;
}

//-----------------------------------------------------------------------------
void FENodalLoad::AddNode(int nid, double scale)
{
	m_item.push_back(nid);
	m_data.Add(scale);
}

//-----------------------------------------------------------------------------
void FENodalLoad::AddNodes(const FENodeSet& ns, double scale)
{
	int N = ns.size();
	for (int i = 0; i<N; ++i) AddNode(ns[i], scale);
}

//-----------------------------------------------------------------------------
void FENodalLoad::SetLoad(double s)
{
	m_scale = s;
}

//-----------------------------------------------------------------------------
//! Return the current value of the nodal load
double FENodalLoad::NodeValue(int n) const
{
	return m_scale*m_data.getValue(n);
}
