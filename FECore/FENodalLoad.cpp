#include "stdafx.h"
#include "FENodalLoad.h"
#include "FENodeSet.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FENodalLoad, FEBoundaryCondition)
	ADD_PARAMETER(m_scale, FE_PARAM_DOUBLE, "scale");
	ADD_PARAMETER(m_data, FE_PARAM_DATA_ARRAY, "value");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FENodalLoad::FENodalLoad(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem), m_data(FE_DOUBLE)
{
	m_scale = 1.0;
	m_dof = -1;
}

//-----------------------------------------------------------------------------
void FENodalLoad::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	FEBoundaryCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_dof << m_item;
	}
	else
	{
		ar >> m_dof >> m_item;
	}
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
void FENodalLoad::SetLoad(double s, int lc)
{
	m_scale = s;
	if (lc >= 0)
	{
		FEParam& p = *FEParamContainer::FindParameterFromData((void*)(&m_scale));
		p.SetLoadCurve(lc, m_scale);
	}
}

//-----------------------------------------------------------------------------
//! Return the current value of the nodal load
double FENodalLoad::NodeValue(int n) const
{
	return m_scale*m_data.getValue(n);
}
