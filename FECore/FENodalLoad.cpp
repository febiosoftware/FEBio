/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FENodalLoad.h"
#include "FENodeSet.h"
#include "DumpStream.h"
#include "FENode.h"
#include "FEMaterialPoint.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FENodalLoad, FEModelLoad)
	ADD_PARAMETER(m_brelative, "relative");
//	ADD_PROPERTY(m_nodeSet, "node_set", FEProperty::Reference);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FENodalLoad::FENodalLoad(FEModel* pfem) : FEModelLoad(pfem), m_dofs(pfem)
{
	m_brelative = false;
	m_nodeSet = nullptr;
}

//-----------------------------------------------------------------------------
void FENodalLoad::Serialize(DumpStream& ar)
{
	FEModelComponent::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_dofs & m_nodeSet;
}

//-----------------------------------------------------------------------------
bool FENodalLoad::Init()
{
	// Make sure we have a node set
	if (m_nodeSet == nullptr) return false;

	// Get the DOF list from the derived class
	m_dofs.Clear();
	if (SetDofList(m_dofs) == false) return false;

	// make sure the dof list is not empty
	if (m_dofs.IsEmpty()) return false;

	// so far so good.
	return FEModelLoad::Init();
}

//-----------------------------------------------------------------------------
//! activation
void FENodalLoad::Activate()
{
	FEModelLoad::Activate();

	if (m_brelative)
	{
		int nodes = m_nodeSet->Size();
		int dofs = m_dofs.Size();
		if ((dofs == 0) || (nodes == 0)) return;

		m_rval.resize(nodes, vector<double>(dofs, 0.0));

		// get the current nodal loads
		for (int i = 0; i < nodes; ++i)
		{
			FENode& node = *m_nodeSet->Node(i);

			for (int j = 0; j < dofs; ++j)
			{
				m_rval[i][j] = -node.get_load(m_dofs[j]);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Get the DOF list
const FEDofList& FENodalLoad::GetDOFList() const
{
	return m_dofs;
}

//-----------------------------------------------------------------------------
void FENodalLoad::SetNodeSet(FENodeSet* ns)
{
	m_nodeSet = ns;
}

//-----------------------------------------------------------------------------
//! get the nodeset
FENodeSet* FENodalLoad::GetNodeSet()
{
	return m_nodeSet;
}

//-----------------------------------------------------------------------------
void FENodalLoad::LoadVector(FEGlobalVector& R)
{
	FENodeSet& nset = *m_nodeSet;
	int dofs = m_dofs.Size();
	vector<double> val(dofs, 0.0);
	int N = nset.Size();
	for (int i = 0; i<N; ++i)
	{
		int nid = nset[i];

		// get the nodal values
		GetNodalValues(i, val);

		// add relative values
		if (m_brelative)
		{
			for (int j = 0; j < dofs; ++j) val[j] += m_rval[i][j];
		}

		// assemble into residual
		for (int j=0; j<dofs; ++j)
			R.Assemble(nid, m_dofs[j], val[j]);
	}
}

//-----------------------------------------------------------------------------
void FENodalLoad::StiffnessMatrix(FELinearSystem& LS)
{
	// Nothing to do here.
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FENodalDOFLoad, FENodalLoad)
	ADD_PARAMETER(m_dof, "dof", 0, "$(dof_list)");
	ADD_PARAMETER(m_scale, "scale")->SetFlags(FE_PARAM_ADDLC | FE_PARAM_VOLATILE);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FENodalDOFLoad::FENodalDOFLoad(FEModel* fem) : FENodalLoad(fem)
{
	m_dof = -1;
	m_scale = 0.0;

	m_dtscale = false;
}

//-----------------------------------------------------------------------------
void FENodalDOFLoad::SetDtScale(bool b)
{
	m_dtscale = b;
}

//-----------------------------------------------------------------------------
void FENodalDOFLoad::SetLoad(double s)
{
	m_scale = s;
}

//-----------------------------------------------------------------------------
bool FENodalDOFLoad::SetDofList(FEDofList& dofList)
{
	return dofList.AddDof(m_dof);
}

//-----------------------------------------------------------------------------
//! Return the current value of the nodal load
void FENodalDOFLoad::GetNodalValues(int n, std::vector<double>& val)
{
	assert(val.size() == 1);
	const FENodeSet& nset = *GetNodeSet();
	int nid = nset[n];
	const FENode& node = *nset.Node(n);

	FEMaterialPoint mp;
	mp.m_r0 = node.m_r0;
	mp.m_index = n;

	val[0] = m_scale(mp);

	if (m_dtscale)
	{
		double dt = GetTimeInfo().timeIncrement;
		val[0] *= dt;
	}
}

double FENodalDOFLoad::NodeValue(int n)
{
	const FENodeSet& nset = *GetNodeSet();
	int nid = nset[n];
	const FENode& node = *nset.Node(n);

	FEMaterialPoint mp;
	mp.m_r0 = node.m_r0;
	mp.m_index = n;

	double val = m_scale(mp);

	if (m_dtscale)
	{
		double dt = GetTimeInfo().timeIncrement;
		val *= dt;
	}

	return val;
}
