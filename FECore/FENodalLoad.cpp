/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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

REGISTER_SUPER_CLASS(FENodalLoad, FENODALLOAD_ID);

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FENodalLoad, FEModelLoad)
	ADD_PROPERTY(m_nodeSet, "node_set", FEProperty::Reference);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FENodalLoad::FENodalLoad(FEModel* pfem) : FEModelLoad(pfem), m_dofs(pfem)
{
	m_nodeSet = nullptr;
}

//-----------------------------------------------------------------------------
void FENodalLoad::Serialize(DumpStream& ar)
{
	FEModelComponent::Serialize(ar);
	if (ar.IsShallow()) return;
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
void FENodalLoad::Residual(FEGlobalVector& R, const FETimeInfo& tp)
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

		// assemble into residual
		for (int j=0; j<dofs; ++j)
			R.Assemble(nid, m_dofs[j], val[j]);
	}
}

//-----------------------------------------------------------------------------
void FENodalLoad::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	// Nothing to do here.
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FENodalDOFLoad, FENodalLoad)
	ADD_PARAMETER(m_dof, "dof", 0, "@dof_list");
	ADD_PARAMETER(m_scale, "scale");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FENodalDOFLoad::FENodalDOFLoad(FEModel* fem) : FENodalLoad(fem)
{
	m_dof = -1;
	m_scale = 0.0;
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
	val[0] = m_scale;
}

double FENodalDOFLoad::NodeValue(int n)
{
	return m_scale;
}
