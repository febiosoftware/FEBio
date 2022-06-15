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
#include "FELinearConstraint.h"
#include "FEMesh.h"
#include "FEModel.h"
#include "DumpStream.h"

FELinearConstraint::DOF::DOF()
{
	node = dof = -1;
	val = 0.0;

	GetParameterList();
	AddParameter(val, "value");
}

void FELinearConstraint::DOF::Serialize(DumpStream& ar)
{
	FEParamContainer::Serialize(ar);
	ar & node & dof & val;
}

//-----------------------------------------------------------------------------
FELinearConstraint::FELinearConstraint() : FEModelComponent(nullptr)
{
	m_parentDof = nullptr;
	m_off = 0.0;
}

//-----------------------------------------------------------------------------
FELinearConstraint::FELinearConstraint(FEModel* pfem) : FEModelComponent(pfem) 
{
	m_parentDof = new DOF;
	m_off = 0.0;
}

//-----------------------------------------------------------------------------
FELinearConstraint::~FELinearConstraint()
{
	Clear();
}

//-----------------------------------------------------------------------------
// return offset
double FELinearConstraint::GetOffset() const
{
	return m_off;
}

//-----------------------------------------------------------------------------
void FELinearConstraint::Clear()
{
	if (m_parentDof) delete m_parentDof; m_parentDof = nullptr;
	for (size_t i = 0; i < m_childDof.size(); ++i) delete m_childDof[i];
	m_childDof.clear();
}

//-----------------------------------------------------------------------------
FELinearConstraint::FELinearConstraint(const FELinearConstraint& LC) : FEModelComponent(LC.GetFEModel())
{
	m_parentDof = nullptr;
	CopyFrom(LC);
}

//-----------------------------------------------------------------------------
void FELinearConstraint::CopyFrom(const FELinearConstraint& LC)
{
	Clear();
	if (LC.m_parentDof)
	{
		m_parentDof = new DOF;
		m_parentDof->node = LC.m_parentDof->node;
		m_parentDof->dof  = LC.m_parentDof->dof;
		m_parentDof->val  = LC.m_parentDof->val;
	}
	m_off = LC.m_off;
	int n = (int)LC.m_childDof.size();
	vector<DOF*>::const_iterator it = LC.m_childDof.begin();
	for (int i = 0; i < n; ++i, ++it)
	{
		DOF* d = new DOF;
		d->node = (*it)->node;
		d->dof  = (*it)->dof;
		d->val  = (*it)->val;
		m_childDof.push_back(d);
	}
}

//-----------------------------------------------------------------------------
void FELinearConstraint::SetParentDof(int dof, int node)
{
	if (m_parentDof == nullptr) m_parentDof = new DOF;
	m_parentDof->dof = dof;
	m_parentDof->node = node;
}

//-----------------------------------------------------------------------------
void FELinearConstraint::SetParentNode(int node)
{
	if (m_parentDof == nullptr) m_parentDof = new DOF;
	m_parentDof->node = node;
}

//-----------------------------------------------------------------------------
void FELinearConstraint::SetParentDof(int dof)
{
	if (m_parentDof == nullptr) m_parentDof = new DOF;
	m_parentDof->dof = dof;
}

//-----------------------------------------------------------------------------
// get the parent dof
int FELinearConstraint::GetParentDof() const
{
	return m_parentDof->dof;
}

//-----------------------------------------------------------------------------
int FELinearConstraint::GetParentNode() const
{
	return m_parentDof->node;
}

//-----------------------------------------------------------------------------
// get the child DOF
const FELinearConstraint::DOF& FELinearConstraint::GetChildDof(int n) const
{
	return *m_childDof[n];
}

//-----------------------------------------------------------------------------
size_t FELinearConstraint::Size() const
{
	return m_childDof.size();
}

//-----------------------------------------------------------------------------
FELinearConstraint::dof_iterator FELinearConstraint::begin()
{
	return m_childDof.begin();
}

//-----------------------------------------------------------------------------
void FELinearConstraint::AddChildDof(int dof, int node, double v)
{
	DOF* d = new DOF;
	d->dof = dof;
	d->node = node;
	d->val = v;
	m_childDof.push_back(d);
}

//-----------------------------------------------------------------------------
void FELinearConstraint::AddChildDof(FELinearConstraint::DOF* dof)
{
	m_childDof.push_back(dof);
}

//-----------------------------------------------------------------------------
// Initialization.
// Make sure the parent dof does not appear as a child dof
bool FELinearConstraint::Init()
{
	if (m_parentDof == nullptr) return false;
	int n = (int)m_childDof.size();
	for (int i=0; i<n; ++i)
	{
		DOF& childNode = *m_childDof[i];
		if ((childNode.node == m_parentDof->node) && (childNode.dof == m_parentDof->dof)) return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
// This is called during model activation (i.e. at the start of an analysis step)
// The parent dof is fixed in order to make sure that they are not assigned an equation number.
void FELinearConstraint::Activate()
{
	FEModelComponent::Activate();
	FEMesh& mesh = GetFEModel()->GetMesh();

	// we need the parent node to be fixed so that no equation is allocated
	FENode& node = mesh.Node(m_parentDof->node);
	node.set_bc(m_parentDof->dof, DOF_FIXED);
}

//-----------------------------------------------------------------------------
void FELinearConstraint::Deactivate()
{
	FEModelComponent::Deactivate();
	FEMesh& mesh = GetFEModel()->GetMesh();

	FENode& node = mesh.Node(m_parentDof->node);
	node.set_bc(m_parentDof->dof, DOF_OPEN);
}

//-----------------------------------------------------------------------------
void FELinearConstraint::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	if (ar.IsSaving())
	{
		m_parentDof->Serialize(ar);
		int n = (int)m_childDof.size();
		ar << n;
		vector<DOF*>::iterator it = m_childDof.begin();
		for (int i = 0; i < n; ++i, ++it) (*it)->Serialize(ar);
	}
	else
	{
		m_childDof.clear();
		if (m_parentDof == nullptr) m_parentDof = new DOF;
		m_parentDof->Serialize(ar);
		int n;
		ar >> n;
		for (int i=0; i<n; ++i)
		{
			DOF* dof = new DOF;
			dof->Serialize(ar);
			m_childDof.push_back(dof);
		}
	}
}
