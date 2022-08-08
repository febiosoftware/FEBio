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
#include "FEPrescribedDOF.h"
#include "FENodeSet.h"
#include "DumpStream.h"
#include "FEMesh.h"
#include "log.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEPrescribedDOF, FEPrescribedNodeSet)
	ADD_PARAMETER(m_scale, "scale");
	ADD_PARAMETER(m_dof  , "dof", 0, "$(dof_list)");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEPrescribedDOF::FEPrescribedDOF(FEModel* pfem) : FEPrescribedNodeSet(pfem)
{
	m_scale = 0.0;
	m_dof = -1;
}

//-----------------------------------------------------------------------------
FEPrescribedDOF::FEPrescribedDOF(FEModel* pfem, int dof, FENodeSet* nset) : FEPrescribedNodeSet(pfem)
{
	m_scale = 0.0;
	SetNodeSet(nset);
	SetDOF(dof);
}

//-----------------------------------------------------------------------------
void FEPrescribedDOF::SetDOF(int ndof)
{
	m_dof = ndof;
}

//-----------------------------------------------------------------------------
bool FEPrescribedDOF::SetDOF(const char* szdof)
{
	int ndof = GetDOFIndex(szdof);
	assert(ndof >= 0);
	if (ndof < 0) return false;
	SetDOF(ndof);
	return true;
}

//-----------------------------------------------------------------------------
// Sets the displacement scale factor. An optional load curve index can be given
// of the load curve that will control the scale factor.
FEPrescribedDOF& FEPrescribedDOF::SetScale(double s, int lc)
{
	m_scale = s;
	if (lc >= 0)
	{
		AttachLoadController(&m_scale, lc);
	}
	return *this;
}

//-----------------------------------------------------------------------------
bool FEPrescribedDOF::Init()
{
	// set the dof first before calling base class
	if (m_dof < 0) return false;
	SetDOFList(m_dof);

	// don't forget to call the base class
	if (FEPrescribedNodeSet::Init() == false) return false;

	// make sure this is not a rigid node
	FEMesh& mesh = GetMesh();
	int NN = mesh.Nodes();
	const FENodeSet& nset = *GetNodeSet();
	for (size_t i = 0; i<nset.Size(); ++i)
	{
		int nid = nset[i];
		if ((nid < 0) || (nid >= NN)) return false;
		if (mesh.Node(nid).m_rid != -1)
		{
			feLogError("Rigid nodes cannot be prescribed.");
			return false;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEPrescribedDOF::GetNodalValues(int n, std::vector<double>& val)
{
	assert(val.size() == 1);
	const FENodeSet& nset = *GetNodeSet();
	int nid = nset[n];
	const FENode& node = *nset.Node(n);

	FEMaterialPoint mp;
	mp.m_r0 = node.m_r0;
	mp.m_index = n;

	val[0] = m_scale(mp);
}

//-----------------------------------------------------------------------------
void FEPrescribedDOF::CopyFrom(FEBoundaryCondition* pbc)
{
	FEPrescribedDOF* ps = dynamic_cast<FEPrescribedDOF*>(pbc); assert(ps);
	m_scale = ps->m_scale;
	CopyParameterListState(ps->GetParameterList());
}
