/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2026 University of Utah, The Trustees of Columbia University in
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
#include "FENodalDOFLoad.h"
#include "FENodeSet.h"
#include "FENode.h"
#include "FEMaterialPoint.h"

BEGIN_FECORE_CLASS(FENodalDOFLoad, FENodalLoad)
	ADD_PARAMETER(m_dof, "dof", 0, "$(dof_list)");
	ADD_PARAMETER(m_scale, "scale")->SetFlags(FE_PARAM_ADDLC | FE_PARAM_VOLATILE);
END_FECORE_CLASS();

FENodalDOFLoad::FENodalDOFLoad(FEModel* fem) : FENodalLoad(fem)
{
	m_dof = -1;
	m_scale = 0.0;

	m_dtscale = false;
}

void FENodalDOFLoad::SetDtScale(bool b)
{
	m_dtscale = b;
}

void FENodalDOFLoad::SetLoad(double s)
{
	m_scale = s;
}

bool FENodalDOFLoad::SetDofList(FEDofList& dofList)
{
	return dofList.AddDof(m_dof);
}

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
