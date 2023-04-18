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
#include "FEBCRigidDeformation.h"
#include "FEBioMech.h"
#include <FECore/FENodeSet.h>
#include <FECore/FENode.h>

BEGIN_FECORE_CLASS(FEBCRigidDeformation, FEPrescribedNodeSet)
	ADD_PARAMETER(m_rt, "pos");
	ADD_PARAMETER(m_qt, "rot");
END_FECORE_CLASS();

FEBCRigidDeformation::FEBCRigidDeformation(FEModel* fem) : FEPrescribedNodeSet(fem)
{
	m_rt = vec3d(0, 0, 0);
	m_qt = vec3d(0, 0, 0);

	// TODO: Can this be done in Init, since there is no error checking
	if (fem)
	{
		FEDofList dofs(fem);
		dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
		SetDOFList(dofs);
	}
}

void FEBCRigidDeformation::Activate()
{
	m_r0 = m_rt;
	FEPrescribedNodeSet::Activate();
}

void FEBCRigidDeformation::CopyFrom(FEBoundaryCondition* pbc)
{
	FEBCRigidDeformation* ps = dynamic_cast<FEBCRigidDeformation*>(pbc); assert(ps);
	m_rt = ps->m_rt;
	m_qt = ps->m_qt;
	CopyParameterListState(ps->GetParameterList());
}

void FEBCRigidDeformation::GetNodalValues(int nodelid, std::vector<double>& val)
{
	vec3d X = GetNodeSet()->Node(nodelid)->m_r0;

	quatd Q(m_qt);

	vec3d x = Q*(X - m_r0) + m_rt;
	vec3d u = x - X;

	val[0] = u.x;
	val[1] = u.y;
	val[2] = u.z;
}
