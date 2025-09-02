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
#include "FEInitialRigidKinematics.h"
#include "FEBioMech.h"
#include <FECore/FEMaterialPoint.h>
#include <FECore/FENode.h>

BEGIN_FECORE_CLASS(FEInitialRigidKinematics, FENodalIC)
	ADD_PARAMETER(m_v, "velocity")->setUnits(UNIT_VELOCITY);
	ADD_PARAMETER(m_w, "angular_velocity")->setUnits(UNIT_ANGULAR_VELOCITY);
	ADD_PARAMETER(m_c, "center_of_rotation")->setUnits(UNIT_LENGTH);
END_FECORE_CLASS();

FEInitialRigidKinematics::FEInitialRigidKinematics(FEModel* fem) : FENodalIC(fem)
{
	m_v = vec3d(0, 0, 0);
	m_w = vec3d(0, 0, 0);
	m_c = vec3d(0, 0, 0);
}

// initialization
bool FEInitialRigidKinematics::Init()
{
	FEDofList dofs(GetFEModel());
	if (dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::VELOCITY)) == false) return false;
	SetDOFList(dofs);
	return FENodalIC::Init();
}

// return the values for node i
void FEInitialRigidKinematics::GetNodalValues(int inode, std::vector<double>& values)
{
	assert(values.size() == 3);

	const FENodeSet& nset = *GetNodeSet();
	const FENode& node = *nset.Node(inode);

	vec3d r = node.m_rt;
	vec3d v = m_v + (m_w ^ (r - m_c));

	values[0] = v.x;
	values[1] = v.y;
	values[2] = v.z;
}
