/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
#include "FEFixedFluidVelocity.h"

BEGIN_FECORE_CLASS(FEFixedFluidVelocity, FEFixedBC)
	ADD_PARAMETER(m_dof[0], "wx_dof")->setLongName("x-velocity");
	ADD_PARAMETER(m_dof[1], "wy_dof")->setLongName("y-velocity");
	ADD_PARAMETER(m_dof[2], "wz_dof")->setLongName("z-velocity");
END_FECORE_CLASS();

FEFixedFluidVelocity::FEFixedFluidVelocity(FEModel* fem) : FEFixedBC(fem)
{
	m_dof[0] = false;
	m_dof[1] = false;
	m_dof[2] = false;
}

bool FEFixedFluidVelocity::Init()
{
	vector<int> dofs;
	if (m_dof[0]) dofs.push_back(GetDOFIndex("wx"));
	if (m_dof[1]) dofs.push_back(GetDOFIndex("wy"));
	if (m_dof[2]) dofs.push_back(GetDOFIndex("wz"));
	SetDOFList(dofs);

	return FEFixedBC::Init();
}
