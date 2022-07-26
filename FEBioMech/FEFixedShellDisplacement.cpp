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
#include "FEFixedShellDisplacement.h"

BEGIN_FECORE_CLASS(FEFixedShellDisplacement, FEFixedBC)
	ADD_PARAMETER(m_dof_sx, "sx_dof")->setLongName("sx-displacement");
	ADD_PARAMETER(m_dof_sy, "sy_dof")->setLongName("sy-displacement");
	ADD_PARAMETER(m_dof_sz, "sz_dof")->setLongName("sz-displacement");
END_FECORE_CLASS();

FEFixedShellDisplacement::FEFixedShellDisplacement(FEModel* fem) : FEFixedBC(fem)
{
	m_dof_sx = false;
	m_dof_sy = false;
	m_dof_sz = false;
}

bool FEFixedShellDisplacement::Init()
{
	vector<int> dofs;
	if (m_dof_sx) dofs.push_back(GetDOFIndex("sx"));
	if (m_dof_sy) dofs.push_back(GetDOFIndex("sy"));
	if (m_dof_sz) dofs.push_back(GetDOFIndex("sz"));

	SetDOFList(dofs);
	return FEFixedBC::Init();
}
