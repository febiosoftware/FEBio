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
#include "FENodalFluidFlux.h"
#include <FECore/FENodeSet.h>
#include <FECore/FEMaterialPoint.h>
#include <FECore/FENode.h>
#include "FEBioMix.h"

BEGIN_FECORE_CLASS(FENodalFluidFlux, FENodalLoad)
	ADD_PARAMETER(m_w, "value");
END_FECORE_CLASS();

FENodalFluidFlux::FENodalFluidFlux(FEModel* fem) : FENodalLoad(fem)
{
	m_w = 0.0;
}

bool FENodalFluidFlux::SetDofList(FEDofList& dofList)
{
	return dofList.AddVariable(FEBioMix::GetVariableName(FEBioMix::FLUID_PRESSURE));
}

void FENodalFluidFlux::GetNodalValues(int inode, std::vector<double>& val)
{
	assert(val.size() == 1);
	const FENodeSet& nset = *GetNodeSet();
	int nid = nset[inode];
	const FENode& node = *nset.Node(inode);

	FEMaterialPoint mp;
	mp.m_r0 = node.m_r0;
	mp.m_index = inode;

	const FETimeInfo& tp = GetTimeInfo();

	// for consistency with evaluation of residual and stiffness matrix
	val[0] = m_w(mp)* tp.timeIncrement;
}
