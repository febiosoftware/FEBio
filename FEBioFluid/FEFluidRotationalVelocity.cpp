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
#include "FEFluidRotationalVelocity.h"
#include "FECore/FEModel.h"
#include "FECore/FEElemElemList.h"
#include "FECore/FEGlobalMatrix.h"
#include "FECore/FEGlobalVector.h"
#include "FECore/log.h"
#include "FECore/LinearSolver.h"

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidRotationalVelocity, FEPrescribedBC)
	ADD_PARAMETER(m_w, "angular_speed");
	ADD_PARAMETER(m_n, "axis"         );
	ADD_PARAMETER(m_p, "origin"       );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidRotationalVelocity::FEFluidRotationalVelocity(FEModel* pfem) : FEPrescribedBC(pfem)
{
    m_w = 0.0;
    m_n = vec3d(0,0,1);
    m_p = vec3d(0,0,0);
    
	vector<int> dofs(3);
    dofs[0] = pfem->GetDOFIndex("wx");
	dofs[1] = pfem->GetDOFIndex("wy");
	dofs[2] = pfem->GetDOFIndex("wz");
	SetDOFList(dofs);
}

//-----------------------------------------------------------------------------
// copy data from another class
void FEFluidRotationalVelocity::CopyFrom(FEPrescribedBC* pbc)
{
	FEFluidRotationalVelocity* rv = dynamic_cast<FEFluidRotationalVelocity*>(pbc);
	m_w = rv->m_w;
	m_n = rv->m_n;
	m_p = rv->m_p;
	CopyParameterListState(rv->GetParameterList());
}

//-----------------------------------------------------------------------------
// return nodal value
void FEFluidRotationalVelocity::NodalValues(int nodelid, std::vector<double>& val)
{
	vec3d v = (m_n ^ m_r[nodelid])*m_w;
	val[0] = v.x;
	val[1] = v.y;
	val[2] = v.z;
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidRotationalVelocity::Init()
{
    FEModelComponent::Init();
    
    m_n.unit();
    
    // evaluate nodal radial positions
	const FENodeSet& nset = GetNodeSet();
    int N = nset.size();
    m_r.resize(N,vec3d(0,0,0));
    for (int i=0; i<N; ++i) {
        vec3d x = nset.Node(i)->m_r0 - m_p;
        m_r[i] = x - m_n*(x*m_n);
    }
    return true;
}
