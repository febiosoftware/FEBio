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
#include "FEFluidRotationalVelocity.h"
#include "FECore/FEElemElemList.h"
#include "FECore/FEGlobalMatrix.h"
#include "FECore/FEGlobalVector.h"
#include "FECore/log.h"
#include "FECore/LinearSolver.h"
#include <FECore/FENode.h>
#include "FEBioFluid.h"

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidRotationalVelocity, FEPrescribedNodeSet)
	ADD_PARAMETER(m_w, "angular_speed");
	ADD_PARAMETER(m_n, "axis"         );
	ADD_PARAMETER(m_p, "origin"       );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidRotationalVelocity::FEFluidRotationalVelocity(FEModel* pfem) : FEPrescribedNodeSet(pfem)
{
    m_w = 0.0;
    m_n = vec3d(0,0,1);
    m_p = vec3d(0,0,0);

	// Set the dof list
	// TODO: Can this be done in Init, since  there is no error checking
	if (pfem)
	{
		FEDofList dofs(pfem);
		dofs.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY));
		SetDOFList(dofs);
	}
}

//-----------------------------------------------------------------------------
// copy data from another class
void FEFluidRotationalVelocity::CopyFrom(FEBoundaryCondition* pbc)
{
	FEFluidRotationalVelocity* rv = dynamic_cast<FEFluidRotationalVelocity*>(pbc);
	m_w = rv->m_w;
	m_n = rv->m_n;
	m_p = rv->m_p;
	CopyParameterListState(rv->GetParameterList());
}

//-----------------------------------------------------------------------------
// return nodal value
void FEFluidRotationalVelocity::GetNodalValues(int nodelid, std::vector<double>& val)
{
	vec3d n(m_n); n.unit();
	vec3d v = (n ^ m_r[nodelid])*m_w;
	val[0] = v.x;
	val[1] = v.y;
	val[2] = v.z;
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidRotationalVelocity::Init()
{
    // evaluate nodal radial positions
	vec3d n(m_n); n.unit();
	const FENodeSet& nset = *GetNodeSet();
    int N = nset.Size();
    m_r.resize(N,vec3d(0,0,0));
    for (int i=0; i<N; ++i) {
        vec3d x = nset.Node(i)->m_r0 - m_p;
        m_r[i] = x - n*(x*n);
    }

    return FEPrescribedNodeSet::Init();
}

//-----------------------------------------------------------------------------
//! serialization
void FEFluidRotationalVelocity::Serialize(DumpStream& ar)
{
	FEPrescribedNodeSet::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_r;
}
