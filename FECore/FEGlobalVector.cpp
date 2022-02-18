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
#include "FEGlobalVector.h"
#include "vec3d.h"
#include "FEModel.h"

//-----------------------------------------------------------------------------
FEGlobalVector::FEGlobalVector(FEModel& fem, vector<double>& R, vector<double>& Fr) : m_fem(fem), m_R(R), m_Fr(Fr)
{
}

//-----------------------------------------------------------------------------
FEGlobalVector::~FEGlobalVector()
{

}

//-----------------------------------------------------------------------------
void FEGlobalVector::Assemble(vector<int>& en, vector<int>& elm, vector<double>& fe, bool bdom)
{
	vector<double>& R = m_R;

	// assemble the element residual into the global residual
	int ndof = (int)fe.size();
	for (int i=0; i<ndof; ++i)
	{
		int I = elm[i];
		if ( I >= 0) {
#pragma omp atomic
			R[I] += fe[i];
		}
// TODO: Find another way to store reaction forces
		else if (-I-2 >= 0) {
#pragma omp atomic
			m_Fr[-I-2] -= fe[i];
		}
	}
}

//-----------------------------------------------------------------------------
//! \todo This function does not add to m_Fr. Is this a problem?
void FEGlobalVector::Assemble(vector<int>& lm, vector<double>& fe)
{
	vector<double>& R = m_R;
	const int n = (int) lm.size();
	for (int i=0; i<n; ++i)
	{
		int nid = lm[i];
		if (nid >= 0) {
#pragma omp atomic
			R[nid] += fe[i];
		}
	}
}

//-----------------------------------------------------------------------------
//! assemble a nodel value
void FEGlobalVector::Assemble(int nodeId, int dof, double f)
{
	// get the equation number
	FENode& node = m_fem.GetMesh().Node(nodeId);
	int n = node.m_ID[dof];

	// assemble into global vector
	if (n >= 0) {
#pragma omp atomic
		m_R[n] += f;
	}
}
