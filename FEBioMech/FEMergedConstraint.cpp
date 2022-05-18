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
#include "FEMergedConstraint.h"
#include <FECore/FEModel.h>
#include <FECore/FELinearConstraintManager.h>
#include <FECore/FELinearConstraint.h>

FEMergedConstraint::FEMergedConstraint(FEModel& fem) : m_fem(fem)
{
}

FEMergedConstraint::~FEMergedConstraint()
{
}

bool FEMergedConstraint::Merge(FEFacetSet* surf1, FEFacetSet* surf2, const vector<int>& dofList)
{
	// extract the nodes from the surfaces
	FENodeList set1 = surf1->GetNodeList();
	FENodeList set2 = surf2->GetNodeList();

	// find for each node on surface1 a corresponding node on surface 2 within tolerance
	// First, make sure that set2 is larger than set1
	if (set1.Size() > set2.Size()) return false;
	int N1 = set1.Size();
	int N2 = set2.Size();

	// make sure there is something to do
	if (N1 == 0) return true;
	if (dofList.size() == 0) return true;

	// alright, let's get going
	vector<int> tag(N1, -1);
	for (int i=0; i<N1; ++i)
	{
		// get the node position
		vec3d ri = set1.Node(i)->m_rt;

		// find the closest node
		int n = 0;
		double Dmin = (set2.Node(0)->m_rt - ri).norm2();
		for (int j=1; j<N2; ++j)
		{
			vec3d rj = set2.Node(j)->m_rt;
			double D2 = (ri - rj).norm2();
			if (D2 < Dmin)
			{
				Dmin = D2;
				n = j;
			}
		}

		// since this interface type assumes that the nodes match identically,
		// we check that the min distance is indeed very small
		if (Dmin > 1e-9) { 
			return false;
		}

		// store this node
		tag[i] = n;
	}

	// next, create the linear constraints 
	// get the linear constraint manager
	int ndofs = (int) dofList.size();
	FELinearConstraintManager& LCM = m_fem.GetLinearConstraintManager();
	for (int i=0; i<N1; ++i)
	{
		for (int j=0; j<ndofs; ++j)
		{
			int dof = dofList[j];
			FELinearConstraint* lc = fecore_alloc(FELinearConstraint, &m_fem);
			lc->SetParentDof(dof, set1[i]);
			lc->AddChildDof(dof, set2[tag[i]], 1.0);

			LCM.AddLinearConstraint(lc);
		}
	}

	// all done
	return true;
}
