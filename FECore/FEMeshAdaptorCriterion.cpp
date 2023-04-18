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
#include "FEMeshAdaptorCriterion.h"
#include "FEMeshAdaptor.h"
#include "FEModel.h"
#include "FEMesh.h"
#include <algorithm>

//=============================================================================
void FEMeshAdaptorSelection::Sort(FEMeshAdaptorSelection::SortFlag sortFlag)
{
	if (sortFlag == SORT_DECREASING)
	{
		std::sort(m_itemList.begin(), m_itemList.end(), [](Item& e1, Item& e2) {
			return (e1.m_elemValue > e2.m_elemValue);
		});
	}
	else
	{
		std::sort(m_itemList.begin(), m_itemList.end(), [](Item& e1, Item& e2) {
			return (e1.m_elemValue < e2.m_elemValue);
		});
	}
}

//=============================================================================
BEGIN_FECORE_CLASS(FEMeshAdaptorCriterion, FEModelComponent)
END_FECORE_CLASS();

FEMeshAdaptorCriterion::FEMeshAdaptorCriterion(FEModel* fem) : FEModelComponent(fem)
{
}

// return a list of elements that satisfy the criterion
FEMeshAdaptorSelection FEMeshAdaptorCriterion::GetElementSelection(FEElementSet* elemSet)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// allocate iterator
	FEElementIterator it(&mesh, elemSet);

	// loop over elements
	FEMeshAdaptorSelection selectedElements;
	for (;it.isValid(); ++it)
	{
		FEElement& el = *it;
		if (el.isActive())
		{
			// evaluate element average
			bool bvalid = true;
			double elemVal = 0.0;
			int ni = el.GaussPoints();
			for (int i = 0; i < ni; ++i)
			{
				double vali = 0.0;
				bool b = GetMaterialPointValue(*el.GetMaterialPoint(i), vali);
				if (b) elemVal += vali;
				bvalid = (bvalid && b);
			}
			elemVal /= (double)ni;

			if (bvalid)
			{
				int nid = el.GetID();
				selectedElements.push_back(nid, elemVal);
			}
		}
	}

	return selectedElements;
}

bool FEMeshAdaptorCriterion::GetMaterialPointValue(FEMaterialPoint& mp, double& elemVal)
{
	return false;
}
