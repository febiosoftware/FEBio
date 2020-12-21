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
#include "FEMeshAdaptorCriterion.h"
#include "FEMeshAdaptor.h"
#include "FEModel.h"
#include "FEMesh.h"
#include <algorithm>

REGISTER_SUPER_CLASS(FEMeshAdaptorCriterion, FEMESHADAPTORCRITERION_ID);

BEGIN_FECORE_CLASS(FEMeshAdaptorCriterion, FECoreBase)
	ADD_PARAMETER(m_sortList, "sort");
	ADD_PARAMETER(m_maxelem, "max_elems");
END_FECORE_CLASS();

FEMeshAdaptorCriterion::FEMeshAdaptorCriterion(FEModel* fem) : FECoreBase(fem)
{
	m_sortList = false;
	m_maxelem = 0;
}

void FEMeshAdaptorCriterion::SetSort(bool b)
{
	m_sortList = b;
}

void FEMeshAdaptorCriterion::SetMaxElements(int m)
{
	m_maxelem = m;
}

// return a list of elements that satisfy the criterion
FEMeshAdaptorSelection FEMeshAdaptorCriterion::GetElementSelection(FEElementSet* elemSet)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// allocate iterator
	FEElementIterator it(&mesh, elemSet);

	// loop over elements
	int nselected = 0;
	int nelem = 0;
	vector< pair<int, double> > elem;
	for (;it.isValid(); ++it)
	{
		FEElement& el = *it;
		if (el.isActive())
		{
			double elemVal = 0;
			bool bselect = Check(el, elemVal);
			if (bselect)
			{
				int nid = nelem;
				if (elemSet) nid = (*elemSet)[nelem];
				elem.push_back(pair<int, double>(nid, elemVal));
				nselected++;
			}
		}
		nelem++;
	}

	FEMeshAdaptorSelection selectedElement;
	if (nselected > 0)
	{
		// sort the list
		if (m_sortList) {
			std::sort(elem.begin(), elem.end(), [](pair<int, double>& e1, pair<int, double>& e2) {
				return e1.second > e2.second;
			});
		}

		int nelem = elem.size();
		if ((m_maxelem > 0) && (nelem > m_maxelem)) nelem = m_maxelem;

		selectedElement.resize(nelem);
		for (int i = 0; i < nelem; ++i)
		{
			selectedElement[i].m_elementIndex = elem[i].first;
			selectedElement[i].m_scaleFactor = 0.;
		}
	}

	return selectedElement;
}

bool FEMeshAdaptorCriterion::Check(FEElement& el, double& elemVal)
{
	return false;
}
