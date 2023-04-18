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
#include "FEElementSelectionCriterion.h"
#include <FECore/FEMesh.h>

BEGIN_FECORE_CLASS(FEElementSelectionCriterion, FEMeshAdaptorCriterion)
	ADD_PARAMETER(m_value, "value");
	ADD_PARAMETER(m_elemList, "element_list");
END_FECORE_CLASS();

FEElementSelectionCriterion::FEElementSelectionCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{
	m_value = 1.0;
}

FEMeshAdaptorSelection FEElementSelectionCriterion::GetElementSelection(FEElementSet* elemSet)
{
	FEMesh& mesh = GetMesh();
	FEMeshAdaptorSelection elemList;
	FEElementIterator it(&mesh, elemSet);
	for (; it.isValid(); ++it)
	{
		FEElement& el = *it;

		if (el.isActive())
		{
			// see if this element is in the element_list
			// TODO: This is really slow. Need to speed this up!
			int eid = el.GetID();
			int n = -1;
			for (int i = 0; i < m_elemList.size(); ++i)
			{
				if (m_elemList[i] == eid)
				{
					n = i;
					break;
				}
			}

			// set the value
			if (n != -1)
			{
				elemList.push_back(eid, m_value);
			}
		}
	}
	return elemList;
}
