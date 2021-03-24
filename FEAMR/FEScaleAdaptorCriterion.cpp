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
#include "FEScaleAdaptorCriterion.h"
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>

BEGIN_FECORE_CLASS(FEScaleAdaptorCriterion, FEMeshAdaptorCriterion)
	ADD_PARAMETER(m_scale, "scale");
END_FECORE_CLASS();

FEScaleAdaptorCriterion::FEScaleAdaptorCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{
	m_scale = 1.0;
}

FEMeshAdaptorSelection FEScaleAdaptorCriterion::GetElementSelection(FEElementSet* elemSet)
{
	// get the mesh
	FEMesh& mesh = GetFEModel()->GetMesh();
	int NE = mesh.Elements();

	// the element list of elements that need to be refined
	FEMeshAdaptorSelection elemList;

	// loop over the elements
	FEElementIterator it(&mesh, elemSet);
	for (int i = 0; it.isValid(); ++it, ++i)
	{
		FEElement& el = *it;
		int ne = el.Nodes();
		int ni = el.GaussPoints();

		double v = 0.0;
		for (int j = 0; j < ni; ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			double s = m_scale(mp);
			v += s;
		}
		v /= (double)ni;

		elemList.push_back(el.GetID(), v);
	}

	// create the element list of elements that need to be refined
	return elemList;
}
