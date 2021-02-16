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
#include "FEDomainErrorCriterion.h"
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FEMeshAdaptor.h> // for projectToNodes

BEGIN_FECORE_CLASS(FEDomainErrorCriterion, FEMeshAdaptorCriterion)
	ADD_PARAMETER(m_pct, "error");
	ADD_PARAMETER(m_hmin, "hmin");
END_FECORE_CLASS();

FEDomainErrorCriterion::FEDomainErrorCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{
	m_pct = 0.0;

	m_hmin = 0;	// criterion not used by default

	// set sort on by default
	SetSort(true);
}

double MinEdgeLengthSqr(FEMesh* pm, FEElement& el)
{
	double h2 = 1e99;
	int ne = el.Nodes();
	for (int i = 0; i < ne; ++i)
	{
		for (int j = 0; j < ne; ++j)
		{
			if (i != j)
			{
				vec3d a = pm->Node(el.m_node[i]).m_r0;
				vec3d b = pm->Node(el.m_node[j]).m_r0;

				double L2 = (a - b).norm2();
				if (L2 < h2) h2 = L2;
			}
		}
	}
	return h2;
}

FEMeshAdaptorSelection FEDomainErrorCriterion::GetElementSelection(FEElementSet* elemSet)
{
	// get the mesh
	FEMesh& mesh = GetFEModel()->GetMesh();
	int NE = mesh.Elements();
	int NN = mesh.Nodes();

	// calculate the recovered nodal stresses
	vector<double> sn(NN);
	projectToNodes(mesh, sn, [this](FEMaterialPoint& mp) {
		return GetMaterialPointValue(mp);
	});

	// the element list of elements that need to be refined
	FEMeshAdaptorSelection elemList;

	// find the min and max stress values
	FEElementIterator it(&mesh, elemSet);
	double smin = 1e99, smax = -1e99;
	for (; it.isValid(); ++it)
	{
		FEElement& el = *it;
		int ni = el.GaussPoints();
		for (int j = 0; j < ni; ++j)
		{
			double sj = GetMaterialPointValue(*el.GetMaterialPoint(j));
			if (sj < smin) smin = sj;
			if (sj > smax) smax = sj;
		}
	}
	if (fabs(smin - smax) < 1e-12) return elemList;

	// calculate errors
	double ev[FEElement::MAX_NODES];
	it.reset();
	for (int i = 0; it.isValid(); ++it, ++i)
	{
		FEElement& el = *it;
		int ne = el.Nodes();
		int ni = el.GaussPoints();

		// check for minimal element size
		bool sizeCheck = true;
		if (m_hmin > 0)
		{
			double h2 = MinEdgeLengthSqr(&mesh, el);
			if (h2 < m_hmin*m_hmin) sizeCheck = false;
		}

		if (sizeCheck)
		{
			// get the nodal values
			for (int j = 0; j < ne; ++j)
			{
				ev[j] = sn[el.m_node[j]];
			}

			// evaluate element error
			double max_err = 0;
			for (int j = 0; j < ni; ++j)
			{
				double sj = GetMaterialPointValue(*el.GetMaterialPoint(j));

				double snj = el.Evaluate(ev, j);

				double err = fabs(sj - snj) / (smax - smin);
				if (err > max_err) max_err = err;
			}

			// see if it's too large
			if (max_err > m_pct)
			{
				double f = m_pct / max_err;
				elemList.push_back(i, f);
			}
		}
	}

	// create the element list of elements that need to be refined
	return elemList;
}
