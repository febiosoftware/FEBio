/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in
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
#include "FEStressErrorCriterion.h"
#include "FEElasticMaterial.h"
#include <FECore/FEModel.h>
#include <FECore/FEElement.h>
#include <FECore/FEElementList.h>

BEGIN_FECORE_CLASS(FEStressErrorCriterion, FEMeshAdaptorCriterion)
	ADD_PARAMETER(m_pct, FE_RANGE_OPEN(0.0, 1.0), "error");
END_FECORE_CLASS();

void project_stresses(FEMesh& mesh, vector<double>& nodeVals);

FEStressErrorCriterion::FEStressErrorCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{
	m_pct = 0.1;

	// set sort on by default
	SetSort(true);
}

std::vector<pair<int, double> > FEStressErrorCriterion::GetElementList()
{
	// get the mesh
	FEMesh& mesh = GetFEModel()->GetMesh();
	int NE = mesh.Elements();
	int NN = mesh.Nodes();

	// calculate the recovered stresses
	vector<double> sn(NN);
	project_stresses(mesh, sn);

	// the element list of elements that need to be refined
	std::vector<pair<int, double> > elemList;

	FEElementList EL(mesh);
	FEElementList::iterator it = EL.begin();
	// find the min and max stress values
	double smin = 1e99, smax = -1e99;
	for (int i = 0; i < NE; ++i, ++it)
	{
		FEElement& el = *it;
		int ni = el.GaussPoints();
		for (int j = 0; j < ni; ++j)
		{
			FEElasticMaterialPoint* ep = el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>();
			double sj = ep->m_s.effective_norm();

			if (sj < smin) smin = sj;
			if (sj > smax) smax = sj;
		}
	}
	if (fabs(smin - smax) < 1e-12) return elemList;

	// calculate errors
	double ev[FEElement::MAX_NODES];
	it = EL.begin();
	for (int i = 0; i < NE; ++i, ++it)
	{
		FEElement& el = *it;
		int ne = el.Nodes();
		int ni = el.GaussPoints();

		// get the nodal values
		for (int j = 0; j < ne; ++j)
		{
			ev[j] = sn[el.m_node[j]];
		}

		// evaluate element error
		double max_err = 0;
		for (int j = 0; j < ni; ++j)
		{
			FEElasticMaterialPoint* ep = el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>();
			double sj = ep->m_s.effective_norm();

			double snj = el.Evaluate(ev, j);

			double err = fabs(sj - snj) / (smax - smin);
			if (err > max_err) max_err = err;
		}

		// see if it's too large
		if (max_err > m_pct)
		{
			double f = m_pct / max_err;
			elemList.push_back(pair<int,double>(i, f));
		}
	}

	// create the element list of elements that need to be refined
	return elemList;
}


void project_stresses(FEMesh& mesh, vector<double>& nodeVals)
{
	// temp storage 
	double si[FEElement::MAX_INTPOINTS];
	double sn[FEElement::MAX_NODES];

	// allocate nodeVals and create valence array (tag)
	int NN = mesh.Nodes();
	vector<int> tag(NN, 0);
	nodeVals.assign(NN, 0.0);

	// loop over all elements
	int NE = mesh.Elements();
	FEElementList EL(mesh);
	FEElementList::iterator it = EL.begin();
	for (int i = 0; i < NE; ++i, ++it)
	{
		FEElement& e = *it;
		int ne = e.Nodes();
		int ni = e.GaussPoints();

		// get the integration point values
		for (int k = 0; k < ni; ++k)
		{
			FEMaterialPoint& mp = *e.GetMaterialPoint(k);
			FEElasticMaterialPoint* ep = mp.ExtractData<FEElasticMaterialPoint>();

			mat3ds& s = ep->m_s;

			double v = s.effective_norm();

			si[k] = v;
		}

		// project to nodes
		e.project_to_nodes(si, sn);

		for (int j = 0; j < ne; ++j)
		{
			nodeVals[e.m_node[j]] += sn[j];
			tag[e.m_node[j]]++;
		}
	}

	for (int i = 0; i < NN; ++i)
	{
		if (tag[i] > 0) nodeVals[i] /= (double)tag[i];
	}
}
