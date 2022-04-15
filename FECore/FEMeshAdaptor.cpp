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
#include "FEMeshAdaptor.h"
#include "FESolidDomain.h"
#include "FEElementList.h"
#include "FEModel.h"

FEMeshAdaptor::FEMeshAdaptor(FEModel* fem) : FEStepComponent(fem)
{
	m_elemSet = nullptr;
}

void FEMeshAdaptor::SetElementSet(FEElementSet* elemSet)
{
	m_elemSet = elemSet;
}

FEElementSet* FEMeshAdaptor::GetElementSet()
{
	return m_elemSet;
}

void FEMeshAdaptor::UpdateModel()
{
	FEModel& fem = *GetFEModel();
	fem.Reactivate();
}

// helper function for projecting integration point data to nodes
void projectToNodes(FEMesh& mesh, std::vector<double>& nodeVals, std::function<double(FEMaterialPoint& mp)> f)
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
			si[k] = f(mp);
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

// helper function for projecting integration point data to nodes
void projectToNodes(FEDomain& dom, std::vector<double>& nodeVals, std::function<double(FEMaterialPoint& mp)> f)
{
	// temp storage 
	double si[FEElement::MAX_INTPOINTS];
	double sn[FEElement::MAX_NODES];

	// allocate nodeVals and create valence array (tag)
	int NN = dom.Nodes();
	vector<int> tag(NN, 0);
	nodeVals.assign(NN, 0.0);

	// loop over all elements
	int NE = dom.Elements();
	for (int i = 0; i < NE; ++i)
	{
		FEElement& e = dom.ElementRef(i);
		int ne = e.Nodes();
		int ni = e.GaussPoints();

		// get the integration point values
		for (int k = 0; k < ni; ++k)
		{
			FEMaterialPoint& mp = *e.GetMaterialPoint(k);
			si[k] = f(mp);
		}

		// project to nodes
		e.project_to_nodes(si, sn);

		for (int j = 0; j < ne; ++j)
		{
			nodeVals[e.m_lnode[j]] += sn[j];
			tag[e.m_node[j]]++;
		}
	}

	for (int i = 0; i < NN; ++i)
	{
		if (tag[i] > 0) nodeVals[i] /= (double)tag[i];
	}
}

