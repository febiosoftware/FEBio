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
#include "FEMeshAdaptor.h"
#include "FEDomain.h"
#include "FESolidDomain.h"
#include <FECore/FEModel.h>
#include <FECore/FEElement.h>
#include <FECore/FEElementList.h>

REGISTER_SUPER_CLASS(FEMeshAdaptor, FEMESHADAPTOR_ID);

FEMeshAdaptor::FEMeshAdaptor(FEModel* fem) : FECoreBase(fem)
{

}

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
std::vector<pair<int, double> > FEMeshAdaptorCriterion::GetElementList()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	int nselected = 0;
	int nelem = 0;
	vector< pair<int, double> > elem;
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		int NE = dom.Elements();
		for (int j = 0; j < NE; ++j, ++nelem)
		{
			FEElement& el = dom.ElementRef(j);
			if (el.isActive())
			{
				double elemVal = 0;
				bool bselect = Check(el, elemVal);
				if (bselect)
				{
					elem.push_back(pair<int, double>(nelem, elemVal));
					nselected++;
				}
			}
		}
	}

	std::vector<pair<int, double> > selectedElement;
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
			selectedElement[i].first = elem[i].first;
			selectedElement[i].second = 0.;
		}
	}

	return selectedElement;
}

bool FEMeshAdaptorCriterion::Check(FEElement& el, double& elemVal)
{
	return false;
}

BEGIN_FECORE_CLASS(FEMaxVolumeCriterion, FEMeshAdaptorCriterion)
	ADD_PARAMETER(m_maxVolume, "max_vol");
END_FECORE_CLASS();

FEMaxVolumeCriterion::FEMaxVolumeCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{
	m_maxVolume = 0.0;
}

bool FEMaxVolumeCriterion::Check(FEElement& el, double& elemVal)
{
	FESolidDomain* dom = dynamic_cast<FESolidDomain*>(el.GetMeshPartition());
	if (dom == nullptr) return false;

	elemVal = dom->Volume(dynamic_cast<FESolidElement&>(el));

	return (elemVal >= m_maxVolume);
}

BEGIN_FECORE_CLASS(FEMaxVariableCriterion, FEMeshAdaptorCriterion)
	ADD_PARAMETER(m_maxValue, "max_value");
	ADD_PARAMETER(m_dof, "dof");
END_FECORE_CLASS();

FEMaxVariableCriterion::FEMaxVariableCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{
	m_maxValue = 0.0;
	m_dof = -1;
}

bool FEMaxVariableCriterion::Check(FEElement& el, double& elemVal)
{
	if (m_dof == -1) return false;

	FESolidDomain* dom = dynamic_cast<FESolidDomain*>(el.GetMeshPartition());
	if (dom == nullptr) return false;

	FEMesh& mesh = *dom->GetMesh();
	double maxVal = -1e99;
	for (int i = 0; i < el.Nodes(); ++i)
	{
		double vi = mesh.Node(el.m_node[i]).get(m_dof);
		if (vi > maxVal) maxVal = vi;
	}
	elemVal = maxVal;
	return (elemVal >= m_maxValue);
}

BEGIN_FECORE_CLASS(FEElementSelectionCriterion, FEMeshAdaptorCriterion)
	ADD_PARAMETER(m_elemList, "element_list");
END_FECORE_CLASS();

FEElementSelectionCriterion::FEElementSelectionCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{

}

std::vector<pair<int, double> > FEElementSelectionCriterion::GetElementList()
{
	vector<pair<int, double> > elemList(m_elemList.size(), pair<int, double>(-1, 0.0));
	FEMesh& mesh = GetFEModel()->GetMesh();
	for (int i = 0; i < (int)m_elemList.size(); ++i)
	{
		elemList[i].first = m_elemList[i] - 1;
	}
	return elemList;
}

BEGIN_FECORE_CLASS(FEDomainErrorCriterion, FEMeshAdaptorCriterion)
	ADD_PARAMETER(m_pct, "error");
END_FECORE_CLASS();

FEDomainErrorCriterion::FEDomainErrorCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{
	m_pct = 0.0;

	// set sort on by default
	SetSort(true);
}

std::vector<pair<int, double> > FEDomainErrorCriterion::GetElementList()
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
			double sj = GetMaterialPointValue(*el.GetMaterialPoint(j));
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
			double sj = GetMaterialPointValue(*el.GetMaterialPoint(j));

			double snj = el.Evaluate(ev, j);

			double err = fabs(sj - snj) / (smax - smin);
			if (err > max_err) max_err = err;
		}

		// see if it's too large
		if (max_err > m_pct)
		{
			double f = m_pct / max_err;
			elemList.push_back(pair<int, double>(i, f));
		}
	}

	// create the element list of elements that need to be refined
	return elemList;
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
