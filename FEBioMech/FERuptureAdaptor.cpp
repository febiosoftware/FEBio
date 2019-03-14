#include "stdafx.h"
#include "FERuptureAdaptor.h"
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FEDomain.h>
#include "FEElasticMaterial.h"
#include <algorithm>

BEGIN_FECORE_CLASS(FERuptureAdaptor, FEMeshAdaptor)
	ADD_PARAMETER(m_maxStress, FE_RANGE_GREATER(0.0), "max_stress");
	ADD_PARAMETER(m_maxElems, "max_elems");
	ADD_PARAMETER(m_maxIters, "max_iters");
	ADD_PARAMETER(m_metric  , "metric");
END_FECORE_CLASS();

FERuptureAdaptor::FERuptureAdaptor(FEModel* fem) : FEMeshAdaptor(fem)
{
	m_maxStress = 0.0;
	m_maxElems = 0;
	m_maxIters = -1;
	m_metric = 0;
}

bool FERuptureAdaptor::Apply()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	int naug = fem.GetTime().augmentation;
	if ((m_maxIters >= 0) && (naug >= m_maxIters)) return true;

	int deactiveElems = 0;

	vector< pair<int, double> > elem;
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		int NE = dom.Elements();
		for (int j = 0; j < NE; ++j)
		{
			FEElement& el = dom.ElementRef(j);
			if (el.isActive())
			{

				bool bdeactivate = false;
				int nint = el.GaussPoints();
				double elemVal = 0;
				for (int n = 0; n < nint; ++n)
				{
					FEMaterialPoint* mp = el.GetMaterialPoint(n);
					FEElasticMaterialPoint* ep = mp->ExtractData<FEElasticMaterialPoint>();
					if (ep)
					{
						mat3ds& s = ep->m_s;

						switch (m_metric)
						{
						case 0: elemVal = s.effective_norm(); break;
						case 1: elemVal = fabs(s.zz()); break;
						default:
							break;
						}
						
						if (elemVal >= m_maxStress)
						{
							bdeactivate = true;
							break;
						}
					}
				}

				if (bdeactivate)
				{
					elem.push_back(pair<int, double>(el.GetID(), elemVal));
					deactiveElems++;
				}
			}
		}
	}
	if (deactiveElems == 0) return true;

	// sort the list
	std::sort(elem.begin(), elem.end(), [](pair<int, double>& e1, pair<int, double>& e2) {
		return e1.second > e2.second;
	});

	// deactivate the elements
	int elems = elem.size();
	if (elems > m_maxElems) elems = m_maxElems;
	for (int i = 0; i < elems; ++i)
	{
		pair<int, double>& e = elem[i];
		FEElement* pe = mesh.FindElementFromID(e.first); assert(pe);
		pe->setInactive();
	}

	// if any nodes were orphaned, we need to deactivate them
	int NN = mesh.Nodes();
	vector<int> tag(NN, 0);
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		int NE = dom.Elements();
		for (int j = 0; j < NE; ++j)
		{
			FEElement& el = dom.ElementRef(j);
			if (el.isActive())
			{
				int neln = el.Nodes();
				for (int n = 0; n < neln; ++n) tag[el.m_node[n]] = 1;
			}
		}
	}

	for (int i = 0; i < NN; ++i)
	{
		FENode& node = mesh.Node(i);
		if (tag[i] == 0)
		{
			vector<int>& id = node.m_ID;
			for (int j = 0; j < id.size(); ++j)
			{
				int n = id[j];
				if (n >= 0) id[j] = -n - 2;
			}
		}
	}

	return (elems == 0);
}
