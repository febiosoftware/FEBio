#include "stdafx.h"
#include "FEMeshAdaptor.h"
#include "FEModel.h"
#include "FEDomain.h"

REGISTER_SUPER_CLASS(FEMeshAdaptor, FEMESHADAPTOR_ID);

FEMeshAdaptor::FEMeshAdaptor(FEModel* fem) : FECoreBase(fem)
{

}

REGISTER_SUPER_CLASS(FEMeshAdaptorCriterion, FEMESHADAPTORCRITERION_ID);

BEGIN_FECORE_CLASS(FEMeshAdaptorCriterion, FECoreBase)
	ADD_PARAMETER(m_sortList, "sort");
	ADD_PARAMETER(m_maxelem, "max_elem");
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
std::vector<int> FEMeshAdaptorCriterion::GetElementList()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	int nselected = 0;
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
				double elemVal = 0;
				bool bselect = Check(el, elemVal);
				if (bselect)
				{
					elem.push_back(pair<int, double>(el.GetID(), elemVal));
					nselected++;
				}
			}
		}
	}

	std::vector<int> selectedElement;
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
			selectedElement[i] = elem[i].first;
		}
	}

	return selectedElement;
}

bool FEMeshAdaptorCriterion::Check(FEElement& el, double& elemVal)
{
	return false;
}
