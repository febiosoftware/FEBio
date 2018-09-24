#include "stdafx.h"
#include "FEModelData.h"
#include "FEModel.h"

BEGIN_PARAMETER_LIST(FEModelData, FECoreBase)
	ADD_PARAMETER(m_data, "value");
END_PARAMETER_LIST();

FEModelData::FEModelData(FEModel* fem, FELogElemData* eval, vector<int>& item) : FECoreBase(FEMODELDATA_ID)
{
	m_fem = fem;
	m_eval = eval;
	m_item = item;
	m_data.resize(item.size(), 0.0);
}

FEModel* FEModelData::GetFEModel()
{
	return m_fem;
}

void FEModelData::Update()
{
	// build the reference table
	if (m_ELT.empty()) BuildELT();

	// get the number of items
	int N = (int)m_item.size();

	if (m_eval == 0)
	{
		m_data.assign(N, 0.0);
		return;
	}

	FEMesh& m = m_fem->GetMesh();

	m_data.resize(N);
	for (int i=0; i<N; ++i)
	{
		int item = m_item[i];
		int index = item - m_offset;
		if ((index >= 0) && (index < m_ELT.size()))
		{
			ELEMREF& e = m_ELT[index];
			assert((e.ndom != -1) && (e.nid != -1));
			FEElement* pe = &m.Domain(e.ndom).ElementRef(e.nid); assert(pe);
			assert(pe->GetID() == item);

			// get the element value
			m_data[i] = m_eval->value(*pe);
		}
	}
}

//-----------------------------------------------------------------------------
void FEModelData::BuildELT()
{
	m_ELT.clear();
	FEMesh& m = m_fem->GetMesh();

	// find the min, max ID
	int minID = -1, maxID = 0;
	for (int i = 0; i<m.Domains(); ++i)
	{
		FEDomain& dom = m.Domain(i);
		int NE = dom.Elements();
		for (int j = 0; j<NE; ++j)
		{
			FEElement& el = dom.ElementRef(j);
			int id = el.GetID(); assert(id > 0);

			if ((minID < 0) || (id < minID)) minID = id;
			if (id > maxID) maxID = id;
		}
	}

	// allocate lookup table
	int nsize = maxID - minID + 1;
	m_ELT.resize(nsize);
	for (int i = 0; i<nsize; ++i)
	{
		m_ELT[i].ndom = -1;
		m_ELT[i].nid = -1;
	}

	// build lookup table
	m_offset = minID;
	for (int i = 0; i<m.Domains(); ++i)
	{
		FEDomain& d = m.Domain(i);
		int ne = d.Elements();
		for (int j = 0; j<ne; ++j)
		{
			FEElement& el = d.ElementRef(j);
			int id = el.GetID() - minID;
			m_ELT[id].ndom = i;
			m_ELT[id].nid = j;
		}
	}
}
