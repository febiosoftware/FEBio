#include "stdafx.h"
#include "ElementDataRecord.h"
#include "FECoreKernel.h"
#include "FEModel.h"

//-----------------------------------------------------------------------------
ElementDataRecord::ElementDataRecord(FEModel* pfem, const char* szfile) : DataRecord(pfem, szfile, FE_DATA_ELEM)
{
	m_offset = 0;
}

//-----------------------------------------------------------------------------
void ElementDataRecord::Parse(const char *szexpr)
{
	char szcopy[MAX_STRING] = {0};
	strcpy(szcopy, szexpr);
	char* sz = szcopy, *ch;
	m_Data.clear();
	strcpy(m_szdata, szexpr);
	do
	{
		ch = strchr(sz, ';');
		if (ch) *ch++ = 0;
		FELogElemData* pdata = fecore_new<FELogElemData>(FEELEMLOGDATA_ID, sz, m_pfem);
		if (pdata) m_Data.push_back(pdata);
		else throw UnknownDataField(sz);
		sz = ch;
	}
	while (ch);
}

//-----------------------------------------------------------------------------
double ElementDataRecord::Evaluate(int item, int ndata)
{
	// make sure we have an ELT
	if (m_ELT.empty()) BuildELT();

	// find the element
	FEMesh& mesh = m_pfem->GetMesh();
	int index = item - m_offset;
	if ((index >= 0) && (index < m_ELT.size()))
	{
		ELEMREF& e = m_ELT[index];
		assert((e.ndom != -1) && (e.nid != -1));
		FEElement* pe = &mesh.Domain(e.ndom).ElementRef(e.nid); assert(pe);
		assert(pe->GetID() == item);

		// get the element value
		return m_Data[ndata]->value(*pe);
	}
	else return 0.0;
}

//-----------------------------------------------------------------------------
void ElementDataRecord::BuildELT()
{
	m_ELT.clear();
	FEMesh& m = m_pfem->GetMesh();

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
	for (int i=0; i<nsize; ++i) 
	{
		m_ELT[i].ndom = -1;
		m_ELT[i].nid  = -1;
	}

	// build lookup table
	m_offset = minID;
	for (int i=0; i<m.Domains(); ++i)
	{
		FEDomain& d = m.Domain(i);
		int ne = d.Elements();
		for (int j=0; j<ne; ++j)
		{
			FEElement& el = d.ElementRef(j);
			int id = el.GetID() - minID;
			m_ELT[id].ndom = i;
			m_ELT[id].nid  = j;
		}
	}
}

//-----------------------------------------------------------------------------
void ElementDataRecord::SelectAllItems()
{
	FEMesh& m = m_pfem->GetMesh();
	int n = m.Elements();
	m_item.resize(n);
	n = 0;
	for (int i=0; i<m.Domains(); ++i)
	{
		FEDomain& dom = m.Domain(i);
		int NE = dom.Elements();
		for (int j=0; j<NE; ++j, n++)
		{
			FEElement& el = dom.ElementRef(j);
			m_item[n] = el.GetID();
		}
	}
}

//-----------------------------------------------------------------------------
// This sets the item list based on a element set.
void ElementDataRecord::SetItemList(FEElementSet* pg)
{
	int n = pg->size();
	assert(n);
	m_item.resize(n);
	for (int i=0; i<n; ++i) m_item[i] = (*pg)[i];
}
