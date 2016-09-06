#include "stdafx.h"
#include "ElementDataRecord.h"
#include "FECoreKernel.h"
#include "FEModel.h"

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
	assert((item >= 1) && (item <= mesh.Elements()));
	ELEMREF e = m_ELT[item-1];
	assert((e.ndom != -1) && (e.nid != -1));
	FEElement* pe = &mesh.Domain(e.ndom).ElementRef(e.nid); assert(pe);

	// get the element value
	return m_Data[ndata]->value(*pe);
}

//-----------------------------------------------------------------------------
void ElementDataRecord::BuildELT()
{
	m_ELT.clear();
	FEMesh& m = m_pfem->GetMesh();
	int NE = m.Elements();
	m_ELT.resize(NE);
	for (int i=0; i<NE; ++i) 
	{
		m_ELT[i].ndom = -1;
		m_ELT[i].nid  = -1;
	}

	int n = 0;
	for (int i=0; i<m.Domains(); ++i)
	{
		FEDomain& d = m.Domain(i);
		int ne = d.Elements();
		for (int j=0; j<ne; ++j, ++n)
		{
			FEElement& el = d.ElementRef(j);
			m_ELT[n].ndom = i;
			m_ELT[n].nid  = j;
		}
	}
}

//-----------------------------------------------------------------------------
void ElementDataRecord::SelectAllItems()
{
	FEMesh& m = m_pfem->GetMesh();
	int n = m.Elements();
	m_item.resize(n);
	for (int i=0; i<n; ++i) m_item[i] = i+1;
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
