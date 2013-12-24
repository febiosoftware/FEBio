#include "stdafx.h"
#include "NodeDataRecord.h"
#include "FEAnalysis.h"
#include "FECoreKernel.h"
#include "FEModel.h"

//-----------------------------------------------------------------------------
void NodeDataRecord::Parse(const char* szexpr)
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
		FENodeLogData* pdata = fecore_new<FENodeLogData>(FENODELOGDATA_ID, sz, m_pfem);
		if (pdata) m_Data.push_back(pdata);
		else throw UnknownDataField(sz);
		sz = ch;
	}
	while (ch);
}

//-----------------------------------------------------------------------------
double NodeDataRecord::Evaluate(int item, int ndata)
{
	FEMesh& mesh = m_pfem->GetMesh();
	int nnode = item - 1;
	if ((nnode < 0) || (nnode >= mesh.Nodes())) return 0;
	return m_Data[ndata]->value(nnode);
}

//-----------------------------------------------------------------------------
void NodeDataRecord::SelectAllItems()
{
	int n = m_pfem->GetMesh().Nodes();
	m_item.resize(n);
	for (int i=0; i<n; ++i) m_item[i] = i+1;
}

//-----------------------------------------------------------------------------
void NodeDataRecord::SetItemList(FENodeSet* pns)
{
	int n = pns->size();
	assert(n);
	m_item.resize(n);
	for (int i=0; i<n; ++i) m_item[i] = (*pns)[i];
}
