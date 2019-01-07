#include "stdafx.h"
#include "NodeDataRecord.h"
#include "FEAnalysis.h"
#include "FECoreKernel.h"
#include "FEModel.h"

REGISTER_SUPER_CLASS(FENodeLogData, FENODELOGDATA_ID);

//-----------------------------------------------------------------------------
FENodeLogData::FENodeLogData(FEModel* fem) : FECoreBase(fem, FENODELOGDATA_ID) {}

//-----------------------------------------------------------------------------
FENodeLogData::~FENodeLogData() {}

//-----------------------------------------------------------------------------
NodeDataRecord::NodeDataRecord(FEModel* pfem, const char* szfile) : DataRecord(pfem, szfile, FE_DATA_NODE) {}

//-----------------------------------------------------------------------------
int NodeDataRecord::Size() const { return (int)m_Data.size(); }

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
		FENodeLogData* pdata = fecore_new<FENodeLogData>(sz, m_pfem);
		if (pdata) m_Data.push_back(pdata);
		else 
		{
			// see if this refers to a DOF of the model
			int ndof = m_pfem->GetDOFIndex(sz);
			if (ndof >= 0)
			{
				// Add an output for a nodal variable
				pdata = new FENodeVarData(m_pfem, ndof);
				m_Data.push_back(pdata);
			}
			else throw UnknownDataField(sz);
		}
		sz = ch;
	}
	while (ch);
}

//-----------------------------------------------------------------------------
double NodeDataRecord::Evaluate(int item, int ndata)
{
	FEMesh& mesh = m_pfem->GetMesh();
	int nnode = item - 1;
	assert((nnode>=0)&&(nnode<mesh.Nodes()));
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
// This sets the item list based on a node set.
// Note that node sets store the nodes in a zero-based list. However, we need
// a one-base list here.
void NodeDataRecord::SetItemList(FENodeSet* pns)
{
	int n = pns->size();
	assert(n);
	m_item.resize(n);
	for (int i=0; i<n; ++i) m_item[i] = (*pns)[i] + 1;
}

//-----------------------------------------------------------------------------
FENodeVarData::FENodeVarData(FEModel* pfem, int ndof) : FENodeLogData(pfem), m_ndof(ndof) {}

//-----------------------------------------------------------------------------
double FENodeVarData::value(int node)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	return mesh.Node(node).get(m_ndof);
}
