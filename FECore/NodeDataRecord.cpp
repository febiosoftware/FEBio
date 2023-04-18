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
#include "NodeDataRecord.h"
#include "FEAnalysis.h"
#include "FECoreKernel.h"
#include "FEModel.h"

//-----------------------------------------------------------------------------
FELogNodeData::FELogNodeData(FEModel* fem) : FELogData(fem) {}

//-----------------------------------------------------------------------------
FELogNodeData::~FELogNodeData() {}

//-----------------------------------------------------------------------------
NodeDataRecord::NodeDataRecord(FEModel* pfem) : DataRecord(pfem, FE_DATA_NODE) {}

//-----------------------------------------------------------------------------
int NodeDataRecord::Size() const { return (int)m_Data.size(); }

//-----------------------------------------------------------------------------
void NodeDataRecord::SetData(const char* szexpr)
{
	char szcopy[MAX_STRING] = {0};
	strcpy(szcopy, szexpr);
	char* sz = szcopy, *ch;
	m_Data.clear();
	strcpy(m_szdata, szexpr);
	FEModel* fem = GetFEModel();
	do
	{
		ch = strchr(sz, ';');
		if (ch) *ch++ = 0;
		FELogNodeData* pdata = fecore_new<FELogNodeData>(sz, fem);
		if (pdata) m_Data.push_back(pdata);
		else 
		{
			// see if this refers to a DOF of the model
			int ndof = fem->GetDOFIndex(sz);
			if (ndof >= 0)
			{
				// Add an output for a nodal variable
				pdata = new FENodeVarData(fem, ndof);
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
	FEMesh& mesh = GetFEModel()->GetMesh();
	int nnode = item - 1;
	assert((nnode>=0)&&(nnode<mesh.Nodes()));
	if ((nnode < 0) || (nnode >= mesh.Nodes())) return 0;

	FENode& node = mesh.Node(nnode);
	return m_Data[ndata]->value(node);
}

//-----------------------------------------------------------------------------
void NodeDataRecord::SelectAllItems()
{
	int n = GetFEModel()->GetMesh().Nodes();
	m_item.resize(n);
	for (int i=0; i<n; ++i) m_item[i] = i+1;
}

//-----------------------------------------------------------------------------
void NodeDataRecord::SetItemList(FEItemList* items, const std::vector<int>& selection)
{
	// TODO: We don't support using a selection of a node set yet. 
	assert(selection.empty());
	FENodeSet* pns = dynamic_cast<FENodeSet*>(items); assert(pns);
	int n = pns->Size();
	m_item.resize(n);
	for (int i = 0; i < n; ++i) m_item[i] = (*pns)[i] + 1;
}

//-----------------------------------------------------------------------------
FENodeVarData::FENodeVarData(FEModel* pfem, int ndof) : FELogNodeData(pfem), m_ndof(ndof) {}

//-----------------------------------------------------------------------------
double FENodeVarData::value(const FENode& node)
{
	return node.get(m_ndof);
}
