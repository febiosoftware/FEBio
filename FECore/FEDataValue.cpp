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
#include "FEDataValue.h"
#include "FELogNodeData.h"
#include "FELogElemData.h"
#include "FEModel.h"
#include "FEMesh.h"

FEDataValue::FEDataValue()
{

}

bool FEDataValue::IsValid() const
{
	return (m_logData != nullptr);
}

void FEDataValue::SetLogData(FELogData* logData)
{
	m_logData = logData;
}

bool FEDataValue::GetValues(const FEItemList* itemList, std::vector<double>& val)
{
	if (m_logData == nullptr) return false;
	if (itemList == nullptr) return false;

	FELogNodeData* nodeData = dynamic_cast<FELogNodeData*>(m_logData);
	if (nodeData)
	{
		const FENodeSet* nset = dynamic_cast<const FENodeSet*>(itemList);
		if (nset == nullptr) return false;

		FEModel* fem = nodeData->GetFEModel();
		FEMesh& mesh = fem->GetMesh();

		int n = nset->Size();
		val.resize(n);
		for (int i = 0; i < n; ++i)
		{
			int nid = (*nset)[i];
			FENode* node = mesh.FindNodeFromID(nid);
			if (node == nullptr) return false;
			val[i] = nodeData->value(*node);
		}
		return true;
	}

	FELogElemData* elemData = dynamic_cast<FELogElemData*>(m_logData);
	if (elemData)
	{
		const FEElementSet* eset = dynamic_cast<const FEElementSet*>(itemList);
		if (eset == nullptr) return false;

		FEModel* fem = elemData->GetFEModel();
		FEMesh& mesh = fem->GetMesh();

		int n = eset->Elements();
		val.resize(n);
		for (int i = 0; i < n; ++i)
		{
			int eid = (*eset)[i];
			FEElement* elem = mesh.FindElementFromID(eid);
			if (elem == nullptr) return false;
			val[i] = elemData->value(*elem);
		}
		return true;
	}

	return false;
}
