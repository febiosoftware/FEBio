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
#include "FaceDataRecord.h"
#include "FECoreKernel.h"
#include "FEModel.h"
#include "FESurface.h"

//-----------------------------------------------------------------------------
FELogFaceData::FELogFaceData(FEModel* fem) : FELogData(fem) {}

//-----------------------------------------------------------------------------
FELogFaceData::~FELogFaceData() {}

//-----------------------------------------------------------------------------
FaceDataRecord::FaceDataRecord(FEModel* pfem) : DataRecord(pfem, FE_DATA_FACE) 
{
	m_surface = nullptr;
}

//-----------------------------------------------------------------------------
int FaceDataRecord::Size() const { return (int)m_Data.size(); }

//-----------------------------------------------------------------------------
void FaceDataRecord::SetData(const char* szexpr)
{
	char szcopy[MAX_STRING] = { 0 };
	strcpy(szcopy, szexpr);
	char* sz = szcopy, *ch;
	m_Data.clear();
	strcpy(m_szdata, szexpr);
	do
	{
		ch = strchr(sz, ';');
		if (ch) *ch++ = 0;
		FELogFaceData* pdata = fecore_new<FELogFaceData>(sz, GetFEModel());
		if (pdata) m_Data.push_back(pdata);
		else throw UnknownDataField(sz);
		sz = ch;
	}
	while (ch);
}

//-----------------------------------------------------------------------------
bool FaceDataRecord::Initialize()
{
	return (m_item.empty() == false);
}

//-----------------------------------------------------------------------------
double FaceDataRecord::Evaluate(int item, int ndata)
{
	int nface = item - 1;
	return m_Data[ndata]->value(m_surface->Element(nface));
}

//-----------------------------------------------------------------------------
void FaceDataRecord::SelectAllItems()
{
	assert(false);
}

void FaceDataRecord::SetItemList(FEItemList* itemList, const std::vector<int>& selection)
{
	FEFacetSet* facetSet = dynamic_cast<FEFacetSet*>(itemList); assert(facetSet);
	m_surface = facetSet->GetSurface(); assert(m_surface);
	int n = m_surface->Elements();
	if (selection.empty())
	{
		m_item.resize(n);
		for (int i = 0; i < n; ++i) m_item[i] = i + 1;
	}
	else
	{
		m_item = selection;
	}
}
