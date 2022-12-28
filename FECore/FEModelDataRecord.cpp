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
#include "FEModelDataRecord.h"

FEModelLogData::FEModelLogData(FEModel* fem) : FELogData(fem) {}
FEModelLogData::~FEModelLogData() {}

FEModelDataRecord::FEModelDataRecord(FEModel* pfem) : DataRecord(pfem, FE_DATA_MODEL) {}
double FEModelDataRecord::Evaluate(int item, int ndata)
{
	assert(item == 0);
	return m_data[ndata]->value();
}

void FEModelDataRecord::ClearData()
{
	for (int i = 0; i < m_data.size(); ++i) delete m_data[i];
	m_data.clear();
}

void FEModelDataRecord::SetData(const char* szexpr)
{
	char szcopy[MAX_STRING] = { 0 };
	strcpy(szcopy, szexpr);
	char* sz = szcopy, * ch;
	ClearData();
	strcpy(m_szdata, szexpr);
	do
	{
		ch = strchr(sz, ';');
		if (ch) *ch++ = 0;
		FEModelLogData* pdata = fecore_new<FEModelLogData>(sz, GetFEModel());
		if (pdata) m_data.push_back(pdata);
		else throw UnknownDataField(sz);
		sz = ch;
	} while (ch);
}

void FEModelDataRecord::SelectAllItems()
{
	std::vector<int> items;
	items.push_back(0);
	SetItemList(items);
}

int FEModelDataRecord::Size() const
{
	return 1;
}
