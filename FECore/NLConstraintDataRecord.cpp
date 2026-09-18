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
#include "NLConstraintDataRecord.h"
#include "FECoreKernel.h"
#include "FEModel.h"

void NLConstraintDataRecord::SetData(const char* szexpr)
{
	std::vector<DataRecordItem> data = ProcessDataString(szexpr);
	if (data.empty()) throw UnknownDataField(szexpr);

	m_Data.clear();
	m_data = szexpr;
	for (int i=0; i<data.size(); ++i)
	{
		FELogNLConstraintData* pdata = fecore_new<FELogNLConstraintData>(data[i].name.c_str(), GetFEModel());
		if (pdata == nullptr) throw UnknownDataField(data[i].name);
		m_Data.push_back(pdata);

		if (!ApplyDataRecordItem(*pdata, data[i]))
			throw UnknownDataField(data[i].name);
	}
}

NLConstraintDataRecord::NLConstraintDataRecord(FEModel* pfem) : DataRecord(pfem, FE_DATA_NLC) {}

int NLConstraintDataRecord::Size() const { return (int)m_Data.size(); }

double NLConstraintDataRecord::Evaluate(int item, int ndata)
{
    FEModel* fem = GetFEModel();
    int nc = item - 1;
    if ((nc < 0) || (nc >= fem->NonlinearConstraints())) return 0;
    
	FENLConstraint& nlc = *fem->NonlinearConstraint(nc);
	return m_Data[ndata]->value(nlc);
}

void NLConstraintDataRecord::SelectAllItems()
{
    FEModel* fem = GetFEModel();
    int n = fem->NonlinearConstraints();
	m_item.resize(n);
	for (int i = 0; i<n; ++i) m_item[i] = i + 1;
}
