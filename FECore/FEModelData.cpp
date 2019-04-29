/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "FEModelData.h"
#include "FEModel.h"
#include "FEDomain.h"

REGISTER_SUPER_CLASS(FEModelData, FEMODELDATA_ID);

BEGIN_FECORE_CLASS(FEModelData, FECoreBase)
	ADD_PARAMETER(m_data, "value");
END_FECORE_CLASS();

FEModelData::FEModelData(FEModel* fem, FELogElemData* eval, vector<int>& item) : FECoreBase(fem)
{
	m_eval = eval;
	m_item = item;
	m_data.resize(item.size(), 0.0);
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

	FEMesh& m = GetFEModel()->GetMesh();

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
	FEMesh& m = GetFEModel()->GetMesh();

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
