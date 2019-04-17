/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FEElementSet.h"
#include "FEElement.h"
#include "FEMesh.h"
#include "FEDomain.h"
#include "DumpStream.h"
#include "FEModel.h"

//-----------------------------------------------------------------------------
FEElementSet::FEElementSet(FEModel* fem) : FEItemList(fem)
{
	m_minID = -1;
	m_maxID = -1;
}

//-----------------------------------------------------------------------------
void FEElementSet::Create(const std::vector<int>& elemList)
{
	m_dom.Clear();
	m_Elem = elemList;
	BuildLUT();
}

//-----------------------------------------------------------------------------
void FEElementSet::Create(FEDomain* dom)
{
	m_dom.Clear();
	m_dom.AddDomain(dom);

	int NE = dom->Elements();
	m_Elem.resize(NE, -1);
	for (int i = 0; i < NE; ++i)
	{
		FEElement& el = dom->ElementRef(i);
		m_Elem[i] = el.GetID();
	}

	BuildLUT();
}

//-----------------------------------------------------------------------------
void FEElementSet::Create(FEDomainList& domList)
{
	int NT = 0;
	m_dom.Clear();
	for (int n = 0; n < domList.Domains(); ++n)
	{
		FEDomain* dom = domList.GetDomain(n);
		m_dom.AddDomain(dom);
		NT += dom->Elements();
	}

	m_Elem.resize(NT, -1);
	NT = 0;
	for (int n = 0; n < domList.Domains(); ++n)
	{
		FEDomain* dom = domList.GetDomain(n);
		int NE = dom->Elements();
		for (int i = 0; i < NE; ++i)
		{
			FEElement& el = dom->ElementRef(i);
			m_Elem[NT + i] = el.GetID();
		}
		NT += NE;
	}

	BuildLUT();
}

//-----------------------------------------------------------------------------
void FEElementSet::BuildLUT()
{
	FEMesh* mesh = GetMesh();
	int N = (int)m_Elem.size();
	m_minID = m_maxID = -1;
	for (int i = 0; i < N; ++i)
	{
		FEElement* pe = mesh->FindElementFromID(m_Elem[i]);
		int id = pe->GetID();

		if ((id < m_minID) || (m_minID == -1)) m_minID = id;
		if ((id > m_maxID) || (m_maxID == -1)) m_maxID = id;
	}

	int lutSize = m_maxID - m_minID + 1;
	m_LUT.resize(lutSize, -1);
	for (int i = 0; i < N; ++i)
	{
		FEElement* pe = mesh->FindElementFromID(m_Elem[i]);
		int id = pe->GetID() - m_minID;
		m_LUT[id] = i;
	}
}

//-----------------------------------------------------------------------------
FEElement& FEElementSet::Element(int i)
{
	FEMesh* mesh = GetMesh();
	return *mesh->FindElementFromID(m_Elem[i]);
}

//-----------------------------------------------------------------------------
void FEElementSet::Serialize(DumpStream& ar)
{
	FEItemList::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_Elem;
	ar & m_LUT;
	ar & m_minID & m_maxID;
}
