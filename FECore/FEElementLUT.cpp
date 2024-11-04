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
#include "FEElementLUT.h"
#include "FEMesh.h"
#include "FEElement.h"
#include "FEDomain.h"

FEElementLUT::FEElementLUT(FEMesh& mesh)
{
	// get the ID ranges
	m_minID = -1;
	m_maxID = -1;
	int NDOM = mesh.Domains();
	for (int i = 0; i < NDOM; ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		int NE = dom.Elements();
		for (int j = 0; j < NE; ++j)
		{
			FEElement& el = dom.ElementRef(j);
			int eid = el.GetID();
			if ((eid < m_minID) || (m_minID == -1)) m_minID = eid;
			if ((eid > m_maxID) || (m_maxID == -1)) m_maxID = eid;
		}
	}

	// allocate size
	int nsize = m_maxID - m_minID + 1;
	m_elem.resize(nsize, (FEElement*)0);
	m_elid.resize(nsize, -1);

	// fill the table
	int index = 0;
	for (int i = 0; i < NDOM; ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		int NE = dom.Elements();
		for (int j = 0; j < NE; ++j, ++index)
		{
			FEElement& el = dom.ElementRef(j);
			int eid = el.GetID();
			m_elem[eid - m_minID] = &el;
			m_elid[eid - m_minID] = index;
		}
	}
}

// Find an element from its ID
FEElement* FEElementLUT::Find(int elemID) const
{
	if ((elemID < m_minID) || (elemID > m_maxID)) return nullptr;
	return m_elem[elemID - m_minID];
}

int FEElementLUT::FindIndex(int elemID) const
{
	if ((elemID < m_minID) || (elemID > m_maxID)) return -1;
	return m_elid[elemID - m_minID];
}
