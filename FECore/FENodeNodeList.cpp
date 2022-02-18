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
#include "FENodeNodeList.h"
#include "FENodeElemList.h"
#include "FEMesh.h"
#include "FEDomain.h"
#include <stdlib.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FENodeNodeList::FENodeNodeList()
{

}

FENodeNodeList::~FENodeNodeList()
{

}

FENodeNodeList* FENodeNodeList::m_pthis = 0;

//////////////////////////////////////////////////////////////////////
// FENodeNodeList
//////////////////////////////////////////////////////////////////////

void FENodeNodeList::Create(FEMesh& mesh)
{
	int i, j, k, n, m;

	// get the nr of nodes
	int NN = mesh.Nodes();

	// create the node-element list
	FENodeElemList EL; 
	EL.Create(mesh);

	// create the nodal tag array
	vector<int> tag; tag.assign(NN, 0);

	// calculate nodal valences
	m_nval.assign(NN, 0);
	m_pn.resize(NN);

	int nsize = 0;
	int* en;
	vector<int> buf(NN);
	int nb;
	for (i=0; i<NN; ++i)
	{
		nb = 0;
		n = EL.Valence(i);
		FEElement** pe = EL.ElementList(i);
		for (j=0; j<n; ++j)
		{
			FEElement* pel = pe[j];
			m = pel->Nodes();
			en = &pel->m_node[0];
			for (k=0; k<m; ++k)
				if ((en[k] != i) && (tag[ en[k] ] == 0))
				{
					++m_nval[i];
					++tag[en[k]];
					buf[nb++] = en[k];
					++nsize;
				}
		}

		// clear the tag array
		for (j=0; j<nb; ++j) tag[ buf[j] ] = 0;
		nb = 0;
	}

	// create the node reference array
	m_nref.resize(nsize);

	// set nref pointers
	m_pn[0] = 0;
	for (i=1; i<NN; ++i)
	{
		m_pn[i] = m_pn[i-1] + m_nval[i-1];
	}

	// reset valence pointers
	for (i=0; i<NN; ++i) m_nval[i] = 0;

	// fill the nref pointers
	for (i=0; i<NN; ++i)
	{
		nb = 0;
		n = EL.Valence(i);
		FEElement** pe = EL.ElementList(i);
		for (j=0; j<n; ++j)
		{
			FEElement* pel = pe[j];
			m = pel->Nodes();
			en = &pel->m_node[0];
			for (k=0; k<m; ++k)
				if ((en[k] != i) && (tag[ en[k] ] == 0))
				{
					m_nref[m_pn[i] + m_nval[i]] = en[k];	

					++tag[en[k]];
					++m_nval[i];
					buf[nb++] = en[k];
					++nsize;
				}
		}

		// clear the tag array
		for (j=0; j<nb; ++j) tag[ buf[j] ] = 0;
		nb = 0;
	}
}

//-----------------------------------------------------------------------------
void FENodeNodeList::Create(FEDomain& dom)
{
	int i, j, k, n, m;

	// get the mesh
	FEMesh& mesh = *dom.GetMesh();

	// get the nr of nodes
	int NN = mesh.Nodes();

	// create the node-element list
	FENodeElemList EL; 
	EL.Create(dom);

	// create the nodal tag array
	vector<int> tag; tag.assign(NN, 0);

	// calculate nodal valences
	m_nval.assign(NN, 0);
	m_pn.resize(NN);

	int nsize = 0;
	int* en;
	vector<int> buf(NN);
	int nb;
	for (i=0; i<NN; ++i)
	{
		nb = 0;
		n = EL.Valence(i);
		FEElement** pe = EL.ElementList(i);
		for (j=0; j<n; ++j)
		{
			FEElement* pel = pe[j];
			m = pel->Nodes();
			en = &pel->m_node[0];
			for (k=0; k<m; ++k)
				if ((en[k] != i) && (tag[ en[k] ] == 0))
				{
					++m_nval[i];
					++tag[en[k]];
					buf[nb++] = en[k];
					++nsize;
				}
		}

		// clear the tag array
		for (j=0; j<nb; ++j) tag[ buf[j] ] = 0;
		nb = 0;
	}

	// create the node reference array
	m_nref.resize(nsize);

	// set nref pointers
	m_pn[0] = 0;
	for (i=1; i<NN; ++i)
	{
		m_pn[i] = m_pn[i-1] + m_nval[i-1];
	}

	// reset valence pointers
	for (i=0; i<NN; ++i) m_nval[i] = 0;

	// fill the nref pointers
	for (i=0; i<NN; ++i)
	{
		nb = 0;
		n = EL.Valence(i);
		FEElement** pe = EL.ElementList(i);
		for (j=0; j<n; ++j)
		{
			FEElement* pel = pe[j];
			m = pel->Nodes();
			en = &pel->m_node[0];
			for (k=0; k<m; ++k)
				if ((en[k] != i) && (tag[ en[k] ] == 0))
				{
					m_nref[m_pn[i] + m_nval[i]] = en[k];	

					++tag[en[k]];
					++m_nval[i];
					buf[nb++] = en[k];
					++nsize;
				}
		}

		// clear the tag array
		for (j=0; j<nb; ++j) tag[ buf[j] ] = 0;
		nb = 0;
	}
}


///////////////////////////////////////////////////////////////////////////////

int FENodeNodeList::compare(const void* e1, const void* e2)
{
	int n1 = *((int*) e1);
	int n2 = *((int*) e2);

	FENodeNodeList& L = *m_pthis;

	return (L.Valence(n1) - L.Valence(n2));
}

void FENodeNodeList::Sort()
{
	int n, *pn;
	m_pthis = this;
	for (int i=0; i<Size(); ++i)
	{
		n = Valence(i);
		pn = NodeList(i);
		qsort(pn, n, sizeof(int), compare);
	}
}
