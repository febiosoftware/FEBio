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
#include "FENodeElemList.h"
#include "FESurface.h"
#include "FEMesh.h"
#include "FEDomain.h"
#include "FEShellDomain.h"
#include "DumpStream.h"

//-----------------------------------------------------------------------------
int FENodeElemList::MaxValence()
{
	int nmax = 0;
	int N = (int) m_nval.size();
	for (int i=0; i<N; ++i) if (m_nval[i] > nmax) nmax = m_nval[i];
	return nmax;
}

//-----------------------------------------------------------------------------
//! This function builds the node-element list for a surface

void FENodeElemList::Create(const FESurface& s)
{
	int i, j, n;

	// get the number of nodes
	int nn = s.Nodes();

	// get the number of elements
	int ne = s.Elements();

	// create nodal valence array
	m_nval.assign(nn, 0);
	m_pn.resize(nn);

	// fill valence table
	int nsize = 0;
	for (i=0; i<ne; ++i)
	{
		const FESurfaceElement& el = s.Element(i);

		for (j=0; j<el.Nodes(); ++j)
		{
			n = el.m_lnode[j];
			m_nval[n]++;
			nsize++;
		}
	}

	// create the element reference array
	m_eref.resize(nsize);
	m_iref.resize(nsize);

	// set eref pointers
	m_pn[0] = 0;
	for (i=1; i<nn; ++i)
	{
		m_pn[i] = m_pn[i-1] + m_nval[i-1];
	}

	// reset valence pointers
	for (i=0; i<nn; ++i) m_nval[i] = 0;

	// fill eref table
	for (i=0; i<ne; ++i)
	{
		const FESurfaceElement& el = s.Element(i);

		for (j=0; j<el.Nodes(); ++j)
		{
			n = el.m_lnode[j];
			m_eref[m_pn[n] + m_nval[n]] = const_cast<FESurfaceElement*>(&el);
			m_iref[m_pn[n] + m_nval[n]] = i;
			m_nval[n]++;
		}
	}
}

//-----------------------------------------------------------------------------
//! This function builds the node-element list for a mesh

void FENodeElemList::Create(FEMesh& mesh)
{
	int i, j, n, nd;

	// get the number of nodes
	int NN = mesh.Nodes();

	// create nodal valence array
	m_nval.assign(NN, 0);
	m_pn.resize(NN);

	// fill valence table
	int nsize = 0;
	for (nd=0; nd<mesh.Domains(); ++nd)
	{
		FEDomain& d = mesh.Domain(nd);
		for (i=0; i<d.Elements(); ++i)
		{
			FEElement& el = d.ElementRef(i);
			for (j=0; j<el.Nodes(); ++j)
			{
				n = el.m_node[j];
				m_nval[n]++;
				nsize++;
			}
		}
	}

	// create the element reference array
	m_eref.resize(nsize);
	m_iref.resize(nsize);

	// set eref pointers
	m_pn[0] = 0;
	for (i=1; i<NN; ++i)
	{
		m_pn[i] = m_pn[i-1] + m_nval[i-1];
	}

	// reset valence pointers
	for (i=0; i<NN; ++i) m_nval[i] = 0;

	// fill eref table
    // Prioritize shell domains over other domains.
    // This is needed when shells are connected to solids
    // and contact interfaces need to use the shell properties
    // for auto-penalty calculation.
	int nindex = 0;
	for (nd=0; nd<mesh.Domains(); ++nd)
	{
		FEDomain& d = mesh.Domain(nd);
        if (d.Class() == FE_DOMAIN_SHELL) {
            for (i=0; i<d.Elements(); ++i, ++nindex)
            {
                FEElement& el = d.ElementRef(i);
                for (j=0; j<el.Nodes(); ++j)
                {
                    n = el.m_node[j];
                    m_eref[m_pn[n] + m_nval[n]] = &el;
					m_iref[m_pn[n] + m_nval[n]] = nindex;
					m_nval[n]++;
                }
            }
        }
	}
    for (nd=0; nd<mesh.Domains(); ++nd)
    {
        FEDomain& d = mesh.Domain(nd);
        if (d.Class() != FE_DOMAIN_SHELL) {
            for (i=0; i<d.Elements(); ++i, ++nindex)
            {
                FEElement& el = d.ElementRef(i);
                for (j=0; j<el.Nodes(); ++j)
                {
                    n = el.m_node[j];
                    m_eref[m_pn[n] + m_nval[n]] = &el;
					m_iref[m_pn[n] + m_nval[n]] = nindex;
					m_nval[n]++;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
//! This function builds the node-element list for a domain

void FENodeElemList::Create(FEDomain& dom)
{
	int i, j, n;

	// get the mesh
	FEMesh& mesh = *dom.GetMesh();

	// get the number of nodes
	int NN = mesh.Nodes();

	// create nodal valence array
	m_nval.assign(NN, 0);
	m_pn.resize(NN);

	// fill valence table
	int nsize = 0;
	for (i=0; i<dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		for (j=0; j<el.Nodes(); ++j)
		{
			n = el.m_node[j];
			m_nval[n]++;
			nsize++;
		}
	}

	// create the element reference array
	m_eref.resize(nsize);
	m_iref.resize(nsize);

	// set eref pointers
	m_pn[0] = 0;
	for (i=1; i<NN; ++i)
	{
		m_pn[i] = m_pn[i-1] + m_nval[i-1];
	}

	// reset valence pointers
	for (i=0; i<NN; ++i) m_nval[i] = 0;

	// fill eref table
	for (i=0; i<dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		for (j=0; j<el.Nodes(); ++j)
		{
			n = el.m_node[j];
			m_eref[m_pn[n] + m_nval[n]] = &el;
			m_iref[m_pn[n] + m_nval[n]] = i;
			m_nval[n]++;
		}
	}
}

//-----------------------------------------------------------------------------
void FENodeElemList::Clear()
{
	m_nval.clear();
	m_eref.clear();
	m_iref.clear();
	m_pn.clear();
}

//-----------------------------------------------------------------------------
//! Save data to dump file

void FENodeElemList::Serialize(DumpStream& ar)
{
	ar & m_nval & m_iref & m_pn;
}

//-----------------------------------------------------------------------------
void FENodeElemTree::Create(FESurface* ps, int k)
{
	int NN = ps->Nodes();
	int NE = ps->Elements();

	// temporary arrays
	vector< vector<int> > nel;
	vector<int> tag;
	nel.resize(NN);
	tag.assign(NE, -1);

	// build the first level
	for (int i=0; i<NE; ++i)
	{
		FESurfaceElement* pe = &ps->Element(i);
		int ne = pe->Nodes();
		for (int j=0; j<ne; ++j) nel[pe->m_lnode[j]].push_back(i);
	}

	// build the other levels
	for (int l=0; l<k; ++l)
	{
		vector<int> ns(NN);
		for (int i=0; i<NN; ++i) ns[i] = (int) nel[i].size();

		for (int i=0; i<NN; ++i)
		{
			int ntag = l*NN + i;
			vector<int>& NI = nel[i];
			int ni = ns[i];
			for (int j=0; j<ni; ++j) tag[NI[j]] = ntag;

			for (int j=0; j<ni; ++j)
			{
				FESurfaceElement& e = ps->Element(NI[j]);
				int ne = e.Nodes();
				for (int n=0; n<ne; ++n)
				{
					if (e.m_lnode[n] != i)
					{
						vector<int>& NJ = nel[e.m_lnode[n]];
						int nj = ns[e.m_lnode[n]];
						for (int m=0; m<nj; ++m)
						{
							if (tag[NJ[m]] < ntag) 
							{
								NI.push_back(NJ[m]);
								tag[NJ[m]] = ntag;
							}
						}
					}
				}
			}
		}
	}

	// assign the element pointers
	m_nel.resize(NN);
	for (int i=0; i<NN; ++i)
	{
		vector<int>& NI = nel[i];
		sort(NI.begin(), NI.end());
		int ni = (int)NI.size();
		m_nel[i].resize(ni);
		for (int j=0; j<ni; ++j) m_nel[i][j] = &ps->Element(NI[j]);
	}
}
