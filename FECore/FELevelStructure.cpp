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
#include "FELevelStructure.h"
#include <stdlib.h>
#include <queue>
#include <stack>
#include "vector.h"
#include <assert.h>
using namespace std;

//-----------------------------------------------------------------------------

FELevelStructure::FELevelStructure()
{
	m_pNL = 0;
}

//-----------------------------------------------------------------------------

FELevelStructure::~FELevelStructure()
{

}

FELevelStructure* FELevelStructure::m_pthis = 0;

//-----------------------------------------------------------------------------
//! the compare function compares the degree of two nodes. This function is used
//! in the FELevelStructure::SortLevels function. Note the use of the static
//! this pointer to identify the level structure that called this function.

int FELevelStructure::compare(const void* e1, const void* e2)
{
	int n1 = *((int*) e1);
	int n2 = *((int*) e2);

	FENodeNodeList& L = *(m_pthis->m_pNL);

	return (L.Valence(n1) - L.Valence(n2));
}

//-----------------------------------------------------------------------------
//! This function creates a level structure rooted at node nroot of the graph L
//! The root is placed in level 0 and then all nodes adjacent to nodes that have
//! been placed in the level structure are placed in the next level.

void FELevelStructure::Create(FENodeNodeList& L, int nroot)
{
	int i, n, m, *pn;

	int N = L.Size();
	m_pNL = &L;

	// create the node array. This array will store
	// for each node to which level it belongs. We initialize
	// the array with -1 to indicate that no node has been
	// assigned to a level
	m_node.assign(N, -1);

	// place all nodes in a level. We use a queue for this
	// since we want nodes to be processed in a first come
	// first serve way. The root node is placed in level 0
	// and then all nodes adjacent to processed nodes are
	// placed in the next level.
	queue<int> NQ;
	int node = nroot;
	NQ.push(node);
	m_node[nroot] = 0;
	int level;
	while (NQ.size() > 0)
	{
		// get a node from the queue
		node = NQ.front(); NQ.pop();

		// get the next level index
		level = m_node[node] + 1;

		// now we have to place all neighbours on the queue
		n = L.Valence(node);
		pn = L.NodeList(node);
		for (i=0; i<n; ++i) 
		{
			// get an adjacent node
			m = pn[i];

			// make sure the node has not been assigned a level yet
			if (m_node[m] < 0)
			{
				// assign the node to a  level
				m_node[m] = level;

				// push the node on the queue
				NQ.push(m);
			}
		}
	}

	// determine the total nr of levels we generated
	int nlevels = 0;
	for (i=0; i<N; ++i)
		if (m_node[i] > nlevels) nlevels = m_node[i];
	++nlevels;

	// allocate levels data
	m_lval.assign(nlevels, 0);
	m_pl.resize(nlevels);
	m_nref.resize(N);

	// At this point it is important to realize that
	// there still may be nodes that are not assigned to 
	// a level. The reason is that these nodes are not connected
	// to the component that the root is part of. When looping
	// over all nodes (as below) it is important to explicitly
	// check whether the node belongs to the level structure or not.
	// Nodes that are not part of it have m_node[i] < 0

	// count the width of each level
	for (i=0; i<N; ++i)
	{
		if (m_node[i] >= 0) ++m_lval[ m_node[i] ];
	}

	// set nref pointers
	m_pl[0] = 0;
	m_nwidth = m_lval[0];
	for (i=1; i<nlevels; ++i)
	{
		m_pl[i] = m_pl[i-1] + m_lval[i-1];
		if (m_lval[i] > m_nwidth) m_nwidth = m_lval[i];
	}

	// reset the valence
	zero(m_lval);

	// fill the nref list and recount the width of each level
	// remember to see if the node is part of the level structure
	for (i=0; i<N; ++i)
	{
		n = m_node[i];
		if (n >= 0)
		{
			m_nref[ m_pl[n] + m_lval[n] ] = i;
			++m_lval[n];
		}
	}
}

//-----------------------------------------------------------------------------
//! This function sorts the nodes in each level from level l0 to l1 in order
//! of increasing degree. Note that the static this pointer gets set since the
//! compare function has to be static and therefore will not receive the this function
//! from the calling class member

void FELevelStructure::SortLevels(int l0, int l1)
{
	int n, *pn;
	m_pthis = this;
	for (int i=l0; i<=l1; ++i)
	{
		n = Valence(i);
		pn = NodeList(i);
		qsort(pn, n, sizeof(int), compare);
	}
}

//-----------------------------------------------------------------------------
//! This function merges two level structures into one. The merging algorithm below
//! usually returns a level structure that has a smaller width than either L1 or L2.
//! The swap variable will be set to true if second index of the node-pair array G 
//! was used of the first component.

void FELevelStructure::Merge(FELevelStructure& L1, FELevelStructure& L2, bool& bswap)
{
	int i, j, m, l1, l2;

	// get the node list that was used
	// to generate L1 and L2
	FENodeNodeList& NL = *L1.m_pNL;

	// store the node list for later use
	m_pNL = L1.m_pNL;

	// get the level depth of L1 and 
	// make sure that it is the same as L2's
	int k = L1.Depth();
	assert( k == L2.Depth()); 

	// get the level widths of L1 and L2
	int w1 = L1.Width();
	int w2 = L2.Width();

	// get the number of nodes
	int N = NL.Size();

	// create the valence array
	m_lval.assign(k, 0);

	// create the nodal array
	m_node.assign(N, -1);

	// In a moment nodes of the graph NL will be marked
	// as "removed". This will split the graph into 
	// several disconnected components. The comp array
	// will store for each node the component to which
	// it belongs.

	// create the component array
	// this array will store for each node
	// the index of the component
	// the nodes in the "0" component are the ones
	// that are marked as "removed" from the graph
	std::vector<int> comp; 
	comp.assign(N, -1);

	// fill the level pair array G will store for each
	// node the index of the level that it has in L1 and
	// the reversed index of level L2.
	// Note that l1 < 0 for nodes that are not part of L1
	// and similar for l2.
	vector<int>	G1(N), G2(N);
	for (i=0; i<N; ++i)
	{
		l1 = L1.m_node[i];
		l2 = L2.m_node[i];

		if ((l1>=0) && (l2>=0))
		{
			G1[i] = l1;
			G2[i] = k-1-l2;

			if (G1[i] == G2[i]) 
			{
				m_node[i] = l1;
				++m_lval[l1];
				comp[i] = 0;
			}
		}
		else
		{
			// node i does not belong to L1 or L2
			G1[i] = -1;
			G2[i] = -1;
		}
	}

	// Next, we find all connected components. A component
	// is created by finding an unprocessed node and attaching all
	// nodes that can be reached from this node. Remember that nodes
	// which are assigned to component 0 are considered as removed
	// from the graph
	int nc = 1;
	stack<int> NS;
	do
	{
		// find a node that has not been assigned to a component
		// (and that belongs to the graph (ie. node[i]>0))
		int node = -1;
		for (i=0; i<N; ++i)
		{
			if ((comp[i] == -1) && (L1.m_node[i] >= 0))
			{
				node = i;
				break;
			}
		}

		// if we didn't find one, we can stop
		if (node == -1) break;

		// assign the node and all adjacent nodes to a component
		NS.push(node);
		comp[node] = nc;
		while (!NS.empty())
		{
			node = NS.top(); NS.pop();

			int n = NL.Valence(node);
			int* pn = NL.NodeList(node);
			for (i=0; i<n; ++i)
				if (comp[pn[i]] == -1) 
				{
					comp[pn[i]] = nc;
					NS.push(pn[i]);
				}
		}

		// increase the component counter
		++nc;
	}
	while (true);

	// if we found more than 1 component than we'll have the process
	// each component seperatly. Note that we will only have one component
	// if all node-pairs of G are in the form (i,i). In that case all nodes
	// are placed in level 0, ie. are considered as "removed"
	if (nc>1)
	{
		// the vectors ni, hi and li will help us determine where to place
		// a node of a component
		std::vector<int> ni(k);
		std::vector<int> hi(k);
		std::vector<int> li(k);

		// Next, we have to create the components explicitly. The Comp arrays
		// will store for each component a list of nodes that belongs to this
		// component. We do this in two steps. First, we determine the sizes
		// of each component so that we can allocate the Comp array. Next, we
		// fill the Comp array by looping over all nodes and placing each node
		// in the correct component.

		// count the elements of the components
		// and determing the maximum component size
		std::vector<int> cc;
		cc.assign(nc, 0);
		int ncmax = 0;
		for (i=0; i<N; ++i) 
		{
			// make sure the node belongs to the graph L
			if (comp[i] >= 0)
			{
				++cc[ comp[i] ];
				if (cc[ comp[i] ] > ncmax) ncmax = cc[comp[i]];
			}
		}

		// create the component arrays
		std::vector< std::vector<int> > Comp(nc);
		for (i=0; i<nc; ++i) 
		{
			Comp[i].resize(cc[i]);
			cc[i] = 0;
		}

		// fill the components
		for (i=0; i<N; ++i)
		{
			m = comp[i];
			if (m>=0)
			{
				Comp[m][cc[m]] = i;
				++cc[m];
			}
		}

		// Next, we loop over all components and place each node of the 
		// component in a level. The level is determined by the node-pair
		// graph G. We can use either G[i][0] or G[i][1] for node i. Which
		// one gets picked determines on some criteria described below.
		// Note that we skip the "0" component since the nodes of this component
		// have already been assigned to a level.
		int ns = -1, nsm;
		for (i=1; i<nc; ++i)
		{
			// we need to loop over all components in order of increasing size.
			// In stead of actually sorting the components we just find the next
			// smallest component by looping again over all components. At the bottom
			// we set cc[i] of the smallest component to -1 so that it can't get
			// picked anymore.
			nsm = ncmax;
			int ncomp = -1;
			for (j=1; j<nc; ++j)
			{
				if ((cc[j] >= ns) && (cc[j] <= nsm))
				{
					ncomp = j;
					nsm = cc[j];
				}
			}
			ns = nsm;
			assert(ncomp > 0);

			// calculate the vectors ni, hi, li
			// ni[i] = number of nodes in level i
			// li[i] = ni + nodes that would be put in level i based on G[node][0]
			// hi[i] = ni + nodes that would be put in level i based on G[node][1]
			for (j=0; j<k; ++j)
			{
				li[j] = hi[j] = ni[j] = m_lval[j];
			}

			std::vector<int>& pn = Comp[ncomp];
			for (j=0; j<cc[ncomp]; ++j)
			{
				m = pn[j];
				++hi[ G1[m] ];
				++li[ G2[m] ];
			}

			// next, we determint h0 and l0, where
			// h0 = max{hi[i] : hi[i] - ni[i] > 0} and similar for l0
			int h0 = -1, l0 = -1;
			for (j=0; j<k; ++j)
			{
				if ((hi[j] > h0) && (hi[j] > ni[j])) h0 = hi[j];
				if ((li[j] > l0) && (li[j] > ni[j])) l0 = li[j];
			}
			assert(h0 >= 0);
			assert(l0 >= 0);

			// Next we can assign the nodes of the component to a level.
			// the criteria are:
			// if (h0 < l0) place the node in level G[i][0]
			// if (l0 < h0) place the node in level G[i][1]
			// if (l0 == h0) use the level of the rooted level structure 
			//   with smallest width
			if ((h0 < l0) || ((h0 == l0) && (w1 <= w2)))
			{
				std::vector<int>& pn = Comp[ncomp];
				for (j=0; j<cc[ncomp]; ++j)
				{
					m = pn[j];
					m_node[m] = G1[m];
					++m_lval[ G1[m] ];
				}
				if (i==1) bswap = false;
			}
			else if ((l0 < h0) || ((h0 == l0) && (w1 > w2)))
			{
				std::vector<int>& pn = Comp[ncomp];
				for (j=0; j<cc[ncomp]; ++j)
				{
					m = pn[j];
					m_node[m] = G2[m];
					++m_lval[ G2[m] ];
				}
				if (i==1) bswap = true;
			}

			// set cc[ncomp] to an impossible value so that it won't get picked any more
			// in the search for the next smallest component
			cc[ncomp] = -1;
		}
	}

	// All nodes are now placed in a level and we now the size of each level
	// We can now finish allocating and initializing the rest of the level structure data

	// allocate levels data
	m_pl.resize(k);
	m_nref.resize(N);

	// set nref pointers
	m_pl[0] = 0;
	m_nwidth = m_lval[0];
	for (i=1; i<k; ++i)
	{
		m_pl[i] = m_pl[i-1] + m_lval[i-1];
		if (m_lval[i] > m_nwidth) m_nwidth = m_lval[i];
	}

	// reset the valence
	zero(m_lval);

	// fill the nref list
	// Remember to check to see if the node belongs to the graph
	for (i=0; i<N; ++i)
	{
		m = m_node[i];
		if (m>=0)
		{
			m_nref[ m_pl[m] + m_lval[m] ] = i;
			++m_lval[m];
		}
	}
}
