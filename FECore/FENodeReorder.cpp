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
#include "FENodeReorder.h"
#include "FEMesh.h"
#include <stack>
using namespace std;

//-----------------------------------------------------------------------------

FENodeReorder::FENodeReorder()
{

}

FENodeReorder::~FENodeReorder()
{

}

//-----------------------------------------------------------------------------
//! This function applies the node-reordering algorithm to a particular mesh.
//! The new node number is stored in P. To be precise, P stores for each new
//! node the old node that corresponds to this node.

void FENodeReorder::Apply(FEMesh& mesh, vector<int>& P)
{
	int i, j, n, l, m;
	int* pn;

	// get the nr of nodes
	int N = mesh.Nodes();

	// initialize the permutation vector
	// Set all entries to -1 to indicate that
	// no node has been given a new number yet.
	vector<int> Q; Q.assign(N, -1);

	// create the node-node list
	// this list stores for each node of the mesh
	// a list of nodes that are adjacent to it.
	FENodeNodeList NL;
	NL.Create(mesh);

	// sort the nodelist in order of increasing degree
	NL.Sort();

	// Levelstructures are used to group all nodes
	// in levels
	FELevelStructure Lv, Lu, Ls, L;

	// this swap variable will tell you if 
	// the node numbering needs to be reversed at the end
	bool bswap;

	// The mesh may comprise of several disconnected
	// components. We therefor have to apply the algorithm
	// several times until all components are processed.
	int neq = 0, neq0;
	while (true)
	{
		// To identify an unprocessed component we try to
		// find an unprocessed node (nv). An unprocessed
		// node will have Q[i]<0. We also require this node 
		// to be of minimal degree.
		int valmin = 0x7fffffff;	// initialize the min valence to a ridiculously large value
		int nv = -1, nu;
		int nval;
		for (i=0; i<N; ++i)
		{
			nval = NL.Valence(i);
			if ((Q[i] < 0) && (nval < valmin)) 
			{ 
				valmin = nval;
				nv = i; 
			}
		}

		// if we havn't found an unprocessed node, we can stop
		if (nv == -1) break;

		// store the starting equation number
		// we need this in case we need to reverse the 
		// node ordering
		neq0 = neq;

		// we now will find the pseudo-diameter of the node graph (i.e. NodeNodeList),
		// which identifies to nodes at the end of the diameter,
		// namely nv, and nu. We will require the level structure for
		// nu to be of minimal width, which is stored in wmin
		int wmin = 0x7fffffff;

		// Let's find the pseudo-diameter
		// if we found one bok will be true
		bool bok;
		do
		{
			// let's be optimistic
			bok = true;

			// create a level structure rooted at nv
			Lv.Create(NL, nv);

			// sort the last level in order of increasing degree
			l = Lv.Depth() - 1;
			Lv.SortLevels(l, l);

			// get the list of nodes at the last level of Lv
			n = Lv.Valence(l);
			pn = Lv.NodeList(l);
			// loop over all the nodes
			for (i=0; i<n; ++i)
			{
				// get a possible candidate for nu
				nu = pn[i];

				// create a level structure rooted at u
				Ls.Create(NL, nu);

				// make sure that the depth of Ls is not
				// greater than that of Lv
				if (Ls.Depth() > Lv.Depth())
				{
					// Oh, oh, the depth of Ls is greater than
					// that of nv, so replace nv with nu and try
					// again.
					nv = nu;
					wmin = 0x7fffffff;	// reset the min width
					bok = false;
					break;
				}
				else
				{
					// store the level stucture with minimal width
					if (Ls.Width() < wmin)
					{
						Lu = Ls;	// store the level structure
						wmin = Ls.Width();	// store the minimum width
					}
				}
			}
		}
		while (!bok);

		// we have found the pseudo-diamter u,v
		// Next, merge the level structures into one
		// This new level structure will usually have a width
		// that is smaller than either Lv or Lu.
		// make sure we start with the node of lowest degree
		if (NL.Valence(nu) < NL.Valence(nv))
		{
			// swap u and v
			nu ^= nv;
			nv ^= nu;
			nu ^= nv;
			L.Merge(Lu, Lv, bswap);
		}
		else
		{
			L.Merge(Lv, Lu, bswap);
			bswap = !bswap;
		}

		// sort all levels of L in order of increasing degree
		L.SortLevels(0, L.Depth()-1);

		// get the width of L. This width is the largest width
		// of any level of L
		int lmax = L.Width();

		// the following stack and vectors will assist us in
		// renumbering the nodes.
		stack<int> NQ;
		vector<int> Vi(0, lmax);
		vector<int> Vip1(0, lmax);

		// Using the level structure L we can now go ahead and define
		// the new node numbering. The basic idea is to loop over all
		// levels and assign the node numbers level by level. In each
		// level the numbers adjacent to the lowest number node in the same
		// level gets numbered first. Then the nodes adjacent to the lowest
		// numbered nodes in the previous level gets numbered next.
		int nw;
		Q[nv] = neq++;
		Vi.push_back(nv);
		for (l=0; l<L.Depth(); ++l)
		{
			// assign numbers to level l
			// The Vi array stores the nodes of level l that have been
			// numbered in order of increasing node number.
			// We loop over all these nodes and assign a 
			// node number to unnumbered adjacent nodes until
			// all nodes of level l have been numbered 
			for(i=0; i<(int) Vi.size(); ++i)
			{
				nw = Vi[i];

				// loop over all unassigned nodes adjacent to nw
				n = NL.Valence(nw);
				pn = NL.NodeList(nw);
				for (j=0; j<n; ++j)
				{
					m = pn[j];
					if ((L.NodeLevel(m) == l) && (Q[m] < 0))
					{
						Q[m] = neq++;
						Vi.push_back(m);
					}
				}

				// If we are about to terminate the loop make sure
				// that all nodes of level l have been assigned 
				// a new node number
				if (i==(Vi.size() - 1))
				{
					// check to see if all nodes in level l are numbered
					n = L.Valence(l);
					pn = L.NodeList(l);
					for (j=0; j<n; ++j)
					{
						m = pn[j];
						if (Q[m] < 0)
						{
							// Aha, we found an unnumbered node in this level
							// so place in the array an continue the loop over i
							Vi.push_back(m);
							Q[m] = neq++;
							break;
						}
					}
				}
			}

			// If this is not the last level, we go ahead and number
			// all nodes in level l+1 that are adjacent to numbered nodes
			// of level l
			if (l<L.Depth()-1)
			{
				Vip1.clear();

				// assign numbers to the next level
				for (i=0; i<(int) Vi.size(); ++i)
				{
					nw = Vi[i];
					n = NL.Valence(nw);
					pn = NL.NodeList(nw);
					for (j=0; j<n; ++j)
					{
						m = pn[j];
						if ((L.NodeLevel(m) == l+1) && (Q[m] < 0))
						{
							Q[m] = neq++;
							Vip1.push_back(m);
						}
					}
				}

				// copy the numbered nodes of level l+1 into Vi
				Vi = Vip1;
			}
		}

		// we have assigned new node numbers 
		// see if we need to invert the numbering
		if (bswap)
		{
			for (i=0; i<N; ++i) 
			{
				if (Q[i] >= neq0)
					Q[i] = neq - 1 - Q[i] + neq0;
			}
		}
	}

	// The Q array stores for each old node its new node number
	// but we actually need the inverse permutation. That is, for
	// each new node, the old node number that is associated with it.
	// This array is stored in P.
	P.resize(N);
	for (i=0; i<N; ++i)
	{
		P[Q[i]] = i;
	}
}
