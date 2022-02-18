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
#include "FESPRProjection.h"
#include "FESolidDomain.h"
#include "FEMesh.h"
using namespace std;

//-------------------------------------------------------------------------------------------------
FESPRProjection::FESPRProjection()
{
	m_p = -1;
}

//-------------------------------------------------------------------------------------------------
void FESPRProjection::SetInterpolationOrder(int p)
{
	m_p = p;
}

//-------------------------------------------------------------------------------------------------
//! Projects the integration point data, stored in d, onto the nodes of the domain.
//! The result is stored in o.
void FESPRProjection::Project(FESolidDomain& dom, const vector< vector<double> >& d, vector<double>& o)
{
	// get the mesh
	FEMesh& mesh = *dom.GetMesh();
	int NN = dom.Nodes();

	// allocate output array
	o.assign(NN, 0.0);

	// check element type
	int NDOF = -1;	// number of degrees of freedom of polynomial
	int NCN  = -1;	// number of corner nodes
	int nshape = dom.GetElementShape();
	switch (nshape)
	{
	case ET_TET4   :
	case ET_TET5   : { NDOF =  4; NCN = 4; } break;
	case ET_TET10  : { NDOF =  4; NCN = 4; } break;
	case ET_TET15  : 
		{
			NDOF = (m_p == 1 ? 4 : 10); 
			NCN = 4; 
		}
		break;
	case ET_TET20 : { NDOF = 10; NCN = 4; } break;
	case ET_HEX8  : { NDOF =  7; NCN = 8; } break;
	case ET_HEX20 : { NDOF = (m_p == 1 ? 7 : 10); NCN = 8; } break;
	case ET_HEX27 : { NDOF = (m_p == 1 ? 7 : 10); NCN = 8; } break;
	default:
		return;
	}

	// we keep a tag array to keep track of which nodes we processed
	int NM = mesh.Nodes();
	vector<int> tag; tag.assign(NM, 0);

	// for higher order elements
	// we need to make sure that we don't process the edge nodes
	// we assume here that the first NCN nodes of the element
	// are the corner nodes and that all other nodes are edge or interior nodes
	int NE = dom.Elements();
	for (int i=0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int ne = el.Nodes();
		for (int j=NCN; j<ne; ++j) tag[el.m_node[j]] = 2;
	}

	// this array will store the results
	vector<double> val;
	val.assign(NM, 0.0);

	// build the node-element-list. This will define our patches
	FENodeElemList NEL;
	NEL.Create(dom);

	// loop over all nodes
	for (int i=0; i<NN; ++i)
	{
		// get the node
		FENode& node = dom.Node(i);
		int in = dom.NodeIndex(i);

		// don't loop over edge nodes (edge or interior nodes have a tag > 1)
		if (tag[in] <= 1)
		{
			// get the nodal position
			vec3d rc = node.m_rt;

			// get the element patch
			int ne = NEL.Valence(in);
			FEElement** ppe = NEL.ElementList(in);
			int* pei = NEL.ElementIndexList(in);

			// setup the A-matrix
			vector<double> pk(NDOF);
			matrix A(NDOF,NDOF); A.zero();
			int m = 0;
			for (int j=0; j<ne; ++j)
			{
				FEElement& el = *(ppe[j]);

				int nint = el.GaussPoints();
				for (int n=0; n<nint; ++n, ++m)
				{
					FEMaterialPoint& mp = *el.GetMaterialPoint(n);
					vec3d r = mp.m_rt - rc;
					pk[0] = 1.0; pk[1] = r.x; pk[2] = r.y; pk[3] = r.z;
					if (NDOF >=  7) { pk[4] = r.x*r.y; pk[5] = r.y*r.z; pk[6] = r.x*r.z; }
					if (NDOF >= 10) { pk[7] = r.x*r.x; pk[8] = r.y*r.y; pk[9] = r.z*r.z; }
					A += outer_product(pk);
				}
			}

			// invert matrix and make sure condition number is good enough
			matrix Ai = A.inverse();

			// make sure we have enough sampling points
			if (m > NDOF + 1)
			{
				vector<double> b; b.assign(NDOF,0.0);
				for (int j=0; j<ne; ++j)
				{
					FEElement& el = *(ppe[j]);
					const vector<double>& ed = d[pei[j]];

					assert(ppe[j] == &dom.Element(pei[j]));

					int nint = el.GaussPoints();
					for (int n=0; n<nint; ++n)
					{
						FEMaterialPoint& mp = *el.GetMaterialPoint(n);
						vec3d r = mp.m_rt - rc;
						pk[0] = 1.0; pk[1] = r.x; pk[2] = r.y; pk[3] = r.z;
						if (NDOF >=  7) { pk[4] = r.x*r.y; pk[5] = r.y*r.z; pk[6] = r.x*r.z; }
						if (NDOF >= 10) { pk[7] = r.x*r.x; pk[8] = r.y*r.y; pk[9] = r.z*r.z; }

						double s = ed[n];
						for (int k=0; k<NDOF; k++) b[k] += s*pk[k];
					}
				}

				// solve the linear system
				vector<double> c = Ai*b;

				// tag this node as processed
				tag[in] = 1;

				// store result
				val[in] = c[0];

				// loop over all unprocessed nodes of this patch
				for (int j=0; j<ne; ++j)
				{
					FEElement& el = *(ppe[j]);
					int en = el.Nodes();
					for (int k=0; k<en; ++k)
					{
						int em = el.m_node[k];
						if (tag[em] != 1)
						{
							vec3d r = mesh.Node(em).m_rt - rc;
							pk[0] = 1.0; pk[1] = r.x; pk[2] = r.y; pk[3] = r.z;
							if (NDOF >=  7) { pk[4] = r.x*r.y; pk[5] = r.y*r.z; pk[6] = r.x*r.z; }
							if (NDOF >= 10) { pk[7] = r.x*r.x; pk[8] = r.y*r.y; pk[9] = r.z*r.z; }

							// calculate the value for this node
							double v = 0;
							for (int l=0; l<NDOF; ++l) v += pk[l]*c[l];

							// for edge nodes, we need to keep track of how often we visit this node
							// Therefore we increment the tag.
							// (remember that the tag started at 2 for edge/interior nodes)
							if (tag[em] >= 2)
							{
								tag[em]++;
								val[em] += v;
							}
							else val[em] = v;
						}
					}
				}
			}
		}
	}

	// copy results to archive
	for (int i=0; i<NN; ++i)
	{
		int in = dom.NodeIndex(i);
		double s = val[in];

		// for edge nodes we need to average
		// (remember that the tag started at 2 for edge/interior nodes)
		if (tag[in] >= 2)
		{
//			assert(tag[in] > 2);	// all edges nodes must be visited at least once!
			int l = tag[in]-2;
			if (l > 0) s /= (double) (l);
		}
		
		o[i] = s;
	}
}
