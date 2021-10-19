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
#include "FEElemElemList.h"
#include "FENodeElemList.h"
#include "FESolidDomain.h"
#include "FESurface.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
FEElemElemList::FEElemElemList(void)
{
}

//-----------------------------------------------------------------------------
FEElemElemList::~FEElemElemList(void)
{
}

//-----------------------------------------------------------------------------
void FEElemElemList::Init()
{
	FEMesh& m = *m_pmesh;

	// allocate storage
	int NE = m.Elements();
	m_ref.resize(NE);

	// count nr of neighbors
	int NN = 0, n = 0, nf;
	m_ref[0] = 0;
	for (int i=0; i<m.Domains(); ++i)
	{
		FEDomain& dom = m.Domain(i);
		for (int j=0; j<dom.Elements(); ++j, ++n)
		{
			FEElement& el = dom.ElementRef(j);
			nf = el.Faces();
			if (n != 0) m_ref[n] = m_ref[n-1] + nf;
			NN += nf;
		}
	}

	m_pel.resize(NN);
	m_peli.resize(NN);

	// TODO: do this for shells as well (if we have to)
}

//-----------------------------------------------------------------------------
bool FEElemElemList::Create(FEMesh* pmesh)
{
	// store a pointer to the mesh
	m_pmesh = pmesh;
	FEMesh& m = *m_pmesh;

	// initialize data structures
	Init();

	// create the node element list
	FENodeElemList NEL;
	NEL.Create(m);

	// loop over all solid elements first
	int en0[FEElement::MAX_NODES], en1[FEElement::MAX_NODES], n0, n1, M = 0;
	int nf0, nf1;
	for (int nd=0; nd<m.Domains(); ++nd)
	{
		FEDomain& dom = m.Domain(nd);
		for (int i=0; i<dom.Elements(); ++i)
		{
			FEElement& el = dom.ElementRef(i);
			
			// get the number of neighbors
			nf0 = el.Faces();

			// loop over all neighbors
			for (int j=0; j<nf0; ++j, ++M)
			{
				// get the face nodes
				n0 = el.GetFace(j, en0);

				// find the neighbor element
				m_pel[M] = 0;
				m_peli[M] = -1;

				// loop over all possible candidates
				int nval = NEL.Valence(en0[0]);
				FEElement** pne = NEL.ElementList(en0[0]);
				int* pnei = NEL.ElementIndexList(en0[0]);
				for (int k=0; k<nval; ++k)
				{
					// make sure we don't compare the current element
					if (pne[k] != &el)
					{
						// get the number of faces
						nf1 = pne[k]->Faces();

						// see if any of these faces match en0
						for (int l=0; l<nf1; ++l)
						{
							n1 = pne[k]->GetFace(l, en1);

							// make sure the faces have the same nr of nodes
							if (n1 == n0)
							{
								// check triangles
								if ((n0 == 3) || (n0 == 6) || (n0 ==7))
								{
									if (((en0[0] == en1[0]) || (en0[0] == en1[1]) || (en0[0] == en1[2])) &&
										((en0[1] == en1[0]) || (en0[1] == en1[1]) || (en0[1] == en1[2])) &&
										((en0[2] == en1[0]) || (en0[2] == en1[1]) || (en0[2] == en1[2])))
									{
										// found it!
										m_pel[M] = pne[k];
										m_peli[M] = pnei[k];
										break;
									}
								}
								// check quads
								else if ((n0 == 4) || (n0 == 8) || (n0 == 9))
								{
									if (((en0[0] == en1[0]) || (en0[0] == en1[1]) || (en0[0] == en1[2]) || (en0[0] == en1[3])) &&
										((en0[1] == en1[0]) || (en0[1] == en1[1]) || (en0[1] == en1[2]) || (en0[1] == en1[3])) &&
										((en0[2] == en1[0]) || (en0[2] == en1[1]) || (en0[2] == en1[2]) || (en0[2] == en1[3])) &&
										((en0[3] == en1[0]) || (en0[3] == en1[1]) || (en0[3] == en1[2]) || (en0[3] == en1[3])))
									{
										// found it!
										m_pel[M] = pne[k];
										m_peli[M] = pnei[k];
										break;
									}
								}
							}

							if (m_pel[M] != 0) break;
						}
					}
				}
			}
		}
	}

	// TODO: do the same for shells

	return true;
}

//-----------------------------------------------------------------------------
//! Find the element neighbors for a surface. In this case, the elements are
//! surface elements (i.e. FESurfaceElement).
bool FEElemElemList::Create(const FESurface* psurf)
{
	// allocate storage
	int NE = psurf->Elements();
	m_ref.resize(NE);

	// count nr of neighbors
	int NN = 0;
	m_ref[0] = 0;
	for (int j=0; j<NE; ++j)
	{
		const FESurfaceElement& el = psurf->Element(j);

		int nf = el.facet_edges();

		if (j != NE-1) m_ref[j+1] = m_ref[j] + nf;
		NN += nf;
	}

	m_pel.resize(NN);

	// create the node element list
	FENodeElemList NEL;
	NEL.Create(*psurf);

	// loop over all facets
	int en0[3], en1[3], M = 0;
	int nf0, nf1;
	for (int i=0; i<NE; ++i)
	{
		const FESurfaceElement& el = psurf->Element(i);
			
		// get the number of neighbors
		nf0 = el.facet_edges();

		// loop over all neighbors
		for (int j=0; j<nf0; ++j, ++M)
		{
			// get the edge nodes
			el.facet_edge(j, en0);

			// find the neighbor element
			m_pel[M] = 0;

			// loop over all possible candidates
			int nval = NEL.Valence(en0[0]);
			FEElement** pne = NEL.ElementList(en0[0]);
			for (int k=0; k<nval; ++k)
			{
				// make sure we don't compare the current element
				if (pne[k] != &el)
				{
					// get the number of edges
					FESurfaceElement& me = dynamic_cast<FESurfaceElement&>(*pne[k]);
					nf1 = me.facet_edges();

					// see if any of these edges match en0
					for (int l=0; l<nf1; ++l)
					{
						me.facet_edge(l, en1);

						if (((en0[0] == en1[0]) || (en0[0] == en1[1])) &&
							((en0[1] == en1[0]) || (en0[1] == en1[1])))
							{
								// found it!
								m_pel[M] = pne[k];
								break;
							}
					}

					if (m_pel[M] != 0) break;
				}
			}
		}
	}

	return true;
}
