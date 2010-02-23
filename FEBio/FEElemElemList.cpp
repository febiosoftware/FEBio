#include "stdafx.h"
#include "FEElemElemList.h"
#include "FENodeElemList.h"

FEElemElemList::FEElemElemList(void)
{
}

FEElemElemList::~FEElemElemList(void)
{
}

void FEElemElemList::Init()
{
	int i;
	FEMesh& m = *m_pmesh;

	FESolidDomain& bd = m.SolidDomain();
	FEShellDomain& sd = m.ShellDomain();

	// get total nr of elements
	int NE = bd.Elements() + sd.Elements();

	// allocate storage
	m_ref.resize(NE);

	// count nr of neighbors
	int NN = 0, n = 0, nf;
	m_ref[0] = 0;
	for (i=0; i<bd.Elements(); ++i, ++n)
	{
		FESolidElement& el = bd.Element(i);
		nf = m.Faces(el);
		if (n != 0) m_ref[n] = m_ref[n-1] + nf;
		NN += nf;
	}

	for (i=0; i<sd.Elements(); ++i, ++n)
	{
		FEShellElement& el = sd.Element(i);
		nf = m.Faces(el);
		if (n != 0) m_ref[n] = m_ref[n-1] + nf;
		NN += nf;
	}

	m_pel.resize(NN);
}

void FEElemElemList::Create(FEMesh* pmesh)
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
	int i, j, k, l;
	int n = 0, en0[4], en1[4], n0, n1, M = 0;
	int nf0, nf1;
	FESolidDomain& bd = m.SolidDomain();
	for (i=0; i<bd.Elements(); ++i, ++n)
	{
		FEElement& el = bd.Element(i);

		// get the number of neighbors
		nf0 = m.Faces(el);

		// loop over all neighbors
		for (j=0; j<nf0; ++j, ++M)
		{
			// get the face nodes
			n0 = m.GetFace(el, j, en0);

			// find the neighbor element
			m_pel[M] = 0;

			// loop over all possible candidates
			int nval = NEL.Valence(en0[0]);
			FEElement** pne = NEL.ElementList(en0[0]);
			for (k=0; k<nval; ++k)
			{
				// make sure we don't compare the current element
				if (pne[k] != &el)
				{
					// get the number of faces
					nf1 = m.Faces(*pne[k]);

					// see if any of these faces match en0
					for (l=0; l<nf1; ++l)
					{
						n1 = m.GetFace(*pne[k], l, en1);

						// make sure the faces have the same nr of nodes
						if (n1 == n0)
						{
							// check triangles
							if (n0 == 3)
							{
								if (((en0[0] == en1[0]) || (en0[0] == en1[1]) || (en0[0] == en1[2])) &&
									((en0[1] == en1[0]) || (en0[1] == en1[1]) || (en0[1] == en1[2])) &&
									((en0[2] == en1[0]) || (en0[2] == en1[1]) || (en0[2] == en1[2])))
								{
									// found it!
									m_pel[M] = pne[k];
									break;
								}
							}
							// check quads
							else if (n0 == 4)
							{
								if (((en0[0] == en1[0]) || (en0[0] == en1[1]) || (en0[0] == en1[2]) || (en0[0] == en1[3])) &&
									((en0[1] == en1[0]) || (en0[1] == en1[1]) || (en0[1] == en1[2]) || (en0[1] == en1[3])) &&
									((en0[2] == en1[0]) || (en0[2] == en1[1]) || (en0[2] == en1[2]) || (en0[2] == en1[3])) &&
									((en0[3] == en1[0]) || (en0[3] == en1[1]) || (en0[3] == en1[2]) || (en0[3] == en1[3])))
								{
									// found it!
									m_pel[M] = pne[k];
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

	// TODO: do the same for shells
}
