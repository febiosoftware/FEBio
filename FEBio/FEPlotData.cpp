#include "stdafx.h"
#include "FEPlotData.h"
#include "fem.h"

//-----------------------------------------------------------------------------
int FEPlotData::VarSize(Var_Type t)
{
	int ndata = 0;
	switch (DataType())
	{
	case FLOAT: ndata = 1; break;
	case VEC3F: ndata = 3; break;
	case MAT3FS: ndata = 6; break;
	}
	assert(ndata);
	return ndata;
}

//-----------------------------------------------------------------------------
void FENodeData::Save(FEM &fem, Archive& ar)
{
	// loop over all node sets
	// write now there is only one, namely the master node set
	// so we just pass the mesh
	int ndata = VarSize(DataType());

	int N = fem.m_mesh.Nodes();
	vector<float> a; a.reserve(ndata*N);
	if (Save(fem.m_mesh, a))
	{
		assert(a.size() == N*ndata);
		ar.WriteChunk(0, a);
	}
}

//-----------------------------------------------------------------------------
void FEDomainData::Save(FEM &fem, Archive& ar)
{
	// loop over all domains
	FEMesh& m = fem.m_mesh;
	int ND = m.Domains();
	for (int i=0; i<ND; ++i)
	{
		// get the domain
		FEDomain& D = m.Domain(i);

		// calculate the size of the data vector
		int nsize = VarSize(DataType());
		switch (m_sfmt)
		{
		case FMT_NODE: nsize *= D.Nodes(); break;
		case FMT_ITEM: nsize *= D.Elements(); break;
		case FMT_MULT:
			{
				// since all elements have the same type within a domain
				// we just grab the number of nodes of the first element 
				// to figure out how much storage we need
				FEElement& e = D.ElementRef(0);
				int n = e.Nodes();
				nsize *= n*D.Elements();
			}
			break;
		default:
			assert(false);
		}

		// fill data vector and save
		vector<float> a; 
		a.reserve(nsize);
		if (Save(D, a))
		{
			assert(a.size() == nsize);
			ar.WriteChunk(i+1, a);
		}
	}
}

//-----------------------------------------------------------------------------
void FESurfaceData::Save(FEM &fem, Archive& ar)
{
	// loop over all surfaces
	FEMesh& m = fem.m_mesh;
	int NS = m.Surfaces();
	for (int i=0; i<NS; ++i)
	{
		FESurface& S = m.Surface(i);

		// Determine data size.
		// Note that for the FMT_MULT case we are 
		// assuming four data entries per facet
		// regardless of the nr of nodes a facet really has
		// this is because for surfaces, all elements are not
		// necessarily of the same type
		int nsize = VarSize(DataType());
		switch (m_sfmt)
		{
		case FMT_NODE: nsize *= S.Nodes(); break;
		case FMT_ITEM: nsize *= S.Elements(); break;
		case FMT_MULT: nsize *= 4*S.Elements(); break;
		default:
			assert(false);
		}

		// save data
		vector<float> a; a.reserve(nsize);
		if (Save(S, a))
		{
			assert(a.size() == nsize);
			ar.WriteChunk(i+1, a);
		}
	}
}
