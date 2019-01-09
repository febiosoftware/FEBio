#include "stdafx.h"
#include "FEFaceList.h"
#include "FEMesh.h"
#include "FEDomain.h"
#include "FEElemElemList.h"
#include "FEElementList.h"

bool FEFaceList::FACE::IsEqual(int* n) const
{
	if (ntype == 3)
	{
		if ((node[0] != n[0]) && (node[0] != n[1]) && (node[0] != n[2])) return false;
		if ((node[1] != n[0]) && (node[1] != n[1]) && (node[1] != n[2])) return false;
		if ((node[2] != n[0]) && (node[2] != n[1]) && (node[2] != n[2])) return false;
	}
	else if (ntype == 4)
	{
		if ((node[0] != n[0]) && (node[0] != n[1]) && (node[0] != n[2]) && (node[0] != n[3])) return false;
		if ((node[1] != n[0]) && (node[1] != n[1]) && (node[1] != n[2]) && (node[1] != n[3])) return false;
		if ((node[2] != n[0]) && (node[2] != n[1]) && (node[2] != n[2]) && (node[2] != n[3])) return false;
		if ((node[3] != n[0]) && (node[3] != n[1]) && (node[3] != n[2]) && (node[3] != n[3])) return false;
	}
	else 
	{
		assert(false);
		return false;
	}

	return true;
}


FEFaceList::FEFaceList() : m_mesh(nullptr)
{

}

int FEFaceList::Faces() const
{
	return (int) m_faceList.size();
}

const FEFaceList::FACE& FEFaceList::Face(int i) const
{
	return m_faceList[i];
}

FEMesh* FEFaceList::GetMesh()
{
	return m_mesh;
}

bool FEFaceList::Create(FEMesh* pmesh)
{
	if (pmesh == nullptr) return false;
	m_mesh = pmesh;
	FEMesh& mesh = *pmesh;

	// create the element neighbor list
	FEElemElemList EEL;
	EEL.Create(pmesh);

	// get the number of elements in this mesh
	int NE = mesh.Elements();

	// count the number of facets we have to create
	int NF = 0;
	FEElementList EL(mesh);
	FEElementList::iterator it = EL.begin();
	for (int i = 0; i<NE; ++i, ++it)
	{
		FEElement& el = *it;
		int nf = mesh.Faces(el);
		for (int j = 0; j<nf; ++j)
		{
			FEElement* pen = EEL.Neighbor(i, j);
			if (pen == 0) ++NF;
			if ((pen != 0) && (el.GetID() < pen->GetID())) ++NF;
		}
	}

	// create the facet list
	m_faceList.resize(NF);

	// build the facets
	int face[FEElement::MAX_NODES];
	NF = 0;
	it = EL.begin();
	for (int i = 0; i<NE; ++i, ++it)
	{
		FEElement& el = *it;
		int nf = mesh.Faces(el);
		for (int j = 0; j<nf; ++j)
		{
			FEElement* pen = EEL.Neighbor(i, j);
			if ((pen == 0) || ((pen != 0) && (el.GetID() < pen->GetID())))
			{
				FACE& se = m_faceList[NF++];
				mesh.GetFace(el, j, face);

				switch (el.Shape())
				{
				case ET_HEX8:
					se.ntype = 4;
					break;
				case ET_TET4:
					se.ntype = 3;
					break;
				default:
					assert(false);
				}

				int nn = se.ntype;
				for (int k = 0; k<nn; ++k)
				{
					se.node[k] = face[k];
				}
			}
		}
	}

	return true;
}

//=============================================================================

FEElementFaceList::FEElementFaceList()
{

}

//NOTE: only works for hex elements
bool FEElementFaceList::Create(FEElementList& elemList, FEFaceList& faceList)
{
	FEMesh& mesh = *faceList.GetMesh();

	const int FTET[4][3] = {
		{ 0, 1, 3},
		{ 1, 2, 3},
		{ 2, 0, 3},
		{ 2, 1, 0}
	};

	const int FHEX[6][4] = { 
		{ 0, 1, 5, 4 },
		{ 1, 2, 6, 5 },
		{ 2, 3, 7, 6 },
		{ 3, 0, 4, 7 },
		{ 3, 2, 1, 0 },
		{ 4, 5, 6, 7 }};

	// build a node face table for FT to facilitate searching
	int NN = mesh.Nodes();
	vector<vector<int> > NFT; NFT.resize(NN);
	for (int i = 0; i<faceList.Faces(); ++i)
	{
		const FEFaceList::FACE& f = faceList.Face(i);
		NFT[f.node[0]].push_back(i);
		NFT[f.node[1]].push_back(i);
		NFT[f.node[2]].push_back(i);
		if (f.ntype == 4) NFT[f.node[3]].push_back(i);
	}

	int NE = mesh.Elements();
	m_EFL.resize(NE);
	int i = 0;
	for (FEElementList::iterator it = elemList.begin(); it != elemList.end(); ++it, ++i)
	{
		const FEElement& el = *it;
		vector<int>& EFLi = m_EFL[i];
		if (el.Shape() == FE_Element_Shape::ET_TET4)
		{
			EFLi.resize(4);
			for (int j = 0; j < 4; ++j)
			{
				int fj[3] = { el.m_node[FTET[j][0]], el.m_node[FTET[j][1]], el.m_node[FTET[j][2]] };
				EFLi[j] = -1;
				vector<int>& nfi = NFT[fj[0]];
				for (int k = 0; k<(int)nfi.size(); ++k)
				{
					const FEFaceList::FACE& fk = faceList.Face(nfi[k]);
					if (fk.IsEqual(fj))
					{
						EFLi[j] = nfi[k];
						break;
					}
				}
				assert(EFLi[j] != -1);
			}
		}
		else if (el.Shape() == FE_Element_Shape::ET_HEX8)
		{
			EFLi.resize(6);
			for (int j = 0; j < 6; ++j)
			{
				int fj[4] = { el.m_node[FHEX[j][0]], el.m_node[FHEX[j][1]], el.m_node[FHEX[j][2]], el.m_node[FHEX[j][3]] };
				EFLi[j] = -1;
				vector<int>& nfi = NFT[fj[0]];
				for (int k = 0; k<(int)nfi.size(); ++k)
				{
					const FEFaceList::FACE& fk = faceList.Face(nfi[k]);
					if (fk.IsEqual(fj))
					{
						EFLi[j] = nfi[k];
						break;
					}
				}
				assert(EFLi[j] != -1);
			}
		}
		else return false;
	}
	return true;
}

int FEElementFaceList::Faces(int elem) const
{
	return (int)m_EFL[elem].size();
}

std::vector<int> FEElementFaceList::FaceList(int elem) const
{
	return m_EFL[elem];
}
