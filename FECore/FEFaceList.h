#pragma once
#include <vector>
#include "fecore_api.h"

class FEMesh;
class FEElementList;
class FEElemElemList;

class FECORE_API FEFaceList
{
public:
	struct FACE
	{
		int		ntype;		// 3 = triangle, 4 = quad
		int		node[4];
		int		nsurf;		// 1 if facet is on surface 
		int		nbr[4];		// neighbor list

		bool IsEqual(int* n) const;

		bool HasEdge(int a, int b) const;
	};

public:
	FEFaceList();
	FEFaceList(const FEFaceList& faceList);

	bool Create(FEMesh& mesh, FEElemElemList& ENL);

	int Faces() const;

	const FACE& operator [] (int i) const;
	const FACE& Face(int i) const;

	FEMesh* GetMesh();

	// Extract the surface only
	FEFaceList GetSurface() const;

	// build the neighbor list
	void BuildNeighbors();

protected:
	FEMesh*				m_mesh;
	std::vector<FACE>	m_faceList;
};

class FECORE_API FENodeFaceList
{
public:
	FENodeFaceList();

	bool Create(FEFaceList& FL);

	int Faces(int node) const;
	std::vector<int> FaceList(int node) const;

private:
	std::vector<std::vector<int> >	m_NFL;
};

class FECORE_API FEElementFaceList
{
public:
	FEElementFaceList();

	bool Create(FEElementList& elemList, FEFaceList& faceList);

	int Faces(int elem) const;
	std::vector<int> FaceList(int elem) const;

private:
	std::vector<std::vector<int> >	m_EFL;
};
