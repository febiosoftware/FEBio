#pragma once
#include <vector>
#include "fecore_api.h"

class FEMesh;
class FEElementList;

class FECORE_API FEFaceList
{
public:
	struct FACE
	{
		int ntype;		// 3 = triangle, 4 = quad
		int	node[4];

		bool IsEqual(int* n) const;
	};

public:
	FEFaceList();

	bool Create(FEMesh* mesh);

	int Faces() const;

	const FACE& Face(int i) const;

	FEMesh* GetMesh();

protected:
	FEMesh*				m_mesh;
	std::vector<FACE>	m_faceList;
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
