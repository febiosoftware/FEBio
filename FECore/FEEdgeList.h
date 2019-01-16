#pragma once
#include <vector>
#include "fecore_api.h"

class FEMesh;
class FEElementList;

class FECORE_API FEEdgeList
{
public:
	struct EDGE
	{
		int ntype;		// 2 = linear, 3 = cubic
		int	node[3];
	};

public:
	FEEdgeList();

	bool Create(FEMesh* mesh);

	int Edges() const;

	const EDGE& operator[] (int i);
	const EDGE& Edge(int i) const;

	FEMesh* GetMesh();

private:
	FEMesh*				m_mesh;
	std::vector<EDGE>	m_edgeList;
};

class FECORE_API FEElementEdgeList
{
public:
	FEElementEdgeList();

	bool Create(FEElementList& elemList, FEEdgeList& edgeList);

	int Edges(int elem) const;
	std::vector<int> EdgeList(int elem) const;

private:
	std::vector<std::vector<int> >	m_EEL;
};
