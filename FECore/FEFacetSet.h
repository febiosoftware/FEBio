#pragma once
#include "fecore_api.h"
#include "FEElement.h"
#include "FENodeSet.h"
#include <vector>
#include <string>

//-----------------------------------------------------------------------------
//! This class defines a set of facets. This can be used in the creation of
//! surfaces.
class FECORE_API FEFacetSet
{
public:
	struct FACET
	{
		int	node[FEElement::MAX_NODES];
		int	ntype;	//	3=tri3, 4=quad4, 6=tri6, 7=tri7, 8=quad8
	};

public:
	FEFacetSet(FEMesh* mesh);

	void SetName(const std::string& name);
	const std::string& GetName() const;

	void Create(int n);

	int Faces() const { return (int)m_Face.size(); }
	FACET& Face(int i);
	const FACET& Face(int i) const;

	void Add(FEFacetSet* pf);

	FENodeSet GetNodeSet();

	void Serialize(DumpStream& ar);

	const FEMesh* GetMesh() const { return m_mesh; }

private:
	std::string			m_name;
	std::vector<FACET>	m_Face;

private:
	FEMesh*			m_mesh;
};

