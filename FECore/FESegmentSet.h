#pragma once
#include "fecore_api.h"
#include "FEElement.h"
#include <vector>
#include <string>

//-----------------------------------------------------------------------------
//! This class defines a set of segments. This can be used in the creation of edges.
class FECORE_API FESegmentSet
{
public:
	struct SEGMENT
	{
		int	node[FEElement::MAX_NODES];
		int	ntype;	//	2=line2

		void Serialize(DumpStream& ar);
	};

public:
	FESegmentSet(FEMesh* pm);

	void SetName(const std::string& name);
	const std::string& GetName() const;

	void Create(int n);

	int Segments() { return (int)m_Seg.size(); }
	SEGMENT& Segment(int i);

	void Serialize(DumpStream& ar);

private:
	vector<SEGMENT>	m_Seg;
	std::string		m_name;
	FEMesh*			m_mesh;
};
