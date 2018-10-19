#pragma once
#include "FEMeshPartition.h"
#include <vector>

class FENodeSet;

//-----------------------------------------------------------------------------
// This class represents an edge of a domain.
class FECORE_API FEEdge : public FEMeshPartition
{
public:
	//! constructor
	FEEdge(FEModel* fem);

	//! destructor
	virtual ~FEEdge();

	//! initialize edge data structure
	virtual bool Init() override;

	//! creates edge
	void Create(int nelems, int elemType = -1) override;

	//! extract node set
	FENodeSet GetNodeSet();

public:

	//! return number of edge elements
	int Elements() const override { return (int)m_Elem.size(); }

	//! return an element of the edge
	FELineElement& Element(int i) { return m_Elem[i]; }

	//! returns reference to element
	FEElement& ElementRef(int n) override { return m_Elem[n]; }
	const FEElement& ElementRef(int n) const override { return m_Elem[n]; }

protected:
	std::vector<FELineElement>	m_Elem;
};
