#pragma once
#include "FEDomain.h"
#include <vector>

//-----------------------------------------------------------------------------
// This class represents an edge of a domain.
class FEEdge : public FEDomain
{
public:
	//! constructor
	FEEdge(FEMesh* pm);

	//! destructor
	virtual ~FEEdge();

	//! initialize edge data structure
	virtual bool Init();

	//! creates surface
	void create(int n);

	//! serialization
	void Serialize(DumpStream& ar);

public:

	//! return number of edge elements
	int Elements() { return m_Elem.size(); }

	//! return an element of the edge
	FELineElement& Element(int i) { return m_Elem[i]; }

	//! returns reference to element
	FEElement& ElementRef(int n) { return m_Elem[n]; }

protected:
	std::vector<FELineElement>	m_Elem;
};
