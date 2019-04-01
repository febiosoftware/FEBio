#pragma once
#include <vector>
#include "fecore_api.h"

class FEMesh;
class FEDomain;

//-----------------------------------------------------------------------------
//! The FENodeNodeList class is a utility class that determines for each node 
//! the adjacent nodes

//! This class analyzes a mesh and finds for each node all nodes that are 
//! adjacent to this node

class FECORE_API FENodeNodeList
{
public:
	//! default constructor
	FENodeNodeList();

	//! desctructor
	virtual ~FENodeNodeList();

	//! create the node-node list for a mesh
	void Create(FEMesh& mesh);

	//! create the node-node list for a domain
	void Create(FEDomain& dom);

	int Size() const { return (int) m_nval.size(); }

	int Valence(int i) { return m_nval[i]; }
	int* NodeList(int i) { return &m_nref[0] + m_pn[i]; }

	void Sort();

protected:
	std::vector<int>	m_nval;	// nodal valences
	std::vector<int>	m_nref;	// adjacent nodes indices
	std::vector<int>	m_pn;	// start index into the nref array

	static FENodeNodeList*	m_pthis;
	static int compare(const void* e1, const void* e2);
};
