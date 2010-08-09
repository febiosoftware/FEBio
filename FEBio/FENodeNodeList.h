// FENodeNodeList.h: interface for the FENodeNodeList class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FENODENODELIST_H__59D213DB_78A1_4318_9E5D_585E70BCB36D__INCLUDED_)
#define AFX_FENODENODELIST_H__59D213DB_78A1_4318_9E5D_585E70BCB36D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/vector.h"

class FEMesh;
class FEDomain;

//-----------------------------------------------------------------------------
//! The FENodeNodeList class is a utility class that determines for each node 
//! the adjacent nodes

//! This class analyzes a mesh and finds for each node all nodes that are 
//! adjacent to this node

class FENodeNodeList  
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

	int Size() { return m_nval.size(); }

	int Valence(int i) { return m_nval[i]; }
	int* NodeList(int i) { return &m_nref[0] + m_pn[i]; }

	void Sort();

protected:
	vector<int>	m_nval;	// nodal valences
	vector<int>	m_nref;	// adjacent nodes indices
	vector<int>	m_pn;	// start index into the nref array

	static FENodeNodeList*	m_pthis;
	static int compare(const void* e1, const void* e2);
};

#endif // !defined(AFX_FENODENODELIST_H__59D213DB_78A1_4318_9E5D_585E70BCB36D__INCLUDED_)
