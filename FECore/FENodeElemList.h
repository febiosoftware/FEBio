// FENodeElemList.h: interface for the FENodeElemList class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FENODEELEMLIST_H__774D7DB9_D0A9_4F6F_AE4D_191126B5D3F8__INCLUDED_)
#define AFX_FENODEELEMLIST_H__774D7DB9_D0A9_4F6F_AE4D_191126B5D3F8__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "DumpStream.h"
#include <vector>

class FESurface;
class FEMesh;
class FEElement;
class FEDomain;

//-----------------------------------------------------------------------------
//! The FENodeElemList class is a utility class that determines for each node 
//! to which element it belongs.

//! This class analyzes a mesh and finds for each node all elements that have
//! this node

class FENodeElemList  
{
public:
	FENodeElemList(){}
	virtual ~FENodeElemList(){}

	//! build the node-element list for a surface
	void Create(FESurface& s);

	//! build the node-selement list for a mesh
	void Create(FEMesh& mesh);

	//! build the node-element list for a domain
	void Create(FEDomain& dom);

	//! serialize data to/from dump file
	void Serialize(DumpStream& ar);

	int MaxValence();
	int Valence(int n) { return m_nval[n]; }
	FEElement** ElementList(int n) { return &m_eref[0] + m_pn[n]; }
	int* ElementIndexList(int n) { return &m_iref[0] + m_pn[n]; }

	int Size() { return (int) m_nval.size(); }

protected:
	std::vector<int>			m_nval;	// nodal valences
	std::vector<FEElement*>		m_eref;	// element pointers
	std::vector<int>			m_iref;	// element indices
	std::vector<int>			m_pn;	// start index into the eref array
};

//-----------------------------------------------------------------------------
//! Like the FEElemElemList, but can create multiple levels
class FENodeElemTree
{
public:
	FENodeElemTree() {}
	virtual ~FENodeElemTree() {}

	void Create(FESurface* ps, int k = 0);

	int Valence(int n) { return (int) m_nel[n].size(); }

	FEElement** ElementList(int n) { return &(m_nel[n][0]);}

	bool empty() { return m_nel.empty(); }

protected:
	std::vector< std::vector<FEElement*> >	m_nel;
};

#endif // !defined(AFX_FENODEELEMLIST_H__774D7DB9_D0A9_4F6F_AE4D_191126B5D3F8__INCLUDED_)
