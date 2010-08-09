// FENodeElemList.h: interface for the FENodeElemList class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FENODEELEMLIST_H__774D7DB9_D0A9_4F6F_AE4D_191126B5D3F8__INCLUDED_)
#define AFX_FENODEELEMLIST_H__774D7DB9_D0A9_4F6F_AE4D_191126B5D3F8__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/vector.h"

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

	int Valence(int n) { return m_nval[n]; }
	FEElement** ElementList(int n) { return &m_eref[0] + m_pn[n]; }

	int Size() { return (int) m_nval.size(); }

protected:
	vector<int>			m_nval;	// nodal valences
	vector<FEElement*>	m_eref;	// element pointers
	vector<int>			m_pn;	// start index into the eref array
};

#endif // !defined(AFX_FENODEELEMLIST_H__774D7DB9_D0A9_4F6F_AE4D_191126B5D3F8__INCLUDED_)
