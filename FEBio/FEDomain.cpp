#include "stdafx.h"
#include "FEDomain.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
// FESolidDomain
//-----------------------------------------------------------------------------
FEElement* FESolidDomain::FindElementFromID(int nid)
{
	for (int i=0; i<m_Elem.size(); ++i)
		if (m_Elem[i].m_nID == nid) return &m_Elem[i];

	return 0;
}

//-----------------------------------------------------------------------------
void FESolidDomain::Reset()
{
	for (int i=0; i<(int) m_Elem.size(); ++i) m_Elem[i].Init(true);
}

//-----------------------------------------------------------------------------
// FEShellDomain
//-----------------------------------------------------------------------------
FEElement* FEShellDomain::FindElementFromID(int nid)
{
	for (int i=0; i<m_Elem.size(); ++i)
		if (m_Elem[i].m_nID == nid) return &m_Elem[i];

	return 0;
}

//-----------------------------------------------------------------------------
void FEShellDomain::Reset()
{
	for (int i=0; i<(int) m_Elem.size(); ++i) m_Elem[i].Init(true);
}

//-----------------------------------------------------------------------------
// FETrussDomain
//-----------------------------------------------------------------------------
FEElement* FETrussDomain::FindElementFromID(int nid)
{
	for (int i=0; i<m_Elem.size(); ++i)
		if (m_Elem[i].m_nID == nid) return &m_Elem[i];

	return 0;
}

//-----------------------------------------------------------------------------
void FETrussDomain::Reset()
{
	for (int i=0; i<(int) m_Elem.size(); ++i) m_Elem[i].Init(true);
}
