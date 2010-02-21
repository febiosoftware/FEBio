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
//-----------------------------------------------------------------------------

void FETrussDomain::UnpackElement(FETrussElement &el, unsigned int nflag)
{
	int i, n;

	vec3d* rt = el.rt();
	vec3d* r0 = el.r0();
	vec3d* vt = el.vt();
	double* pt = el.pt();

	int N = el.Nodes();
	int* lm = el.LM();

	for (i=0; i<N; ++i)
	{
		n = el.m_node[i];
		FENode& node = m_pMesh->Node(n);

		int* id = node.m_ID;

		// first the displacement dofs
		lm[3*i  ] = id[0];
		lm[3*i+1] = id[1];
		lm[3*i+2] = id[2];

		// now the pressure dofs
		lm[3*N+i] = id[6];

		// rigid rotational dofs
		lm[4*N + 3*i  ] = id[7];
		lm[4*N + 3*i+1] = id[8];
		lm[4*N + 3*i+2] = id[9];

		// fill the rest with -1
		lm[7*N + 3*i  ] = -1;
		lm[7*N + 3*i+1] = -1;
		lm[7*N + 3*i+2] = -1;

		lm[10*N + i] = id[10];
	}

	// copy nodal data to element arrays
	for (i=0; i<N; ++i)
	{
		n = el.m_node[i];

		FENode& node = m_pMesh->Node(n);

		// initial coordinates (= material coordinates)
		r0[i] = node.m_r0;

		// current coordinates (= spatial coordinates)
		rt[i] = node.m_rt;

		// current nodal pressures
		pt[i] = node.m_pt;

		// current nodal velocities
		vt[i] = node.m_vt;
	}

	// unpack the traits data
	el.UnpackTraitsData(nflag);
}
