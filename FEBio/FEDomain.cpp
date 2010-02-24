#include "stdafx.h"
#include "FEDomain.h"
#include "FEMesh.h"
#include "log.h"
#include "FESolidSolver.h"

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

bool FESolidDomain::Init(FEM &fem)
{
	bool bmerr = false;

	// get the logfile
	Logfile& log = GetLogfile();

	for (int i=0; i<m_Elem.size(); ++i)
	{
		// unpack element data
		FESolidElement& el = m_Elem[i];

		try
		{
			UnpackElement(el);
		}
		catch (NegativeJacobian e)
		{
			log.printbox("F A T A L   E R R O R", "A negative jacobian was detected at\n integration point %d of element %d.\nDid you use the right node numbering?", e.m_iel, e.m_ng);
			return false;
		}

		if (dynamic_cast<FESolidSolver*>(fem.m_pStep->m_psolver))
		{
			// get the elements material
			FEElasticMaterial* pme = fem.GetElasticMaterial(el.GetMatID());

			// set the local element coordinates
			if (pme)
			{
				if (pme->m_pmap)
				{
					for (int n=0; n<el.GaussPoints(); ++n)
					{
						FEElasticMaterialPoint& pt = *el.m_State[n]->ExtractData<FEElasticMaterialPoint>();
						pt.Q = pme->m_pmap->LocalElementCoord(el, n);
					}
				}
				else
				{
					if (fem.GetDebugFlag())
					{
						// If we get here, then the element has a user-defined fiber axis
						// we should check to see if it has indeed been specified.
						// TODO: This assumes that pt.Q will not get intialized to
						//		 a valid value. I should find another way for checking since I
						//		 would like pt.Q always to be initialized to a decent value.
						if (dynamic_cast<FETransverselyIsotropic*>(pme))
						{
							FEElasticMaterialPoint& pt = *el.m_State[0]->ExtractData<FEElasticMaterialPoint>();
							mat3d& m = pt.Q;
							if (fabs(m.det() - 1) > 1e-7)
							{
								// this element did not get specified a user-defined fiber direction
								log.printbox("ERROR", "Solid element %d was not assigned a fiber direction.", i+1);
								bmerr = true;
							}
						}
					}
				}
			}
		}
	}

	return (bmerr == false);
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
bool FEShellDomain::Init(FEM& fem)
{
	bool bmerr = false;

	// get the logfile
	Logfile& log = GetLogfile();

	for (int i=0; i<m_Elem.size(); ++i)
	{
		// unpack element data
		FEShellElement& el = m_Elem[i];

		// get the elements material
		FEElasticMaterial* pme = fem.GetElasticMaterial(el.GetMatID());

		// set the local element coordinates
		if (pme)
		{
			if (pme->m_pmap)
			{
				for (int n=0; n<el.GaussPoints(); ++n)
				{
					FEElasticMaterialPoint& pt = *el.m_State[n]->ExtractData<FEElasticMaterialPoint>();
					pt.Q = pme->m_pmap->LocalElementCoord(el, n);
				}
			}
			else
			{
				if (fem.GetDebugFlag())
				{
					// If we get here, then the element has a user-defined fiber direction
					// we should check to see if it has indeed been specified.
					// TODO: This assumes that pt.Q will not get intialized to
					//		 a valid value. I should find another way for checking since I
					//		 would like pt.Q always to be initialized to a decent value.
					if (dynamic_cast<FETransverselyIsotropic*>(pme))
					{
						FEElasticMaterialPoint& pt = *el.m_State[0]->ExtractData<FEElasticMaterialPoint>();
						mat3d& m = pt.Q;
						if (fabs(m.det() - 1) > 1e-7)
						{
							// this element did not get specified a user-defined fiber direction
							log.printbox("ERROR", "Shell element %d was not assigned a fiber direction.", i+1);
							bmerr = true;
						}
					}
				}
			}
		}
	}
	return (bmerr == false);
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
