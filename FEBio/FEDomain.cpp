#include "stdafx.h"
#include "FEDomain.h"
#include "FEMesh.h"
#include "log.h"
#include "FESolidSolver.h"

//-----------------------------------------------------------------------------
// FEDomain
//-----------------------------------------------------------------------------
FEElement* FEDomain::FindElementFromID(int nid)
{
	for (int i=0; i<Elements(); ++i)
	{
		FEElement& el = ElementRef(i);
		if (el.m_nID == nid) return &el;
	}

	return 0;
}

//-----------------------------------------------------------------------------

void FEDomain::InitMaterialPointData()
{
	if (m_pMat == 0) return;

	for (int i=0; i<Elements(); ++i)
	{
		FEElement& el = ElementRef(i);
		for (int k=0; k<el.GaussPoints(); ++k) el.SetMaterialPointData(m_pMat->CreateMaterialPointData(), k);
	}
}

//-----------------------------------------------------------------------------
// FEElasticSolidDomain
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void FEElasticSolidDomain::Reset()
{
	for (int i=0; i<(int) m_Elem.size(); ++i) m_Elem[i].Init(true);
}

//-----------------------------------------------------------------------------

bool FEElasticSolidDomain::Initialize(FEM &fem)
{
	bool bmerr = false;

	// get the logfile
	Logfile& log = GetLogfile();

	for (size_t i=0; i<m_Elem.size(); ++i)
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

void FEElasticSolidDomain::InitElements()
{
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		int n = el.GaussPoints();
		for (int j=0; j<n; ++j) el.m_State[j]->Init(false);
	}
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain::Serialize(FEM& fem, Archive &ar)
{
	if (ar.IsSaving())
	{
		for (size_t i=0; i<m_Elem.size(); ++i)
		{
			FESolidElement& el = m_Elem[i];
			int nmat = el.GetMatID();
			ar << el.Type();
			
			ar << nmat;
			ar << el.m_nrigid;
			ar << el.m_nID;
			ar << el.m_node;

			ar << el.m_eJ;
			ar << el.m_ep;
			ar << el.m_Lk;

			for (int j=0; j<el.GaussPoints(); ++j) el.m_State[j]->Serialize(ar);
		}
	}
	else
	{
		int n, mat;
		for (size_t i=0; i<m_Elem.size(); ++i)
		{
			FESolidElement& el = m_Elem[i];
			ar >> n;

			el.SetType(n);

			ar >> mat; el.SetMatID(mat);
			ar >> el.m_nrigid;
			ar >> el.m_nID;
			ar >> el.m_node;

			ar >> el.m_eJ;
			ar >> el.m_ep;
			ar >> el.m_Lk;

			for (int j=0; j<el.GaussPoints(); ++j)
			{
				el.SetMaterialPointData(fem.GetMaterial(el.GetMatID())->CreateMaterialPointData(), j);
				el.m_State[j]->Serialize(ar);
			}
		}
	}
}

//-----------------------------------------------------------------------------
// FEElasticShellDomain
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void FEElasticShellDomain::Reset()
{
	for (int i=0; i<(int) m_Elem.size(); ++i) m_Elem[i].Init(true);
}

//-----------------------------------------------------------------------------
bool FEElasticShellDomain::Initialize(FEM& fem)
{
	bool bmerr = false;

	// get the logfile
	Logfile& log = GetLogfile();

	for (size_t i=0; i<m_Elem.size(); ++i)
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
void FEElasticShellDomain::Serialize(FEM& fem, Archive &ar)
{
	if (ar.IsSaving())
	{
		for (size_t i=0; i<m_Elem.size(); ++i)
		{
			FEShellElement& el = m_Elem[i];
			ar << el.Type();

			ar << el.m_eJ;
			ar << el.m_ep;

			ar << el.GetMatID();
			ar << el.m_nrigid;
			ar << el.m_nID;
			ar << el.m_node;

			ar << el.m_h0;
			ar << el.m_Lk;

			for (int j=0; j<el.GaussPoints(); ++j) el.m_State[j]->Serialize(ar);
		}
	}
	else
	{
		int n, mat;

		for (size_t i=0; i<m_Elem.size(); ++i)
		{
			FEShellElement& el = m_Elem[i];
			ar >> n;

			el.SetType(n);

			ar >> el.m_eJ;
			ar >> el.m_ep;

			ar >> mat; el.SetMatID(mat);
			ar >> el.m_nrigid;
			ar >> el.m_nID;
			ar >> el.m_node;

			ar >> el.m_h0;
			ar >> el.m_Lk;

			for (int j=0; j<el.GaussPoints(); ++j)
			{
				el.SetMaterialPointData(fem.GetMaterial(el.GetMatID())->CreateMaterialPointData(), j);
				el.m_State[j]->Serialize(ar);
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FEElasticShellDomain::InitElements()
{
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FEShellElement& el = m_Elem[i];
		int n = el.GaussPoints();
		for (int j=0; j<n; ++j) el.m_State[j]->Init(false);
	}
}

//-----------------------------------------------------------------------------
// FEElasticTrussDomain
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void FEElasticTrussDomain::Reset()
{
	for (int i=0; i<(int) m_Elem.size(); ++i) m_Elem[i].Init(true);
}
//-----------------------------------------------------------------------------

void FEElasticTrussDomain::UnpackElement(FEElement &el, unsigned int nflag)
{
	int i, n;

	vec3d* rt = el.rt();
	vec3d* r0 = el.r0();
	vec3d* vt = el.vt();
	double* pt = el.pt();

	int N = el.Nodes();
	vector<int>& lm = el.LM();

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

//-----------------------------------------------------------------------------

void FEElasticTrussDomain::InitElements()
{
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FETrussElement& el = m_Elem[i];
		el.m_State[0]->Init(false);
	}
}

//-----------------------------------------------------------------------------
// FEHeatSolidDomain
//-----------------------------------------------------------------------------
