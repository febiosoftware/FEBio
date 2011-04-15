#include "stdafx.h"
#include "FEDomain.h"
#include "FEMesh.h"
#include "log.h"
#include "FESolidSolver.h"
#include "FEUncoupledMaterial.h"
#include "FETransverselyIsotropic.h"

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
// FESolidDomain
//-----------------------------------------------------------------------------
bool FESolidDomain::Initialize(FEM &fem)
{
	int i, j;
	FEMesh& m = fem.m_mesh;
	int N = m.Nodes();
	vector<int> tag; tag.assign(N, -1);

	int NE = Elements();
	int n = 0;
	m_Node.reserve(N);
	for (i=0; i<NE; ++i)
	{
		FESolidElement& e = Element(i);
		int ne = e.Nodes();
		for (j=0; j<ne; ++j)
		{
			int nj = e.m_node[j];
			if (tag[nj] == -1) 
			{
				tag[nj] = n++;
				m_Node.push_back(nj);
			}
		}
	}
	assert(m_Node.size() == n);
	return true;
}

//-----------------------------------------------------------------------------
void FESolidDomain::Serialize(DumpFile &ar)
{
	FEM& fem = *ar.GetFEM();
	if (ar.IsSaving())
	{
		ar << m_Node;

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
		ar >> m_Node;

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
FENode& FESolidDomain::Node(int i) 
{
	return m_pMesh->Node(m_Node[i]); 
}

//-----------------------------------------------------------------------------
void solve_3x3(double A[3][3], double b[3], double x[3])
{
	double D = A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[1][0]*A[2][1]*A[0][2] \
			 - A[1][1]*A[2][0]*A[0][2] - A[2][2]*A[1][0]*A[0][1] - A[0][0]*A[2][1]*A[1][2];

	assert(D != 0);

	double Ai[3][3];
	Ai[0][0] = A[1][1]*A[2][2] - A[2][1]*A[1][2];
	Ai[0][1] = A[2][1]*A[0][2] - A[0][1]*A[2][2];
	Ai[0][2] = A[0][1]*A[1][2] - A[1][1]*A[0][2];

	Ai[1][0] = A[2][0]*A[1][2] - A[1][0]*A[2][2];
	Ai[1][1] = A[0][0]*A[2][2] - A[2][0]*A[0][2];
	Ai[1][2] = A[1][0]*A[0][2] - A[0][0]*A[1][2];

	Ai[2][0] = A[1][0]*A[2][1] - A[2][0]*A[1][1];
	Ai[2][1] = A[2][0]*A[0][1] - A[0][0]*A[2][1];
	Ai[2][2] = A[0][0]*A[1][1] - A[0][1]*A[1][0];

	x[0] = (Ai[0][0]*b[0] + Ai[0][1]*b[1] + Ai[0][2]*b[2])/D;
	x[1] = (Ai[1][0]*b[0] + Ai[1][1]*b[1] + Ai[1][2]*b[2])/D;
	x[2] = (Ai[2][0]*b[0] + Ai[2][1]*b[1] + Ai[2][2]*b[2])/D;


#ifdef _DEBUG
	double r[3];
	r[0] = b[0] - (A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2]); 
	r[1] = b[1] - (A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2]); 
	r[2] = b[2] - (A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2]);

	double nr = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
#endif
}

//-----------------------------------------------------------------------------
//! This function finds the element in which point y lies and returns
//! the isoparametric coordinates in r if an element is found
//! (This has only been implemeneted for hexes!)
FESolidElement* FESolidDomain::FindElement(vec3d y, double r[3])
{
	int i, j;
	int NE = Elements();
	vec3d x[8];
	for (i=0; i<NE; ++i)
	{
		// get the next element
		FESolidElement& e = Element(i);
		assert(e.Type() == FE_HEX);

		// get the element nodal coordinates
		int neln = e.Nodes();
		for (j=0; j<neln; ++j) x[j] = m_pMesh->Node(e.m_node[j]).m_rt;

		// first, as a quick check, we see if y lies in the bounding box defined by x
		FE_BOUNDING_BOX box;
		box.r0 = box.r1 = x[0];
		for (j=1; j<neln; ++j) box += x[j];

		if (box.IsInside(y))
		{
			// If the point y lies inside the box, we apply a Newton method to find
			// the isoparametric coordinates r
			r[0] = r[1] = r[2] = 0;
			const double tol = 1e-5;
			double dr[3], norm;
			double H[8], G[8][3];
			do
			{
				H[0] = 0.125*(1 - r[0])*(1 - r[1])*(1 - r[2]);
				H[1] = 0.125*(1 + r[0])*(1 - r[1])*(1 - r[2]);
				H[2] = 0.125*(1 + r[0])*(1 + r[1])*(1 - r[2]);
				H[3] = 0.125*(1 - r[0])*(1 + r[1])*(1 - r[2]);
				H[4] = 0.125*(1 - r[0])*(1 - r[1])*(1 + r[2]);
				H[5] = 0.125*(1 + r[0])*(1 - r[1])*(1 + r[2]);
				H[6] = 0.125*(1 + r[0])*(1 + r[1])*(1 + r[2]);
				H[7] = 0.125*(1 - r[0])*(1 + r[1])*(1 + r[2]);

				G[0][0] = -0.125*(1 - r[1])*(1 - r[2]); G[0][1] = -0.125*(1 - r[0])*(1 - r[2]); G[0][2] = -0.125*(1 - r[0])*(1 - r[1]);
				G[1][0] =  0.125*(1 - r[1])*(1 - r[2]); G[1][1] = -0.125*(1 + r[0])*(1 - r[2]); G[1][2] = -0.125*(1 + r[0])*(1 - r[1]);
				G[2][0] =  0.125*(1 + r[1])*(1 - r[2]); G[2][1] =  0.125*(1 + r[0])*(1 - r[2]); G[2][2] = -0.125*(1 + r[0])*(1 + r[1]);
				G[3][0] = -0.125*(1 + r[1])*(1 - r[2]); G[3][1] =  0.125*(1 - r[0])*(1 - r[2]); G[3][2] = -0.125*(1 - r[0])*(1 + r[1]);
				G[4][0] = -0.125*(1 - r[1])*(1 + r[2]); G[4][1] = -0.125*(1 - r[0])*(1 + r[2]); G[4][2] =  0.125*(1 - r[0])*(1 - r[1]);
				G[5][0] =  0.125*(1 - r[1])*(1 + r[2]); G[5][1] = -0.125*(1 + r[0])*(1 + r[2]); G[5][2] =  0.125*(1 + r[0])*(1 - r[1]);
				G[6][0] =  0.125*(1 + r[1])*(1 + r[2]); G[6][1] =  0.125*(1 + r[0])*(1 + r[2]); G[6][2] =  0.125*(1 + r[0])*(1 + r[1]);
				G[7][0] = -0.125*(1 + r[1])*(1 + r[2]); G[7][1] =  0.125*(1 - r[0])*(1 + r[2]); G[7][2] =  0.125*(1 - r[0])*(1 + r[1]);

				double R[3] = {0}, A[3][3] = {0};
				for (j=0; j<8; ++j)
				{
					R[0] += x[j].x*H[j];
					R[1] += x[j].y*H[j];
					R[2] += x[j].z*H[j];

					A[0][0] -= x[j].x*G[j][0]; A[0][1] -= x[j].x*G[j][1]; A[0][2] -= x[j].x*G[j][2];
					A[1][0] -= x[j].y*G[j][0]; A[1][1] -= x[j].y*G[j][1]; A[1][2] -= x[j].y*G[j][2];
					A[2][0] -= x[j].z*G[j][0]; A[2][1] -= x[j].z*G[j][1]; A[2][2] -= x[j].z*G[j][2];
				}
				R[0] = y.x - R[0];
				R[1] = y.y - R[1];
				R[2] = y.z - R[2];

				solve_3x3(A, R, dr);
				r[0] -= dr[0];
				r[1] -= dr[1];
				r[2] -= dr[2];

				norm = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
			}
			while (norm > tol);

			// see if the point r lies inside the element
			const double eps = 1.0001;
			if ((r[0] >= -eps) && (r[0] <= eps) &&
				(r[1] >= -eps) && (r[1] <= eps) && 
				(r[2] >= -eps) && (r[2] <= eps)) return &e;
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------
// FEElasticSolidDomain
//-----------------------------------------------------------------------------
void FEElasticSolidDomain::Reset()
{
	for (int i=0; i<(int) m_Elem.size(); ++i) m_Elem[i].Init(true);
}

//-----------------------------------------------------------------------------
bool FEElasticSolidDomain::Initialize(FEM &fem)
{
	// initialize base class
	FESolidDomain::Initialize(fem);

	bool bmerr = false;

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
			clog.printbox("F A T A L   E R R O R", "A negative jacobian was detected at\n integration point %d of element %d.\nDid you use the right node numbering?", e.m_iel, e.m_ng);
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
								clog.printbox("ERROR", "Solid element %d was not assigned a fiber direction.", i+1);
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
// FE3FieldElasticSolidDomain
//-----------------------------------------------------------------------------
bool FE3FieldElasticSolidDomain::Initialize(FEM &fem)
{
	// make sure the domain material uses an uncoupled formulation
	if (dynamic_cast<FEUncoupledMaterial*>(m_pMat) == 0) return false;
	return FEElasticSolidDomain::Initialize(fem);
}

//-----------------------------------------------------------------------------
// FEShellDomain
//-----------------------------------------------------------------------------
bool FEShellDomain::Initialize(FEM &fem)
{
	int i, j;
	FEMesh& m = fem.m_mesh;
	int N = m.Nodes();
	vector<int> tag; tag.assign(N, -1);

	int NE = Elements();
	int n = 0;
	for (i=0; i<NE; ++i)
	{
		FEShellElement& e = Element(i);
		int ne = e.Nodes();
		for (j=0; j<ne; ++j)
		{
			int nj = e.m_node[j];
			if (tag[nj] == -1) tag[nj] = n++;
		}
	}
	m_Node.reserve(n);
	for (i=0; i<N; ++i) if (tag[i] >= 0) m_Node.push_back(i);
	assert(m_Node.size() == n);
	return true;
}

//-----------------------------------------------------------------------------
FENode& FEShellDomain::Node(int i) 
{
	return m_pMesh->Node(m_Node[i]); 
}

//-----------------------------------------------------------------------------
void FEShellDomain::Serialize(DumpFile &ar)
{
	FEM& fem = *ar.GetFEM();
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
// FEElasticShellDomain
//-----------------------------------------------------------------------------
void FEElasticShellDomain::Reset()
{
	for (int i=0; i<(int) m_Elem.size(); ++i) m_Elem[i].Init(true);
}

//-----------------------------------------------------------------------------
bool FEElasticShellDomain::Initialize(FEM& fem)
{
	// initialize base class
	FEShellDomain::Initialize(fem);

	bool bmerr = false;

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
							clog.printbox("ERROR", "Shell element %d was not assigned a fiber direction.", i+1);
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
// FETrussDomain
//-----------------------------------------------------------------------------
bool FETrussDomain::Initialize(FEM &fem)
{
	int i, j;
	FEMesh& m = fem.m_mesh;
	int N = m.Nodes();
	vector<int> tag; tag.assign(N, -1);

	int NE = Elements();
	int n = 0;
	for (i=0; i<NE; ++i)
	{
		FETrussElement& e = Element(i);
		int ne = e.Nodes();
		for (j=0; j<ne; ++j)
		{
			int nj = e.m_node[j];
			if (tag[nj] == -1) tag[nj] = n++;
		}
	}
	m_Node.reserve(n);
	for (i=0; i<N; ++i) if (tag[i] >= 0) m_Node.push_back(i);
	assert(m_Node.size() == n);
	return true;
}

//-----------------------------------------------------------------------------
FENode& FETrussDomain::Node(int i) 
{
	return m_pMesh->Node(m_Node[i]); 
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
// FEDiscreteDomain
//-----------------------------------------------------------------------------

bool FEDiscreteDomain::Initialize(FEM &fem)
{
	int i, j;
	FEMesh& m = fem.m_mesh;
	int N = m.Nodes();
	vector<int> tag; tag.assign(N, -1);

	int NE = Elements();
	int n = 0;
	for (i=0; i<NE; ++i)
	{
		FEDiscreteElement& e = m_Elem[i];
		int ne = e.Nodes();
		for (j=0; j<ne; ++j)
		{
			int nj = e.m_node[j];
			if (tag[nj] == -1) tag[nj] = n++;
		}
	}
	m_Node.reserve(n);
	for (i=0; i<N; ++i) if (tag[i] >= 0) m_Node.push_back(i);
	assert(m_Node.size() == n);
	return true;
}

//-----------------------------------------------------------------------------
void FEDiscreteDomain::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_Node;
		int nel = (int) m_Elem.size();
		for (int i=0; i<nel; ++i)
		{
			FEDiscreteElement& el = m_Elem[i];
			int nmat = el.GetMatID();
			ar << (int) el.Type();
			
			ar << nmat;
			ar << el.m_nrigid;
			ar << el.m_nID;
			ar << el.m_node;

			for (int j=0; j<el.GaussPoints(); ++j) el.m_State[j]->Serialize(ar);
		}
	}
	else
	{
		FEM& fem = *ar.GetFEM();
		ar >> m_Node;
		int n, mat;
		for (size_t i=0; i<m_Elem.size(); ++i)
		{
			FEDiscreteElement& el = m_Elem[i];
			ar >> n;

			el.SetType(n);

			ar >> mat; el.SetMatID(mat);
			ar >> el.m_nrigid;
			ar >> el.m_nID;
			ar >> el.m_node;

			for (int j=0; j<el.GaussPoints(); ++j)
			{
				el.SetMaterialPointData(fem.GetMaterial(el.GetMatID())->CreateMaterialPointData(), j);
				el.m_State[j]->Serialize(ar);
			}
		}
	}
}

//-----------------------------------------------------------------------------
FENode& FEDiscreteDomain::Node(int i) 
{
	return m_pMesh->Node(m_Node[i]); 
}

//-----------------------------------------------------------------------------
void FEDiscreteDomain::UnpackElement(FEElement &el, unsigned int nflag)
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
