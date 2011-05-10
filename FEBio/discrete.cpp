#include "stdafx.h"
#include "FESolidSolver.h"
#include "FEDiscreteDomain.h"
#include "FEDiscreteMaterial.h"
#include "log.h"


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

//-----------------------------------------------------------------------------
//! Calculates the forces due to discrete elements (i.e. springs)

void FEDiscreteDomain::Residual(FESolidSolver* psolver, vector<double>& R)
{
	FEM& fem = psolver->m_fem;
	FEMesh& mesh = *m_pMesh;

	vector<double> fe(6);
	vec3d u1, u2;

	vector<int> en(2), lm(6);

	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		// get the discrete element
		FEDiscreteElement& el = m_Elem[i];

		// get the nodes
		FENode& n1 = mesh.Node(el.m_node[0]);
		FENode& n2 = mesh.Node(el.m_node[1]);

		// get the nodal positions
		vec3d& r01 = n1.m_r0;
		vec3d& r02 = n2.m_r0;
		vec3d& rt1 = n1.m_rt;
		vec3d& rt2 = n2.m_rt;

		vec3d e = rt2 - rt1; e.unit();

		// calculate spring lengths
		double L0 = (r02 - r01).norm();
		double Lt = (rt2 - rt1).norm();
		double DL = Lt - L0;
		
		// evaluate the spring stiffness
		FEDiscreteMaterial* pm = dynamic_cast<FEDiscreteMaterial*>(fem.GetMaterial(el.GetMatID()));
		assert(pm);
		double F = pm->force(DL);

		// set up the force vector
		fe[0] =   F*e.x;
		fe[1] =   F*e.y;
		fe[2] =   F*e.z;
		fe[3] =  -F*e.x;
		fe[4] =  -F*e.y;
		fe[5] =  -F*e.z;

		// setup the node vector
		en[0] = el.m_node[0];
		en[1] = el.m_node[1];

		// set up the LM vector
		lm[0] = n1.m_ID[0];
		lm[1] = n1.m_ID[1];
		lm[2] = n1.m_ID[2];
		lm[3] = n2.m_ID[0];
		lm[4] = n2.m_ID[1];
		lm[5] = n2.m_ID[2];

		// assemble element
		psolver->AssembleResidual(en, lm, fe, R);
	}
}

//-----------------------------------------------------------------------------
//! Calculates the discrete element stiffness

void FEDiscreteDomain::StiffnessMatrix(FESolidSolver* psolver)
{
	FEM& fem = psolver->m_fem;
	FEMesh& mesh = *m_pMesh;

	matrix ke(6,6);
	ke.zero();

	vector<int> en(2), lm(6);

	// loop over all discrete elements
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		// get the discrete element
		FEDiscreteElement& el = m_Elem[i];

		// get the nodes of the element
		FENode& n1 = mesh.Node(el.m_node[0]);
		FENode& n2 = mesh.Node(el.m_node[1]);

		// get the nodal positions
		vec3d& r01 = n1.m_r0;
		vec3d& r02 = n2.m_r0;
		vec3d& rt1 = n1.m_rt;
		vec3d& rt2 = n2.m_rt;

		vec3d e = rt2 - rt1; 

		// calculate nodal displacements
		vec3d u1 = rt1 - r01;
		vec3d u2 = rt2 - r02;

		// calculate spring lengths
		double L0 = (r02 - r01).norm();
		double Lt = (rt2 - rt1).norm();
		double DL = Lt - L0;

		// evaluate the stiffness
		FEDiscreteMaterial* pm = dynamic_cast<FEDiscreteMaterial*>(fem.GetMaterial(el.GetMatID()));
		assert(pm);
		double F = pm->force(DL);
		double E = pm->stiffness(DL);

		if (Lt == 0) { F = 0; Lt = 1; e = vec3d(1,1,1); }

		double A[3][3] = {0};
		A[0][0] = (E - F/Lt)*e.x*e.x + F/Lt;
		A[1][1] = (E - F/Lt)*e.y*e.y + F/Lt;
		A[2][2] = (E - F/Lt)*e.z*e.z + F/Lt;

		A[0][1] = A[1][0] = (E - F/Lt)*e.x*e.y;
		A[1][2] = A[2][1] = (E - F/Lt)*e.y*e.z;
		A[0][2] = A[2][0] = (E - F/Lt)*e.x*e.z;

		ke[0][0] = A[0][0]; ke[0][1] = A[0][1]; ke[0][2] = A[0][2];
		ke[1][0] = A[1][0]; ke[1][1] = A[1][1]; ke[1][2] = A[1][2];
		ke[2][0] = A[2][0]; ke[2][1] = A[2][1]; ke[2][2] = A[2][2];

		ke[0][3] = -A[0][0]; ke[0][4] = -A[0][1]; ke[0][5] = -A[0][2];
		ke[1][3] = -A[1][0]; ke[1][4] = -A[1][1]; ke[1][5] = -A[1][2];
		ke[2][3] = -A[2][0]; ke[2][4] = -A[2][1]; ke[2][5] = -A[2][2];

		ke[3][0] = -A[0][0]; ke[3][1] = -A[0][1]; ke[3][2] = -A[0][2];
		ke[4][0] = -A[1][0]; ke[4][1] = -A[1][1]; ke[4][2] = -A[1][2];
		ke[5][0] = -A[2][0]; ke[5][1] = -A[2][1]; ke[5][2] = -A[2][2];

		ke[3][3] = A[0][0]; ke[3][4] = A[0][1]; ke[3][5] = A[0][2];
		ke[4][3] = A[1][0]; ke[4][4] = A[1][1]; ke[4][5] = A[1][2];
		ke[5][3] = A[2][0]; ke[5][4] = A[2][1]; ke[5][5] = A[2][2];

		// setup the node vector
		en[0] = el.m_node[0];
		en[1] = el.m_node[1];

		// set up the LM vector
		lm[0] = n1.m_ID[0];
		lm[1] = n1.m_ID[1];
		lm[2] = n1.m_ID[2];
		lm[3] = n2.m_ID[0];
		lm[4] = n2.m_ID[1];
		lm[5] = n2.m_ID[2];

		// assemble the element into the global system
		psolver->AssembleStiffness(en, lm, ke);
	}
}
