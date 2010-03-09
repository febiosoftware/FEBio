#include "stdafx.h"
#include "ut4.h"
#include "FEMesh.h"
#include "FESolidSolver.h"

//-----------------------------------------------------------------------------
//! Constructor for the UT4Domain
FEUT4Domain::FEUT4Domain(FEMesh*pm) : FEElasticSolidDomain(pm)
{
	// overwrite the domain type
	m_ntype = FE_UT4_DOMAIN; 

	// set the default stabilization factor
	m_alpha = 0.05; 
}

//-----------------------------------------------------------------------------
//! Initialization for the UT4Domain.
//! Note that we first initialize the base class before initializing the domain.
bool FEUT4Domain::Initialize(FEM& fem)
{
	// first call the base class
	if (FEElasticSolidDomain::Initialize(fem) == false) return false;

	// next, we need to identify all the nodes that belong to this domain
	// we do this by looping over all the elements and tagging the nodes
	int NN = m_pMesh->Nodes();
	m_tag.assign(NN, 0);

	int i, NE = Elements();
	for (i=0; i<NE; ++i)
	{
		FESolidElement& el = m_Elem[i];
		assert(el.Type() == FE_TET);
		for (int j=0; j<4; ++j) m_tag[el.m_node[j]] = 1;
	}

	// count the tags
	int N = 0;
	for (i=0; i<NN; ++i) if (m_tag[i] == 1) ++N;
	if (N < 4) return false;

	// allocate node structure
	m_Node.resize(N);

	// initialize nodes
	N = 0;
	for (i=0; i<NN; ++i)
		if (m_tag[i] == 1)
		{
			UT4NODE& n = m_Node[N];
			m_tag[i] = N++;
			n.inode = i;
			n.Vi = n.vi = 0.0;
		}
		else m_tag[i] = -1;

	// calculate the reference nodal volume
	// we do this here since this volume never changes
	double Ve;
	for (i=0; i<NE; ++i)
	{
		FESolidElement& el = m_Elem[i];
		UnpackElement(el);

		// calculate the initial volume
		Ve = TetVolume(el.r0());

		// now assign one-quart to each node
		for (int j=0; j<4; ++j) m_Node[ m_tag[el.m_node[j]]].Vi += 0.25*Ve;
	}


	return true;
}

//-----------------------------------------------------------------------------
//! Update the nodal and element stresses
void FEUT4Domain::UpdateStresses(FEM &fem)
{
	// updating the element stresses is easy, since we only
	// need to call the base class
	FEElasticSolidDomain::UpdateStresses(fem);

	// next we update the nodal data
	// first, let's calculate the new nodal volumes
	int i, NE = Elements();
	for (i=0; i<(int) m_Node.size(); ++i) m_Node[i].vi = 0;
	double ve;
	for (i=0; i<NE; ++i)
	{
		FESolidElement& el = m_Elem[i];
		UnpackElement(el);

		// calculate the current volume
		ve = TetVolume(el.rt());

		// now assign one-quart to each node
		for (int j=0; j<4; ++j) m_Node[ m_tag[el.m_node[j]]].vi += 0.25*ve;
	}
}

//-----------------------------------------------------------------------------
//! Calculate volume for tets using the determinant rule
double FEUT4Domain::TetVolume(vec3d* r)
{
	double A[4][4];
	A[0][0] = r[0].x; A[0][1] = r[1].x; A[0][2] = r[2].x; A[0][3] = r[3].x;
	A[1][0] = r[0].y; A[1][1] = r[1].y; A[1][2] = r[2].y; A[1][3] = r[3].y;
	A[2][0] = r[0].z; A[2][1] = r[1].z; A[2][2] = r[2].z; A[2][3] = r[3].z;
	A[3][0] =    1.0; A[3][1] =    1.0; A[3][2] =    1.0; A[3][3] =    1.0;

	double V = 0;
	V += mat3d(A[0][1], A[0][2], A[0][3], A[1][1], A[1][2], A[1][3], A[2][1], A[2][2], A[2][3]).det();
	V -= mat3d(A[0][0], A[0][2], A[0][3], A[1][0], A[1][2], A[1][3], A[2][0], A[2][2], A[2][3]).det();
	V += mat3d(A[0][0], A[0][1], A[0][3], A[1][0], A[1][1], A[1][3], A[2][0], A[2][1], A[2][3]).det();
	V -= mat3d(A[0][0], A[0][1], A[0][2], A[1][0], A[1][1], A[1][2], A[2][0], A[2][1], A[2][2]).det();

	return V/6.0;
}

//-----------------------------------------------------------------------------
//! The residual is defined as the sum of the nodal residual and the
//! element residual
void FEUT4Domain::Residual(FESolidSolver *psolver, vector<double> &R)
{
	// Calculate the nodal contribution
	NodalResidual(psolver, R);

	// Calculate the element contribution
	ElementResidual(psolver, R);
}

//-----------------------------------------------------------------------------
//! This function calculates the nodal contribution to the residual
void FEUT4Domain::NodalResidual(FESolidSolver* psolver, vector<double>& R)
{

}

//-----------------------------------------------------------------------------
//! This function calculates the element contribution to the residual
void FEUT4Domain::ElementResidual(FESolidSolver* psolver, vector<double>& R)
{
	FEM& fem = psolver->m_fem;

	// element force vector
	vector<double> fe;

	int NE = m_Elem.size();
	for (int i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];

		// unpack the element
		UnpackElement(el);

		// get the element force vector and initialize it to zero
		int ndof = 3*el.Nodes();
		fe.assign(ndof, 0);

		// calculate internal force vector
		ElementInternalForces(el, fe);

		// apply body forces
		if (fem.UseBodyForces())
		{
			// Note that this function is in the base class
			// I'm not sure, but I don't think I have to overwrite this one
			BodyForces(fem, el, fe);
		}

		// assemble element 'fe'-vector into global R vector
		psolver->AssembleResidual(el.m_node, el.LM(), fe, R);
	}
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FEUT4Domain::ElementInternalForces(FESolidElement& el, vector<double>& fe)
{
	int i, n;

	// jacobian matrix, inverse jacobian matrix and determinants
	double Ji[3][3], detJt;

	double Gx, Gy, Gz;

	const double* Gr, *Gs, *Gt;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	double*	gw = el.GaussWeights();

	// repeat for all integration points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// calculate the jacobian
		el.invjact(Ji, n);
		detJt = el.detJt(n);

		detJt *= gw[n];

		// get the stress vector for this integration point
		mat3ds s = pt.s;

		// take the deviatoric component and multiply it
		// with the stabilization factor
//		s = s.dev()*m_alpha;
		
		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);

		for (i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

			// calculate internal force
			// the '-' sign is so that the internal forces get subtracted
			// from the global residual vector
			fe[3*i  ] -= ( Gx*s.xx() +
				           Gy*s.xy() +
					       Gz*s.xz() )*detJt;

			fe[3*i+1] -= ( Gy*s.yy() +
				           Gx*s.xy() +
					       Gz*s.yz() )*detJt;

			fe[3*i+2] -= ( Gz*s.zz() +
				           Gy*s.yz() +
					       Gx*s.xz() )*detJt;
		}
	}
}
