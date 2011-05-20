#include "stdafx.h"
#include "FEHeatSolidDomain.h"
#include "FEHeatSolver.h"
#include "FEIsotropicFourier.h"

//-----------------------------------------------------------------------------
//! Unpack the element. That is, copy element data in traits structure

void FEHeatSolidDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.resize(N*MAX_NDOFS);

	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];
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
		
		// concentration dofs
		lm[11*N+i] = id[11];
	}
}

//-----------------------------------------------------------------------------
//! Unpack the element. That is, copy element data in traits structure

void FEHeatSolidDomain::UnpackElement(FEElement& el, unsigned int nflag)
{
	int i, n;

	vec3d* rt = el.rt();
	vec3d* r0 = el.r0();
	vec3d* vt = el.vt();
	double* pt = el.pt();
	double* ct = el.ct();

	int N = el.Nodes();

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
		
		// current nodal concentrations
		ct[i] = node.m_ct;
	}

	// unpack the traits data
	el.UnpackTraitsData(nflag);
}

//-----------------------------------------------------------------------------
void FEHeatSolidDomain::HeatStiffnessMatrix(FEHeatSolver* psolver)
{
	int i, j, k;
	vector<int> lm(8);
	vector<int> elm;

	FEM& fem = psolver->m_fem;

	for (i=0; i<(int) m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		UnpackElement(el);

		int ne = el.Nodes();

		// build the element stiffness matrix
		matrix ke(ne, ne);
		ConductionStiffness(fem, el, ke);

		// set up the LM matrix
		UnpackLM(el, elm);
		for (j=0; j<ne; ++j) lm[j] = elm[10*ne + j];

		if (fem.m_pStep->m_nanalysis == FE_DYNAMIC) 
		{
			matrix kc(ne, ne);
			CapacitanceStiffness(fem, el, kc);

			// add capacitance matrix to conduction stiffness
			ke += kc;

			// subtract from RHS
			for (j=0; j<ne; ++j)
			{
				if (lm[j] >= 0)
				{
					double q = 0;
					for (k=0; k<ne; ++k)
					{
						if (lm[k] >= 0) q += kc[j][k]*psolver->m_Tp[lm[k]];
						else if (-lm[k]-2 >= 0) q += kc[j][k]*psolver->m_Tp[-lm[k]-2];
					}

					psolver->m_R[lm[j]] += q;
				}
			}
		}

		// assemble into global matrix
		psolver->AssembleStiffness(ke, lm);
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the element stiffness matrix for a particular
//! element.
//!
void FEHeatSolidDomain::ConductionStiffness(FEM& fem, FESolidElement& el, matrix& ke)
{
	int i, j, n;

	int ne = el.Nodes();
	int ni = el.GaussPoints();

	// global derivatives of shape functions
	// NOTE: hard-coding of hex elements!
	// Gx = dH/dx
	double Gx[8], Gy[8], Gz[8];
	double Gr, Gs, Gt;
	double Gi[3], Gj[3];
	double DB[3];

	// jacobian
	double Ji[3][3], detJt;

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// zero stiffness matrix
	ke.zero();

	// conductivity matrix
	double D[3][3];

	FEIsotropicFourier& mat = dynamic_cast<FEIsotropicFourier&>(*fem.GetMaterial(el.GetMatID()));

	// loop over all integration points
	for (n=0; n<ni; ++n)
	{
		// calculate jacobian
		el.invjact(Ji, n);
		detJt = el.detJt(n);

		// evaluate the conductivity
		mat.Conductivity(D);

		for (i=0; i<ne; ++i)
		{
			Gr = el.Gr(n)[i];
			Gs = el.Gs(n)[i];
			Gt = el.Gt(n)[i];

			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx[i] = Ji[0][0]*Gr+Ji[1][0]*Gs+Ji[2][0]*Gt;
			Gy[i] = Ji[0][1]*Gr+Ji[1][1]*Gs+Ji[2][1]*Gt;
			Gz[i] = Ji[0][2]*Gr+Ji[1][2]*Gs+Ji[2][2]*Gt;
		}		

		for (i=0; i<ne; ++i)
		{
			Gi[0] = Gx[i];
			Gi[1] = Gy[i];
			Gi[2] = Gz[i];

			for (j=0; j<ne; ++j)
			{
				Gj[0] = Gx[j];
				Gj[1] = Gy[j];
				Gj[2] = Gz[j];

				DB[0] = D[0][0]*Gj[0] + D[0][1]*Gj[1] + D[0][2]*Gj[2];
				DB[1] = D[1][0]*Gj[0] + D[1][1]*Gj[1] + D[1][2]*Gj[2];
				DB[2] = D[2][0]*Gj[0] + D[2][1]*Gj[1] + D[2][2]*Gj[2];

				ke[i][j] += (Gi[0]*DB[0] + Gi[1]*DB[1] + Gi[2]*DB[2] )*detJt*gw[n];
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEHeatSolidDomain::CapacitanceStiffness(FEM& fem, FESolidElement &el, matrix &ke)
{
	int i, j, n;

	int ne = el.Nodes();
	int ni = el.GaussPoints();

	// shape functions
	// NOTE: hard-coding of hex elements!
	double* H;

	// jacobian
	double Ji[3][3], detJt;

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// zero stiffness matrix
	ke.zero();

	double dt = fem.m_pStep->m_dt;

	FEIsotropicFourier& mat = dynamic_cast<FEIsotropicFourier&>(*fem.GetMaterial(el.GetMatID()));
	double alpha = mat.m_c*mat.m_rho / dt;

	// loop over all integration points
	for (n=0; n<ni; ++n)
	{
		// calculate jacobian
		el.invjact(Ji, n);
		detJt = el.detJt(n);
		H = el.H(n);

		for (i=0; i<ne; ++i)
		{
			for (j=0; j<ne; ++j)
			{
				ke[i][j] += H[i]*H[j]*alpha*detJt*gw[n];
			}
		}
	}
}
