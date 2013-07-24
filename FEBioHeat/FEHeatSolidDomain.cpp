#include "FEHeatSolidDomain.h"
#include "FEHeatTransferMaterial.h"
#include <FECore/FEMesh.h>

//-----------------------------------------------------------------------------
FEDomain* FEHeatSolidDomain::Clone()
{
	FEHeatSolidDomain* pd = new FEHeatSolidDomain(m_pMesh, m_pMat);
	pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
	return pd;
}

//-----------------------------------------------------------------------------
//! Unpack the element. That is, copy element data in traits structure

void FEHeatSolidDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.resize(N);

	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];
		FENode& node = m_pMesh->Node(n);

		int* id = node.m_ID;

		// get temperature equation number
		lm[i] = id[DOF_T];
	}
}

//-----------------------------------------------------------------------------
// Calculate the heat conduction matrix
void FEHeatSolidDomain::ConductionMatrix(FESolver* psolver)
{
	vector<int> lm;

	// loop over all elements in domain
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		int ne = el.Nodes();

		// build the element stiffness matrix
		matrix ke(ne, ne);
		ElementConduction(el, ke);

		// set up the LM matrix
		UnpackLM(el, lm);

		// assemble into global matrix
		psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}

//-----------------------------------------------------------------------------
// Calculate the capacitance matrix
void FEHeatSolidDomain::CapacitanceMatrix(FESolver* psolver, double dt)
{
	vector<int> lm;

	// loop over all elements in domain
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		int ne = el.Nodes();

		// element capacitance matrix
		matrix kc(ne, ne);
		ElementCapacitance(el, kc, dt);

		// set up the LM matrix
		UnpackLM(el, lm);

		// assemble into global matrix
		psolver->AssembleStiffness(el.m_node, lm, kc);
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the element stiffness matrix for a particular
//! element.
//!
void FEHeatSolidDomain::ElementConduction(FESolidElement& el, matrix& ke)
{
	int i, j, n;

	int ne = el.Nodes();
	int ni = el.GaussPoints();

	// global derivatives of shape functions
	// Gx = dH/dx
	const int EN = FEElement::MAX_NODES;
	double Gx[EN], Gy[EN], Gz[EN];

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

	FEHeatTransferMaterial& mat = dynamic_cast<FEHeatTransferMaterial&>(*m_pMat);

	// loop over all integration points
	for (n=0; n<ni; ++n)
	{
		// calculate jacobian
		detJt = invjact(el, Ji, n);

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
void FEHeatSolidDomain::ElementCapacitance(FESolidElement &el, matrix &ke, double dt)
{
	int i, j, n;

	int ne = el.Nodes();
	int ni = el.GaussPoints();

	// shape functions
	double* H;

	// jacobian
	double Ji[3][3], detJt;

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// zero stiffness matrix
	ke.zero();

	FEHeatTransferMaterial& mat = dynamic_cast<FEHeatTransferMaterial&>(*m_pMat);
	double c = mat.Capacitance();
	double d = mat.Density();
	double alpha = c*d/dt;

	// loop over all integration points
	for (n=0; n<ni; ++n)
	{
		// calculate jacobian
		detJt = invjact(el, Ji, n);

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
