#include "FEHeatSolidDomain.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
//! constructor
FEHeatSolidDomain::FEHeatSolidDomain(FEModel* pfem) : FESolidDomain(&pfem->GetMesh()), FEHeatDomain(pfem)
{
	m_pMat = 0;

	// list the degrees of freedom
	vector<int> dof;
	dof.push_back(pfem->GetDOFIndex("T"));
	SetDOF(dof);
}

//-----------------------------------------------------------------------------
void FEHeatSolidDomain::SetMaterial(FEMaterial* pmat)
{
	m_pMat = dynamic_cast<FEHeatTransferMaterial*>(pmat);
	assert(m_pMat);
}

//-----------------------------------------------------------------------------
// Calculate the heat conduction matrix
void FEHeatSolidDomain::ConductionMatrix(FELinearSystem& ls)
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
		ls.AssembleLHS(lm, ke);
	}
}

//-----------------------------------------------------------------------------
// Calculate the capacitance matrix
void FEHeatSolidDomain::CapacitanceMatrix(FELinearSystem& ls, double dt)
{
	vector<int> lm;
	vector<double> fe;

	vector<int> dofs = GetDOFList();
	assert(dofs.size()==1);
	int dofT = dofs[0];

	// get the mesh
	FEMesh& mesh = *GetMesh();

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
		ls.AssembleLHS(lm, kc);

		// we also need to assemble this in the right-hand side
		fe.resize(ne, 0.0);
		for (int i = 0; i<ne; ++i)
		{
			fe[i] = mesh.Node(el.m_node[i]).get(dofT);
		}

		// multiply with me
		fe = kc*fe;

		// assemble this vector to the right-hand side
		ls.AssembleRHS(lm, fe);
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

	// loop over all integration points
	for (n=0; n<ni; ++n)
	{
		// calculate jacobian
		detJt = invjact(el, Ji, n);

		// evaluate the conductivity
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		mat3ds D = m_pMat->Conductivity(mp);

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

				DB[0] = D(0,0)*Gj[0] + D(0,1)*Gj[1] + D(0,2)*Gj[2];
				DB[1] = D(1,0)*Gj[0] + D(1,1)*Gj[1] + D(1,2)*Gj[2];
				DB[2] = D(2,0)*Gj[0] + D(2,1)*Gj[1] + D(2,2)*Gj[2];

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

	double c = m_pMat->Capacitance();
	double d = m_pMat->Density();
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

//-----------------------------------------------------------------------------
// This function updates the domain state data. 
// This function is called during FESolver::Update after the nodal variables are udpated.
void FEHeatSolidDomain::Update(const FETimeInfo& tp)
{
	FEMesh& mesh = *GetMesh();
	int NE = Elements();
	const int dof_T = m_dof[DOF_T];
	for (int j=0; j<NE; ++j)
	{
		FESolidElement& el = Element(j);
		int ni = el.GaussPoints();
		int ne = el.Nodes();

		// get the nodal temperatures
		double T[FEElement::MAX_NODES];
		for (int n=0; n<ne; ++n) T[n] = mesh.Node(el.m_node[n]).get(dof_T);

		// calculate heat flux for each integration point
		for (int n=0; n<ni; ++n)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);
			FEHeatMaterialPoint* pt = (mp.ExtractData<FEHeatMaterialPoint>());
			assert(pt);

			vec3d gradT = gradient(el, T, n);
			mat3ds D = m_pMat->Conductivity(mp);
			pt->m_q = -(D*gradT);
		}
	}
}
