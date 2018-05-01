#include "FEHeatSolidDomain.h"
#include "FECore/FEModel.h"
#include <FECore/Integrate.h>

//-----------------------------------------------------------------------------
//! constructor
FEHeatSolidDomain::FEHeatSolidDomain(FEModel* pfem) : FESolidDomain(pfem), FEHeatDomain(pfem)
{
	m_pMat = 0;

	// list the degrees of freedom
	vector<int> dof;
	dof.push_back(pfem->GetDOFIndex("T"));
	SetDOFList(dof);
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
void FEHeatSolidDomain::HeatSource(FEGlobalVector& R, FEHeatSource& hs)
{
	vector<double> fe;
	vector<int> lm;
	int NE = Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = Element(i);
		int ne = el.Nodes();
		fe.resize(ne);
		ElementHeatSource(hs, el, fe);
		UnpackLM(el, lm);
		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
//! calculate element contribution to heat source term
void FEHeatSolidDomain::ElementHeatSource(FEHeatSource& hs, FESolidElement& el, vector<double>& fe)
{
	zero(fe);
	double* w = el.GaussWeights();
	int ne = el.Nodes();
	int ni = el.GaussPoints();
	for (int n = 0; n<ni; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		double* H = el.H(n);
		double J = detJt(el, n);
		for (int i = 0; i<ne; ++i)
		{
			double Q = hs.value(mp);
			fe[i] += Q*H[i] * J*w[n];
		}
	}
}

//-----------------------------------------------------------------------------
// Valuator class for the conductivity
class FEConductivity : public FEMaterialPointValue<mat3ds>
{
public:
	FEConductivity(FEHeatTransferMaterial* mat) : m_mat(mat){}

	mat3ds operator()(FEMaterialPoint& mp) { return m_mat->Conductivity(mp); }

private:
	FEHeatTransferMaterial* m_mat;
};

//-----------------------------------------------------------------------------
//! This function calculates the element stiffness matrix for a particular
//! element.
//!
void FEHeatSolidDomain::ElementConduction(FESolidElement& el, matrix& ke)
{
	// zero the matrix
	ke.zero();

	// set up the conductivity valuator
	FEConductivity D(m_pMat);

	// do the element integration
	IntegrateBDB(*this, el, D, ke);
}

//-----------------------------------------------------------------------------
void FEHeatSolidDomain::ElementCapacitance(FESolidElement &el, matrix &ke, double dt)
{
	// zero stiffness matrix
	ke.zero();

	// calculate data
	double c = m_pMat->Capacitance();
	double d = m_pMat->Density();
	double alpha = c*d/dt;

	// do the element integration
	IntegrateNCN(*this, el, alpha, ke);
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
