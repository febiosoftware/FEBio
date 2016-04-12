#include "FECoupledHeatSolidSolver.h"
#include "FEBioMech/FELinearSolidDomain.h"
#include "FEBioMech/FELinearElastic.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>

//-----------------------------------------------------------------------------
//! class constructor
FECoupledHeatSolidSolver::FECoupledHeatSolidSolver(FEModel* pfem) : FESolver(pfem), m_Heat(pfem), m_Solid(pfem)
{
}

//-----------------------------------------------------------------------------
void FECoupledHeatSolidSolver::Clean()
{
	m_Heat.Clean();
	m_Solid.Clean();
}

//-----------------------------------------------------------------------------
//! Initialization
bool FECoupledHeatSolidSolver::Init()
{
	if (FESolver::Init() == false) return false;

	// Initialize heat solver
	if (m_Heat.Init() == false) return false;

	// Initialize solid solver
	if (m_Solid.Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
//! Initialize equation system.
//! This solver doesn't manage a linear system, but the two child problems
//! do so we just call the corresponding function.
bool FECoupledHeatSolidSolver::InitEquations()
{
	// Initialize equations for heat problem
	if (m_Heat.InitEquations() == false) return false;

	// Initialize equations for solid problem
	if (m_Solid.InitEquations() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
bool FECoupledHeatSolidSolver::SolveStep(double time)
{
	// First we solve the heat problem
	if (m_Heat.SolveStep(time) == false) return false;

	// Now we project the nodal temperatures
	// to the integration points, and use them to set up an
	// initial stress for the linear solid solver
	CalculateInitialStresses();

	// Next, we solve the linear solid problem
	if (m_Solid.SolveStep(time) == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
//! \todo This won't work anymore since the domains are stored on the FEAnalysis
void FECoupledHeatSolidSolver::CalculateInitialStresses()
{
	const int dof_T = m_fem.GetDOFS().GetDOF("T");
	if (dof_T == -1) { assert(false); return; }

	FEAnalysis* pstep = m_fem.GetCurrentStep();

	FEHeatSolidDomain&   dh = dynamic_cast<FEHeatSolidDomain  &>(*pstep->Domain(0));
	FELinearSolidDomain& ds = dynamic_cast<FELinearSolidDomain&>(*pstep->Domain(0));

	FELinearElastic* pmat = dynamic_cast<FELinearElastic*>(ds.GetMaterial());
	assert(pmat);

	const double alpha = 1.0;

	double tn[9], tj;

	FEMesh& mesh = *ds.GetMesh();

	for (int i=0; i<ds.Elements(); ++i)
	{
		FESolidElement& el = ds.Element(i);
		int nint = el.GaussPoints();
		int neln = el.Nodes();

		// get the nodal temperatures
		for (int j=0; j<neln; ++j) tn[j] = mesh.Node(el.m_node[j]).get(dof_T);

		// loop over integration points
		for (int j=0; j<nint; ++j)
		{
			// get the material point
			// TODO: do I need to do anything to the material point?
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);

			// get the material tangent at this point
			tens4ds C = pmat->Tangent(mp);

			// evaluate the temperature at this point
			tj = el.Evaluate(tn, j);

			// set up the thermal strain
			double g = -alpha*tj;
			mat3ds e = mat3ds(g,g,g,0,0,0);

			// set the initial stress
			FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
			pt.m_s0 = C.dot(e);
		}
	}
}

//-----------------------------------------------------------------------------
void FECoupledHeatSolidSolver::Update(std::vector<double> &u)
{
	// Nothing to do here.
}

//-----------------------------------------------------------------------------
void FECoupledHeatSolidSolver::Serialize(DumpStream &ar)
{
	FESolver::Serialize(ar);
	m_Heat.Serialize(ar);
	m_Solid.Serialize(ar);
}
