#include "stdafx.h"
#include "FESolidSolver.h"
#include "fem.h"

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FESolidSolver::Init
// Allocates and initializes the data structures used by the FESolidSolver
//

bool FESolidSolver::Init()
{
	// initialize base class
	if (FESolver::Init() == false) return false;

	// get nr of equations
	int neq = m_fem.m_neq;

	// allocate vectors
	m_Fn.assign(neq, 0);
	m_Fd.assign(neq, 0);
	m_Fr.assign(neq, 0);
	m_Ui.assign(neq, 0);
	m_Ut.assign(neq, 0);

	// allocate poro-vectors
	if (m_fem.m_npeq > 0)
	{
		m_pi.assign(m_fem.m_npeq, 0);
		m_Pi.assign(m_fem.m_npeq, 0);
	}

	// we need to fill the total displacement vector m_Ut
	// TODO: I need to find an easier way to do this
	FEMesh& mesh = m_fem.m_mesh;
	int i, n;
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		// displacement dofs
		n = node.m_ID[0]; if (n >= 0) m_Ut[n] = node.m_rt.x - node.m_r0.x;
		n = node.m_ID[1]; if (n >= 0) m_Ut[n] = node.m_rt.y - node.m_r0.y;
		n = node.m_ID[2]; if (n >= 0) m_Ut[n] = node.m_rt.z - node.m_r0.z;

		// rotational dofs
		n = node.m_ID[3]; if (n >= 0) m_Ut[n] = node.m_Dt.x - node.m_D0.x;
		n = node.m_ID[4]; if (n >= 0) m_Ut[n] = node.m_Dt.y - node.m_D0.y;
		n = node.m_ID[5]; if (n >= 0) m_Ut[n] = node.m_Dt.z - node.m_D0.z;

		// pressure dofs
		n = node.m_ID[6]; if (n >= 0) m_Ut[n] = node.m_pt;
	}

	// initialize BFGS data
	m_bfgs.Init(neq, this, m_plinsolve);

	// set the create stiffness matrix flag
	m_breshape = true;

	return true;
}
