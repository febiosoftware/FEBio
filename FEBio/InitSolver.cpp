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
	m_Fn.create(neq); m_Fn.zero();
	m_Fd.create(neq); m_Fd.zero();
	m_Fr.create(neq); m_Fr.zero();
	m_ui.create(neq); m_ui.zero();
	m_Ui.create(neq); m_Ui.zero();
	m_Ut.create(neq); m_Ut.zero();
	m_R0.create(neq); m_R0.zero();
	m_R1.create(neq); m_R1.zero();

	// allocate poro-vectors
	if (m_fem.m_npeq > 0)
	{
		m_pi.create(m_fem.m_npeq); m_pi.zero();
		m_Pi.create(m_fem.m_npeq); m_Ui.zero();
	}

	// if we have traction constraints
	// we store the internal forces seperately
	m_Ti.create(neq); m_Ti.zero();

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

	// allocate storage for BFGS update vectors
	m_V.Create(m_maxups, neq);
	m_W.Create(m_maxups, neq);

	m_D .create(neq);
	m_G .create(neq);
	m_H .create(neq);

	// set the create stiffness matrix flag
	m_breshape = true;

	return true;
}
