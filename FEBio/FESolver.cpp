// FESolver.cpp: implementation of the FESolver class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FESolver.h"
#include "fem.h"

///////////////////////////////////////////////////////////////////////////////
// FESolver Construction/Destruction
//

FESolver::FESolver(FEM& fem) : m_fem(fem), m_log(fem.m_log)
{
	// default values
	m_Rtol = 1e10;
	m_Dtol = 0.001;
	m_Etol = 0.01;
	m_LStol = 0.9;
	m_LSmin = 0.01;
	m_LSiter = 5;
	m_maxups = 10;
	m_maxref = 15;
	m_cmax   = 1e5;

	m_niter = 0;

	// Stiffness matrix and linear solver are allocated in Init()
	m_pK = 0;
	m_pM = 0;
	m_psolver = 0;
}

FESolver::~FESolver()
{
	delete m_pK;		// clean up stiffnes matrix data
	delete m_pM;		// clean up mass matrix data
	delete m_psolver;	// clean up linear solver data
}


///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FESolver::CreateStiffness
//  Creates the global stiffness matrix
//

bool FESolver::CreateStiffness(bool breset)
{
	// clean up the solver
	m_psolver->Destroy();

	// clean up the stiffness matrix
	m_pK->Clear();

	// create the stiffness matrix
	m_log.printf("===== reforming stiffness matrix:\n");
	if (m_pK->Create(m_fem, breset) == false) 
	{
		m_log.printf("FATAL ERROR: An error occured while building the stiffness matrix\n\n");
		return false;
	}
	else
	{
		// output some information about the direct linear solver
		int neq = m_fem.m_neq;
		int nnz = m_pK->NonZeroes();
		m_log.printf("\tNr of equations ........................... : %d\n", neq);
		m_log.printf("\tNr of nonzeroes in stiffness matrix ....... : %d\n", nnz);
		m_log.printf("\n");
	}

	// Do the preprocessing of the solver
	m_SolverTime.start();
	{
		if (!m_psolver->PreProcess(*m_pK)) throw FatalError();
	}
	m_SolverTime.stop();

	// done!
	return true;
}

//-----------------------------------------------------------------------------
//! Save data to dump file

void FESolver::Serialize(Archive& ar)
{
	if (ar.IsSaving())
	{
		ar << m_Dtol << m_Etol << m_Rtol << m_LSmin << m_LStol << m_LSiter;
		ar << m_maxups;
		ar << m_maxref;
		ar << m_cmax;

		ar << m_nrhs;
		ar << m_niter;
		ar << m_nref;
		ar << m_nups;
		ar << m_naug;
	}
	else
	{
		ar >> m_Dtol >> m_Etol >> m_Rtol >> m_LSmin >> m_LStol >> m_LSiter;
		ar >> m_maxups;
		ar >> m_maxref;
		ar >> m_cmax;

		ar >> m_nrhs;
		ar >> m_niter;
		ar >> m_nref;
		ar >> m_nups;
		ar >> m_naug;
	}
}
