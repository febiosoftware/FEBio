#include "stdafx.h"
#include "FESolidSolver.h"
#include "fem.h"

///////////////////////////////////////////////////////////////////////////////
// FESolidSolver Construction/Destruction
//

FESolidSolver::FESolidSolver(FEM& fem) : FESolver(fem)
{
	// default values
	m_Rtol = 0;	// deactivate residual convergence 
	m_Dtol = 0.001;
	m_Etol = 0.01;
	m_Ptol = 0.01;

	m_niter = 0;
}

//-----------------------------------------------------------------------------
//! Save data to dump file

void FESolidSolver::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_Dtol << m_Etol << m_Rtol << m_Ptol;
		ar << m_bfgs.m_LStol << m_bfgs.m_LSiter << m_bfgs.m_LSmin;
		ar << m_bfgs.m_maxups;
		ar << m_bfgs.m_maxref;
		ar << m_bfgs.m_cmax;

		ar << m_nrhs;
		ar << m_niter;
		ar << m_nref << m_ntotref;
		ar << m_bfgs.m_nups;
		ar << m_naug;
	}
	else
	{
		ar >> m_Dtol >> m_Etol >> m_Rtol >> m_Ptol;
		ar >> m_bfgs.m_LStol >> m_bfgs.m_LSiter >> m_bfgs.m_LSmin;
		ar >> m_bfgs.m_maxups;
		ar >> m_bfgs.m_maxref;
		ar >> m_bfgs.m_cmax;

		ar >> m_nrhs;
		ar >> m_niter;
		ar >> m_nref >> m_ntotref;
		ar >> m_bfgs.m_nups;
		ar >> m_naug;
	}
}

