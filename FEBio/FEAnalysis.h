#pragma once

#include "FECore/DumpFile.h"
#include "FECore/FEBoundaryCondition.h"
#include <vector>
using namespace std;

// forward declaration of FEM
class FEM;
class FESolver;

//-----------------------------------------------------------------------------
//! The FEAnalysis contains the data that describes a complete, single analysis.

class FEAnalysis
{
public:
	//! constructor
	FEAnalysis(FEM& fem);

	//! descructor
	virtual ~FEAnalysis(void);

	//! sets the plot level
	void SetPlotLevel(int n) { m_nplot = n; }

	//! get the plot level
	int GetPlotLevel() { return m_nplot; }

	//! Sets the print level
	void SetPrintLevel(int n) { m_nprint = n; }

	//! get the print level
	int GetPrintLevel() { return m_nprint; }

	//! initialize step
	bool Init();

	//! Solve the analysis step
	bool Solve();

	//! wrap it up
	void Finish();

	//! Serialize data from and to a binary archive
	void Serialize(DumpFile& ar);

	//! add a boundary condition to the analysis
	void AddBoundaryCondition(FEBoundaryCondition* pbc) { m_BC.push_back(pbc); }

protected:
	//! Do a running restart
	void Retry();

	//! Update Time step
	void AutoTimeStep(int niter);

public:
	FEM&	m_fem;	//!< reference to parent fem object

	// --- Control Data ---
	//{
		int		m_nModule;		//!< module type
		int		m_nanalysis;	//!< analysis type
		int		m_istiffpr;		//!< calculate pressure stiffness
		bool	m_baugment;		//!< use Lagrangian augmentation
		double	m_hg;			//!< hourglass parameter
	//}

	// --- Time Step Data ---
	//{
		int		m_ntime;		//!< nr of timesteps
		double	m_final_time;	//!< end time for this time step
		double	m_dt;			//!< time step size
		double	m_dt0;			//!< initial time step size
		double	m_tend;			//!< end time
		bool	m_bautostep;	//!< use auto stepper?
		int		m_iteopt;		//!< optimum nr of iterations
		double	m_dtmin;		//!< min time step size
		double	m_dtmax;		//!< max time step size
		int		m_nmplc;		//!< must point load curve number
		int		m_naggr;		//!< aggressivness parameter
	//}

	// --- Boundary conditions data ---
	//{
		vector<FEBoundaryCondition*>	m_BC;	//!< array of boundary conditions
	//}

	// --- Quasi-Newton Solver Variables ---
	//{
		// the FE solver
		FESolver*	m_psolver;

		// Newton - Raphson iteration data
		int		m_nretries;		//!< nr of retries tried so far
		int		m_maxretries;	//!< max nr of retries allowed per time step

		int		m_ntotrhs;		//!< total nr of right hand side evaluations
		int		m_ntotref;		//!< total nr of stiffness reformations
		int		m_ntotiter;		//!< total nr of non-linear iterations
		int		m_ntimesteps;	//!< time steps completed
	//}

	// --- I/O Data ---
	//{
		int		m_nprint;	//!< print level
		int		m_nplot;	//!< plot level
		bool	m_bDump;	//!< create a restart file or not
	//}
};
