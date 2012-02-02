#pragma once
#include "DumpFile.h"

class FEModel;
class FENLSolver;

//-----------------------------------------------------------------------------
//! Base class for finite element analysis
class FEAnalysis
{
public:
	FEAnalysis(FEModel& fem) : m_fem(fem) {}

	virtual ~FEAnalysis(){}

	//! Data initialization
	virtual bool Init() { return false; }

	//! Solve the analysis step
	virtual bool Solve() { return false; }

	//! wrap it up
	virtual void Finish() {}

	//! Serialize data from and to a binary archive
	virtual void Serialize(DumpFile& ar) {}

public:
	//! sets the plot level
	void SetPlotLevel(int n) { m_nplot = n; }

	//! get the plot level
	int GetPlotLevel() { return m_nplot; }

	//! Sets the print level
	void SetPrintLevel(int n) { m_nprint = n; }

	//! get the print level
	int GetPrintLevel() { return m_nprint; }

public:
	// --- The FE Model
	//{
		FEModel&	m_fem;
	//}

	// --- Control Data ---
	//{
		int		m_nModule;		//!< module type
		int		m_nanalysis;	//!< analysis type
		int		m_istiffpr;		//!< calculate pressure stiffness (TODO remove)
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
		double	m_ddt;			//!< used by auto-time stepper
		int		m_nmplc;		//!< must point load curve number
		int		m_naggr;		//!< aggressivness parameter
	//}

	// --- Quasi-Newton Solver Variables ---
	//{
		// the FE solver
		FENLSolver*	m_psolver;

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
