#pragma once
#include "DumpFile.h"
#include "FEBoundaryCondition.h"
#include "FEContactInterface.h"
#include <vector>

//-----------------------------------------------------------------------------
class FEModel;
class FENLSolver;
class FEDomain;

//-----------------------------------------------------------------------------
// This class is used as a method to communicate progress with the outside world
class Progress
{
public:
	Progress(){}
	virtual ~Progress() {}

	virtual void SetProgress(double f) = 0;
};

namespace FECore {

//-----------------------------------------------------------------------------
//! Base class for finite element analysis
class FEAnalysis
{
public:
	//! constructor
	FEAnalysis(FEModel& fem, int ntype) : m_fem(fem), m_ntype(ntype) {}

	//! destructor
	virtual ~FEAnalysis(){}

	//! Data initialization
	virtual bool Init() { return false; }

	//! Solve the analysis step
	virtual bool Solve(Progress& prg) { return false; }

	//! wrap it up
	virtual void Finish() {}

	//! Serialize data from and to a binary archive
	virtual void Serialize(DumpFile& ar) {}

	//! get the step type
	int GetType () { return m_ntype; }

	//! set the step type
	void SetType(int ntype) { m_ntype = ntype; }

public:
	//! Get active domains
	int Domains() { return m_Dom.size(); }

	//! Get active domain
	FEDomain* Domain(int i);

	//! Add a domain
	void AddDomain(int i) { m_Dom.push_back(i); }

	//! clear all domains
	void ClearDomains() { m_Dom.clear(); }

public:
	//! add a boundary condition to the analysis
	void AddBoundaryCondition(FEBoundaryCondition* pbc) { m_BC.push_back(pbc); }

	//! return number of boundary conditions
	int BoundaryConditions() { return (int) m_BC.size(); }

	//! add a boundary condition to the analysis
	void AddContactInterface(FEContactInterface* pci) { m_CI.push_back(pci); }

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
		int		m_nanalysis;	//!< analysis type
		int		m_istiffpr;		//!< calculate pressure stiffness (TODO remove)
		bool	m_baugment;		//!< use Lagrangian augmentation (TODO: move to solver class?)
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

protected:
	std::vector<int>					m_Dom;	//!< list of active domains for this analysis
	std::vector<FEBoundaryCondition*>	m_BC;	//!< array of boundary conditions
	std::vector<FEContactInterface* >	m_CI;	//!< active contact interfaces

protected:
	int		m_ntype;		// step type
};

} // namespace FECore
