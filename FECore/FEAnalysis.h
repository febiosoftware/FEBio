#pragma once
#include "DumpFile.h"
#include "FECoreBase.h"
#include "FEBoundaryCondition.h"
#include "FESurfacePairInteraction.h"
#include "FENLConstraint.h"
#include <vector>

//-----------------------------------------------------------------------------
class FEModel;
class FESolver;
class FEDomain;

namespace FECore {

//-----------------------------------------------------------------------------
//! Base class for finite element analysis
class FEAnalysis : public FECoreBase
{
public:
	//! constructor
	FEAnalysis(FEModel* pfem, int ntype);

	//! destructor
	virtual ~FEAnalysis(){}

	//! Data initialization
	virtual bool Init();

	//! Reset analysis data
	virtual void Reset();

	//! Solve the analysis step
	virtual bool Solve();

	//! wrap it up
	virtual void Finish();

	//! Serialize data from and to a binary archive
	virtual void Serialize(DumpFile& ar);

	//! get the step type
	int GetType () { return m_ntype; }

	//! set the step type
	void SetType(int ntype) { m_ntype = ntype; }

public:
	FEModel& GetFEModel() { return m_fem; }

	FESolver* GetFESolver() { return m_psolver; }

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

	//! add a surface pair interaction to the analysis
	void AddSurfacePairInteraction(FESurfacePairInteraction* pci) { m_CI.push_back(pci); }

	//! Add a non-linear constraint to the analysis
	void AddConstraint(FENLConstraint* pnlc) { m_NLC.push_back(pnlc); }

public:
	//! sets the plot level
	void SetPlotLevel(int n) { m_nplot = n; }

	//! sets the plot stride
	void SetPlotStride(int n) { m_nplot_stride = n; }

	//! get the plot level
	int GetPlotLevel() { return m_nplot; }

	//! Sets the print level
	void SetPrintLevel(int n) { m_nprint = n; }

	//! get the print level
	int GetPrintLevel() { return m_nprint; }

protected:
	//! initialize constraint data
	bool InitConstraints();

	//! Do a running restart
	void Retry();

	//! Update Time step
	void AutoTimeStep(int niter);

	//! Adjust for must points
	double CheckMustPoints(double t, double dt);

public:
	// --- The FE Model
	//{
		FEModel&	m_fem;	//!< reference to FE model
	//}

	// --- Control Data ---
	//{
		int		m_nanalysis;	//!< analysis type
		int		m_istiffpr;		//!< calculate pressure stiffness \todo remove
		bool	m_baugment;		//!< use Lagrangian augmentation \todo move to solver class?
	//}

	// --- Time Step Data ---
	//{
		int		m_ntime;		//!< nr of timesteps
		double	m_final_time;	//!< end time for this time step
		double	m_dt;			//!< time step size
		double	m_dt0;			//!< initial time step size
		double	m_tstart;		//!< start time
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
		FESolver*	m_psolver;	//!< pointer to solver class that will solve this step.

		int		m_nretries;		//!< nr of retries tried so far
		int		m_maxretries;	//!< max nr of retries allowed per time step

		int		m_ntotrhs;		//!< total nr of right hand side evaluations
		int		m_ntotref;		//!< total nr of stiffness reformations
		int		m_ntotiter;		//!< total nr of non-linear iterations
		int		m_ntimesteps;	//!< time steps completed
	//}

	// --- I/O Data ---
	//{
		int		m_nprint;		//!< print level
		int		m_nplot;		//!< plot level
		int		m_nplot_stride;	//!< stride for plotting
		bool	m_bDump;		//!< create a restart file or not
	//}

protected:
	std::vector<int>						m_Dom;	//!< list of active domains for this analysis
	std::vector<FEBoundaryCondition*>		m_BC;	//!< array of boundary conditions
	std::vector<FESurfacePairInteraction* >	m_CI;	//!< active surface pair interactions
	std::vector<FENLConstraint*>			m_NLC;	//!< non-linear constraints

protected:
	int		m_nmust;		//!< current must-point
	int		m_ntype;		//!< step type
};

} // namespace FECore
