#pragma once
#include "FESurfacePairInteraction.h"
#include "FENLConstraint.h"
#include <vector>

//-----------------------------------------------------------------------------
class FEModel;
class FESolver;
class FEDomain;
class DumpStream;
class FEBoundaryCondition;

//-----------------------------------------------------------------------------
//! Base class for finite element analysis
class FEAnalysis : public FECoreBase
{
public:
	//! constructor
	FEAnalysis(FEModel* pfem);

	//! destructor
	virtual ~FEAnalysis();

	//! Initialization
	virtual bool Init();

	//! Activation
	virtual bool Activate();

	//! See if this step is active
	bool IsActive();

	//! Reset analysis data
	virtual void Reset();

	//! Solve the analysis step
	virtual bool Solve();

	//! wrap it up
	virtual void Deactivate();

	//! Serialize data from and to a binary archive
	virtual void Serialize(DumpStream& ar);

public:
	FEModel& GetFEModel() { return m_fem; }

	void SetFESolver(FESolver* psolver);

	FESolver* GetFESolver() { return m_psolver; }

public:
	//! Get active domains
	int Domains() { return (int)m_Dom.size(); }

	//! Get active domain
	FEDomain* Domain(int i);

	//! Add a domain
	void AddDomain(int i) { m_Dom.push_back(i); }

	//! clear all domains
	void ClearDomains() { m_Dom.clear(); }

public:
	//! add a model component
	void AddModelComponent(FEModelComponent* pmc);

	//! return number of model components
	int ModelComponents() const;

public:
	//! sets the plot level
	void SetPlotLevel(int n);

	//! sets the plot stride
	void SetPlotStride(int n);

	//! sets the plot range
	void SetPlotRange(int n0, int n1);

	//! sets the zero-state plot flag
	void SetPlotZeroState(bool b);

	//! get the plot level
	int GetPlotLevel() { return m_nplot; }

	//! Sets the print level
	void SetPrintLevel(int n) { m_nprint = n; }

	//! get the print level
	int GetPrintLevel() { return m_nprint; }

	//! Set the output level
	void SetOutputLevel(int n) { m_noutput = n; }

	//! Get the output level
	int GetOutputLevel() { return m_noutput; }

	//! Set the dump level (for cold restarts)
	void SetDumpLevel(int n) { m_ndump = n; }

	//! get the dump level
	int GetDumpLevel() { return m_ndump; }

protected:
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
	//}

	// --- Time Step Data ---
	//{
		int		m_ntime;		//!< nr of timesteps
		double	m_final_time;	//!< end time for this time step
		double	m_dt;			//!< time step size
		double	m_dt0;			//!< initial time step size
		double	m_dtp;			//!< previous time step size
		double	m_tstart;		//!< start time
		double	m_tend;			//!< end time
		bool	m_bautostep;	//!< use auto stepper?
		int		m_iteopt;		//!< optimum nr of iterations
		double	m_dtmin;		//!< min time step size
		double	m_dtmax;		//!< max time step size
		double	m_ddt;			//!< used by auto-time stepper
		int		m_nmplc;		//!< must point load curve number
		int		m_naggr;		//!< aggressivness parameter
		int		m_nmust;		//!< current must-point
		int		m_next_must;	//!< next must-point to visit
	//}

	// --- Quasi-Newton Solver Variables ---
	//{
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
		int		m_noutput;		//!< data output level
		int		m_nplot_stride;	//!< stride for plotting
		int		m_nplotRange[2];	//!< plot range
		bool	m_bplotZero;		//!< Force plotting of time step "zero"
		int		m_ndump;		//!< create a restart file or not
	//}

private:
	// the FE solver
	FESolver*	m_psolver;	//!< pointer to solver class that will solve this step.
	bool		m_bactive;	//!< activation flag

protected:
	std::vector<int>				m_Dom;	//!< list of active domains for this analysis
	std::vector<FEModelComponent*>	m_MC;	//!< array of model components active during this step

	DECLARE_PARAMETER_LIST();
};
