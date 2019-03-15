#pragma once
#include "FESurfacePairConstraint.h"
#include "FENLConstraint.h"
#include "FETimeStepController.h"
#include <vector>

//-----------------------------------------------------------------------------
class FEModel;
class FESolver;
class FEDomain;
class DumpStream;
class FEBoundaryCondition;

//-----------------------------------------------------------------------------
//! Base class for finite element analysis
class FECORE_API FEAnalysis : public FECoreBase
{
	FECORE_SUPER_CLASS

public:
	//! constructor
	FEAnalysis(FEModel* pfem);

	//! destructor
	virtual ~FEAnalysis();

	//! Initialization
	virtual bool Init() override;

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
	virtual void Serialize(DumpStream& ar) override;

	//! copy data from another analysis
	void CopyFrom(FEAnalysis* step);

public:
	void SetFESolver(FESolver* psolver);

	FESolver* GetFESolver();

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
	int GetPlotLevel();

	//! Set the output level
	void SetOutputLevel(int n);

	//! Get the output level
	int GetOutputLevel();

private:
	// Call the FE Solver
	// Returns an error code
	// 0 = all is well, continue
	// 1 = solver has failed, but try auto-time step
	// 2 = abort
	int CallFESolver();

	// initialize the solver
	bool InitSolver();

public:
	// --- Control Data ---
	//{
		int		m_nanalysis;	//!< analysis type
	//}

	// --- Time Step Data ---
	//{
		int		m_ntime;		//!< nr of timesteps
		double	m_final_time;	//!< end time for this time step
		double	m_dt;			//!< current time step 
		double	m_dt0;			//!< initial time step size
		double	m_tstart;		//!< start time
		double	m_tend;			//!< end time
		bool	m_bautostep;	//!< use auto stepper?

		FETimeStepController m_timeController;
	//}

	// --- Quasi-Newton Solver Variables ---
	//{
		int		m_ntotrhs;		//!< total nr of right hand side evaluations
		int		m_ntotref;		//!< total nr of stiffness reformations
		int		m_ntotiter;		//!< total nr of non-linear iterations
		int		m_ntimesteps;	//!< time steps completed
	//}

	// --- I/O Data ---
	//{
		int		m_nplot;		//!< plot level
		int		m_noutput;		//!< data output level
		int		m_nplot_stride;	//!< stride for plotting
		int		m_nplotRange[2];	//!< plot range
		bool	m_bplotZero;		//!< Force plotting of time step "zero"
	//}

private:
	// the FE solver
	FESolver*	m_psolver;	//!< pointer to solver class that will solve this step.
	bool		m_bactive;	//!< activation flag

protected:
	std::vector<int>				m_Dom;	//!< list of active domains for this analysis
	std::vector<FEModelComponent*>	m_MC;	//!< array of model components active during this step

	DECLARE_FECORE_CLASS();
};
