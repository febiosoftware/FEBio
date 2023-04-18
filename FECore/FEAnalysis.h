/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include "FECoreBase.h"
#include "FECoreClass.h"
#include <vector>

//-----------------------------------------------------------------------------
class FEModel;
class FESolver;
class FEDomain;
class DumpStream;
class FEStepComponent;
class FETimeStepController;

//-----------------------------------------------------------------------------
//! Base class for finite element analysis
class FECORE_API FEAnalysis : public FECoreBase
{
	FECORE_SUPER_CLASS(FEANALYSIS_ID)
	FECORE_BASE_CLASS(FEAnalysis);

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
	//! add a step component
	void AddStepComponent(FEStepComponent* pmc);

	//! return number of model components
	int StepComponents() const;

	//! get a step component
	FEStepComponent* GetStepComponent(int i);

public:
	//! sets the plot level
	void SetPlotLevel(int n);

	//! sets the plot stride
	void SetPlotStride(int n);

	//! sets the plot range
	void SetPlotRange(int n0, int n1);

	//! sets the zero-state plot flag
	void SetPlotZeroState(bool b);

	//! sets the plot hint
	void SetPlotHint(int plotHint);

	//! get the plot hint
	int GetPlotHint() const;

	//! get the plot level
	int GetPlotLevel();

	//! Set the output level
	void SetOutputLevel(int n);

	//! Get the output level
	int GetOutputLevel();

	// initialize the solver
	bool InitSolver();

	// Call the FE Solver to solve the time step
	// Returns an error code
	// 0 = all is well, continue
	// 1 = solver has failed, but try auto-time step
	// 2 = abort
	int SolveTimeStep();

public:
	// --- Control Data ---
	//{
		int		m_nanalysis;		//!< analysis type
		bool	m_badaptorReSolve;	//!< resolve analysis after mesh adaptor phase
	//}

	// --- Time Step Data ---
	//{
		int		m_ntime;		//!< nr of timesteps
		double	m_final_time;	//!< end time for this time step
		double	m_dt;			//!< current time step 
		double	m_dt0;			//!< initial time step size
		double	m_tstart;		//!< start time
		double	m_tend;			//!< end time

		FETimeStepController* m_timeController;
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
		int		m_plotHint;			//!< the plot mode
	//}

private:
	// the FE solver
	FESolver*	m_psolver;	//!< pointer to solver class that will solve this step.
	bool		m_bactive;	//!< activation flag

protected:
	std::vector<int>				m_Dom;	//!< list of active domains for this analysis
	std::vector<FEStepComponent*>	m_MC;	//!< array of model components active during this step

	DECLARE_FECORE_CLASS();
};
