#pragma once
#include <FEBioMech/FESolidSolver.h>

//-----------------------------------------------------------------------------
// This class implements a solver for non-linear, thermo-elastic problems
class FEThermoElasticSolver : public FESolidSolver
{
public:
	// constructor
	FEThermoElasticSolver(FEModel* pfem);
	virtual ~FEThermoElasticSolver();

	// Initialization
	bool Init() override;

	// initialize equations
	bool InitEquations() override;

	// prepare first QN-iteration
	void PrepStep(const FETimeInfo& timeInfo) override;

	// Solver a time-step using a Quasi-Newton method
	bool Quasin(double time) override;

	// Residual
	bool Residual(vector<double>& R) override;

	//! incremental update
	void Update(vector<double>& ui) override;

protected:
	void GetDisplacementData(vector<double>& di, const vector<double>& ui);
	void GetTemperatureData(vector<double>& ti, const vector<double>& ui);

public:
	double	m_Ttol;	//!< temperature convergence tolerance

protected:
	int		m_ndeq;	//!< number of displacement equations
	int		m_nteq;	//!< number of temperature equations

	int		m_dofT;	//!< degree of freedom index for temperature

	// thermo data
	vector<double>	m_di;	//!< displacement increment vector
	vector<double>	m_Di;	//!< total displacement vector for iteration
	vector<double>	m_ti;	//!< temperature increment vector
	vector<double>	m_Ti;	//!< Total temperature vector for iteration

	DECLARE_PARAMETER_LIST();
};
