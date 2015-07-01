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
	bool Init();

	// initialize equations
	bool InitEquations();

	// prepare first QN-iteration
	void PrepStep(double time);

	// Solver a time-step using a Quasi-Newton method
	bool Quasin(double time);

	// Residual
	bool Residual(vector<double>& R);

	//! incremental update
	void Update(vector<double>& ui);

protected:
	void GetDisplacementData(vector<double>& di, const vector<double>& ui);
	void GetTemperatureData(vector<double>& ti, const vector<double>& ui);

public:
	double	m_Ttol;	//!< temperature convergence tolerance

protected:
	int		m_ndeq;	//!< number of displacement equations
	int		m_nteq;	//!< number of temperature equations

	// thermo data
	vector<double>	m_di;	//!< displacement increment vector
	vector<double>	m_Di;	//!< total displacement vector for iteration
	vector<double>	m_ti;	//!< temperature increment vector
	vector<double>	m_Ti;	//!< Total temperature vector for iteration

	DECLARE_PARAMETER_LIST();
};
