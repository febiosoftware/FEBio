#pragma once
#include "FESolidSolver.h"

//-----------------------------------------------------------------------------
// This class adds additional functionality to the FESolidSolver to solve
// biphasic problems. 
class FEBiphasicSolver : public FESolidSolver
{
public:
	//! constructor
	FEBiphasicSolver(FEModel& fem);
	virtual ~FEBiphasicSolver() {}

	//! Initialize data structures
	bool Init();

	//! Initialize linear equation system
	bool InitEquations();

	//! prepares the data for the first QN iteration
	virtual void PrepStep(double time);

	//! Performs a Newton-Raphson iteration
	bool Quasin(double time);

	//! serialize data to/from dump file
	void Serialize(DumpFile& ar);

	//! update contact
	void UpdateContact();

public:
	//! Calculates residual (overridden from FESolidSolver)
	virtual bool Residual(vector<double>& R);

protected:
	void GetDisplacementData(vector<double>& di, vector<double>& ui);
	void GetPressureData(vector<double>& pi, vector<double>& ui);

public:
	// additional convergence norms
	double	m_Ptol;			//!< pressure tolerance

	// equation numbers
	int		m_ndeq;				//!< number of equations related to displacement dofs
	int		m_npeq;				//!< number of equations related to pressure dofs
	int		m_nceq[MAX_CDOFS];	//!< number of equations related to concentration dofs

	// poro data
	vector<double>	m_di;	//!< displacement increment vector
	vector<double>	m_Di;	//!< total displacement vector for iteration
	vector<double>	m_pi;	//!< pressure increment vector
	vector<double>	m_Pi;	//!< Total pressure vector for iteration

	// solute data
	vector<double>	m_ci[MAX_CDOFS];	//!< concentration increment vector
	vector<double>	m_Ci[MAX_CDOFS];	//!< Total concentration vector for iteration
};
