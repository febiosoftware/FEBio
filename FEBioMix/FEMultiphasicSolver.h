#pragma once
#include "FEBiphasicSolver.h"

//-----------------------------------------------------------------------------
// This class adds additional functionality to the FESolidSolver to solve
// solute problems. 
class FEMultiphasicSolver : public FESolidSolver
{
public:
	//! con/descructor
	FEMultiphasicSolver(FEModel* pfem);
	virtual ~FEMultiphasicSolver(){}

	//! Initialize data structures
	bool Init();

	//! Initialize equations
	bool InitEquations();

	//! prepares the data for the first QN iteration
	virtual void PrepStep(double time);

	//! Performs a Newton-Raphson iteration
	bool Quasin(double time);

	//! serialize data to/from dump file
	void Serialize(DumpFile& ar);

public:
	//! update contact
	virtual void UpdateContact();

	//! update kinematics
	virtual void UpdateKinematics(vector<double>& ui);

	//! Update poroelastic data
	void UpdatePoro(vector<double>& ui);

	//! Update solute data
	void UpdateSolute(vector<double>& ui);

public:
	//! Calculates residual (overridden from FESolidSolver)
	virtual bool Residual(vector<double>& R);

	//! calculates the global stiffness matrix (overridden from FESolidSolver)
	virtual bool StiffnessMatrix(const FETimePoint& tp);

protected:
	void GetDisplacementData(vector<double>& di, vector<double>& ui);
	void GetPressureData(vector<double>& pi, vector<double>& ui);
	void GetConcentrationData(vector<double>& ci, vector<double>& ui, const int sol);

public:	// Parameters
	double	m_Ctol;			//!< concentration tolerance
	double	m_Ptol;			//!< pressure tolerance

public:
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
	vector< vector<double> >	m_ci;	//!< concentration increment vector
	vector< vector<double> >	m_Ci;	//!< Total concentration vector for iteration

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
