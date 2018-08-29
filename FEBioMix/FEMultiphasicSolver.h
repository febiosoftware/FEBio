#pragma once
#include "FEBiphasicSolver.h"

//-----------------------------------------------------------------------------
// This class adds additional functionality to the FESolidSolver2 to solve
// solute problems. 
class FEMultiphasicSolver : public FESolidSolver2
{
public:
	//! con/descructor
	FEMultiphasicSolver(FEModel* pfem);
	virtual ~FEMultiphasicSolver(){}

	//! Initialize data structures
	bool Init() override;

	//! Initialize equations
	bool InitEquations() override;

	//! prepares the data for the first QN iteration
	void PrepStep() override;

	//! Performs a Newton-Raphson iteration
	bool Quasin() override;

	//! serialize data to/from dump file
	void Serialize(DumpStream& ar) override;

public:
	//! update contact
	virtual void UpdateContact() override;

	//! update kinematics
	virtual void UpdateKinematics(vector<double>& ui) override;

	//! Update poroelastic data
	void UpdatePoro(vector<double>& ui);

	//! Update solute data
	void UpdateSolute(vector<double>& ui);

public:
	//! Calculates concentrated nodal forces (overridden from FESolidSolver2)
	//! (This function is called from FESolidSolver2::PrepStep)
	virtual void NodalForces(vector<double>& F, const FETimeInfo& tp) override;

	//! Calculates residual (overridden from FESolidSolver2)
	virtual bool Residual(vector<double>& R) override;

	//! calculates the global stiffness matrix (overridden from FESolidSolver2)
	virtual bool StiffnessMatrix() override;

protected:
	void GetDisplacementData(vector<double>& di, vector<double>& ui);
	void GetPressureData(vector<double>& pi, vector<double>& ui);
	void GetConcentrationData(vector<double>& ci, vector<double>& ui, const int sol);

public:	// Parameters
	double	m_Ctol;				//!< concentration tolerance
	double	m_Ptol;				//!< pressure tolerance
	bool	m_forcePositive;	//!< force conentrations to remain positive

public:
	// equation numbers
	int		m_ndeq;				//!< number of equations related to displacement dofs
	int		m_npeq;				//!< number of equations related to pressure dofs
	vector<int>		m_nceq;	//!< number of equations related to concentration dofs

	// poro data
	vector<double>	m_di;	//!< displacement increment vector
	vector<double>	m_Di;	//!< total displacement vector for iteration
	vector<double>	m_pi;	//!< pressure increment vector
	vector<double>	m_Pi;	//!< Total pressure vector for iteration

	// solute data
	vector< vector<double> >	m_ci;	//!< concentration increment vector
	vector< vector<double> >	m_Ci;	//!< Total concentration vector for iteration

protected:
	int	m_dofP;	//!< pressure dof index
    int	m_dofQ;	//!< shell pressure dof index
    int	m_dofC;	//!< concentration dof index
    int	m_dofD;	//!< shell concentration dof
    
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
