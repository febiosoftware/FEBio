#pragma once
#include "FEBiphasicSolver.h"

//-----------------------------------------------------------------------------
// This class adds additional functionality to the FEBiphasicSolver to solve
// solute problems. 
class FEBIOMIX_API FEBiphasicSoluteSolver : public FEBiphasicSolver
{
public:
	//! con/descructor
	FEBiphasicSoluteSolver(FEModel* pfem);
	virtual ~FEBiphasicSoluteSolver(){}

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
	//! Calculates concentrated nodal forces (overridden from FESolidSolver2)
	//! (This function is called from FESolidSolver2::PrepStep)
	virtual void NodalForces(vector<double>& F, const FETimeInfo& tp) override;

	//! Calculates residual (overridden from FEBiphasicSolver)
	virtual bool Residual(vector<double>& R) override;

	//! calculates the global stiffness matrix (overridden from FESolidSolver2)
	virtual bool StiffnessMatrix() override;

	//! update kinematics
	virtual void UpdateKinematics(vector<double>& ui) override;

	//! Update solute data
	void UpdateSolute(vector<double>& ui);

protected:
	void GetConcentrationData(vector<double>& ci, vector<double>& ui, const int sol);

public:	// Parameters
	double	m_Ctol;			//!< concentration tolerance

public:
	// solute data
	vector< vector<double> >	m_ci;	//!< concentration increment vector
	vector< vector<double> >	m_Ci;	//!< Total concentration vector for iteration

	int	m_dofC;	//!< concentration dof
    int	m_dofD;	//!< shell concentration dof

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
