#pragma once
#include "FECore/FESolver.h"
#include "FECore/FEGlobalVector.h"
#include "FECore/FETypes.h"

//-----------------------------------------------------------------------------
//! This class implements a solver for solid mechanics problems that uses
//! the conjugate gradient method to solve the nonlinear finite element equations
class FECGSolidSolver : public FESolver
{
public:
	//! constructor
	FECGSolidSolver(FEModel* pfem);

	//! initialization
	bool Init();

	//! clean up
	void Clean();

	//! Performs a CG step
	bool SolveStep(double time);

	//! update nodal positions, velocities, accelerations, etc.
	void UpdateKinematics(vector<double>& ui);

	//! assemble global stiffness matrix
	//! \todo Get rid of this function
	virtual void AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke) {}

	// Initialize linear equation system (TODO: Is this the right place to do this?)
	// \todo Can I make this part of the Init function?
	virtual bool InitEquations();

protected:
	//! Update the stresses
	void UpdateStresses();

	//! Update contact
	void UpdateContact();

	//! update constraints
	void UpdateConstraints();

	//! update rigid bodies
	void UpdateRigidBodies(vector<double>& ui);

	//! Evaluate the residual
	bool Residual(vector<double>& R);

	//! assemble into the residual
	void AssembleResidual(int node_id, int dof, double f, vector<double>& R);

	//! contact forces
	void ContactForces(FEGlobalVector& R);

	//! the non-linear constraint forces
	void NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp);

	//! nodal forces
	void NodalForces(vector<double>& F, const FETimeInfo& tp);

	//! Inertial forces
	void InertialForces(FEGlobalVector& R);

	//! helper function for setting up the solution phase
	void PrepStep(const FETimeInfo& timeInfo);

	//! modified linesearch for Hager-Zhang solver
	double LineSearchCG(double s);

	//! do an augmentation
	bool Augment();

public:
	double	m_Dtol;
	double	m_Etol;
	double	m_Rtol;
	double	m_Rmin;
	double	m_LStol;
	double	m_LSmin;
	int		m_LSiter;

	// Newmark parameters (for dynamic analyses)
	double	m_beta;			//!< Newmark parameter beta (displacement integration)
	double	m_gamma;		//!< Newmark parameter gamme (velocity integration)

private:
	vector<double>	m_R0;
	vector<double>	m_R1;
	vector<double>	m_Ui;
	vector<double>	m_ui;
	vector<double>	m_Ut;
	vector<double>	m_Fn;
	vector<double>	m_Fd;
	vector<double>	m_Fr;

	int				m_neq;
	int				m_nreq;
	bool		m_baugment;

protected:
	int		m_dofX;
	int		m_dofY;
	int		m_dofZ;
	int		m_dofVX;
	int		m_dofVY;
	int		m_dofVZ;
	int		m_dofU;
	int		m_dofV;
	int		m_dofW;
	int		m_dofRU;
	int		m_dofRV;
	int		m_dofRW;

	DECLARE_PARAMETER_LIST();
};
