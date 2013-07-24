#pragma once
#include <FECore\FESolver.h>
#include <FECore\FEGlobalVector.h>

//-----------------------------------------------------------------------------
//! This class implements a nonlinear explicit solver for solid mechanics
//! problems.
class FEExplicitSolidSolver : public FESolver
{
public:
	//! constructor
	FEExplicitSolidSolver(FEModel& fem);

	//! destructor
	virtual ~FEExplicitSolidSolver() {}

public:
	//! Data initialization
	bool Init();

	//! clean up
	void Clean();

	//! Solve an analysis step
	bool SolveStep(double time);

	//! Update data
	void Update(vector<double>& ui);

	//! Serialize data
	void Serialize(DumpFile& ar);

public:
	//! assemble the element residual into the global residual
//	void AssembleResidual(vector<int>& en, vector<int>& elm, vector<double>& fe, vector<double>& R);

public:

	// initialize equations
	bool InitEquations();

	//! update kinematics
	void UpdateKinematics(vector<double>& ui);

	//! Update rigid bodies 
	void UpdateRigidBodies(vector<double>& ui);

	//! Update stresses
	void UpdateStresses();

	//! solve the step
	bool DoSolve(double time);

	void PrepStep(double time);

	void NodalForces(vector<double>& F);

	bool Residual(vector<double>& R);

	void NonLinearConstraintForces(FEGlobalVector& R);

	void InertialForces(FEGlobalVector& R);
	
	void ContactForces(FEGlobalVector& R);

private:
	//! \todo I have to overload this but I need to remove this.
	virtual void AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke) { assert(false); }

public:
	double		m_dyn_damping;	//!< velocity damping for the explicit solver

public:
	// equation numbers
	int		m_neq;			//!< number of equations
	int		m_nreq;			//!< start of rigid body equations

	vector<double> m_inv_mass;	//!< inverse mass vector for explicit analysis
	vector<double> m_Fn;	//!< concentrated nodal force vector
	vector<double> m_Fr;	//!< nodal reaction forces
	vector<double> m_Ui;	//!< Total displacement vector for iteration
	vector<double> m_Ut;	//!< Total dispalcement vector at time t (incl all previous timesteps)
	vector<double> m_Fd;	//!< residual correction due to prescribed displacements

	vector<double> m_ui;	//!< displacement increment vector

	vector<double> m_R0;	//!< residual at iteration i-1
	vector<double> m_R1;	//!< residual at iteration i

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
