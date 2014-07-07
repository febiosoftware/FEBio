#pragma once

#include "FECore/FESolver.h"
#include "FEStiffnessMatrix.h"
#include "FECore/FETypes.h"
#include "FECore/FERigidBody.h"

//-----------------------------------------------------------------------------
//! The FESolidSolver2 class solves large deformation solid mechanics problems
//! It can deal with quasi-static and dynamic problems
//! 
class FESolidSolver2 : public FESolver
{
public:
	//! constructor
	FESolidSolver2(FEModel* pfem);

	//! destructor
	virtual ~FESolidSolver2();

	//! Clean up
	virtual void Clean();

	//! serialize data to/from dump file
	void Serialize(DumpFile& ar);

	//! Initializes data structures
	bool Init();

	//! solves a single time step
	bool SolveStep(double time);

	//! Initialize linear equation system
	bool InitEquations();

public:
	//! assemble the element residual into the global residual
//	void AssembleResidual(vector<int>& en, vector<int>& elm, vector<double>& fe, vector<double>& R);

	//! adjust the residual matrix for prescribed displacements
	void AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke);

	//! assemble global stiffness matrix \todo this is only used by rigid joints
	void AssembleStiffness(vector<int>& elm, matrix& ke);

public:
	//{ --- evaluation and update ---
		//! Perform an update
		void Update(vector<double>& ui);

		//! Evaluate system, i.e. calculate residual
		void Evaluate(vector<double>& R) { Residual(R); }
	//}

	//{ --- Solution functions ---

		//! prepares the data for the first QN iteration
		virtual void PrepStep(double time);

		//! Performs a Newton-Raphson iteration
		virtual bool Quasin(double time);

		//! update nodal positions, velocities, accelerations, etc.
		virtual void UpdateKinematics(vector<double>& ui);

    //! update DOF increments
        virtual void UpdateIncrements(vector<double>& Ui, vector<double>& ui, bool emap);
    
		//! Update Stresses
		void UpdateStresses();

		//! Update rigid body data
		void UpdateRigidBodies(vector<double>& ui);

		//! update contact data
		virtual void UpdateContact();

		//! update constraint data
		virtual void UpdateConstraints();

		//! Lagrangian augmentation
		bool Augment();
	//}

	//{ --- Stiffness matrix routines ---
		//! return pointer to stiffness matrix
		FEStiffnessMatrix* GetStiffnessMatrix() { return m_pK; }

		//! recalculates the shape of the stiffness matrix
		bool CreateStiffness(bool breset);

		//! calculates the global stiffness matrix
		virtual bool StiffnessMatrix(const FETimePoint& tp);

		//! contact stiffness
		void ContactStiffness();

		//! reform the stiffness matrix
		bool ReformStiffness();

		//! calculate the rigid stiffnes matrices
		void RigidStiffness(vector<int>& en, vector<int>& elm, matrix& ke);

		//! calculates stiffness contributon of nonlinear constraints
		void NonLinearConstraintStiffness();
	//}

	//{ --- Residual routines ---

		//! Calculates concentrated nodal forces
		void NodalForces(vector<double>& F);

		//! Calculate inertial forces for dynamic problems
		void InertialForces(FEGlobalVector& R);

		//! Calculate the contact forces
		void ContactForces(FEGlobalVector& R);

		//! Calculates residual
		virtual bool Residual(vector<double>& R);

		//! Calculate nonlinear constraint forces
		void NonLinearConstraintForces(FEGlobalVector& R);
	//}

	//{ --- Rigid body routines ---

		//! calculate rigid body mass matrix
		void RigidMassMatrix(FERigidBody& RB);

		//! calculate rigid inertial forces
		void RigidInertialForces(FERigidBody& RB, FEGlobalVector& R);
	//}

public:
	// convergence tolerances
	double	m_Rtol;			//!< residual tolerance
	double	m_Dtol;			//!< displacement tolerance
	double	m_Etol;			//!< energy tolerance
	double	m_Rmin;			//!< min residual value

	// equation numbers
	int		m_nreq;			//!< start of rigid body equations

public:
	vector<double> m_Fn;	//!< concentrated nodal force vector
	vector<double> m_Fr;	//!< nodal reaction forces
	vector<double> m_Ui;	//!< Total displacement vector for iteration
	vector<double> m_Ut;	//!< Total dispalcement vector at time t (incl all previous timesteps)
	vector<double> m_Fd;	//!< residual correction due to prescribed displacements

	// matrix reshape flag
	bool	m_breshape;		//!< Matrix reshape flag

    // Newmark parameters (for dynamic analyses)
	double	m_alpha;		//!< Newmark parameter alpha (force integration)
	double	m_beta;			//!< Newmark parameter beta (displacement integration)
	double	m_gamma;		//!< Newmark parameter gamme (velocity integration)
    
public:
	LinearSolver*		m_plinsolve;	//!< the linear solver

	// global stiffness matrix
	FEStiffnessMatrix*	m_pK;		//!< global stiffness matrix
	int					m_neq;		//!< number of equations

	BFGSSolver	m_bfgs;			//!< BFGS solver parameters

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
