#pragma once

#include "FESolver.h"

//-----------------------------------------------------------------------------
//! The FESolidSolver class solves large deformation solid mechanics problems
//! It can deal with quasi-static, dynamic, and poro-elastic problems
//! 
class FESolidSolver : public FESolver
{
public:
	//! constructor
	FESolidSolver(FEM& fem);

	//! Initializes data structures
	bool Init();

	//! solves a single time step
	bool SolveStep(double time);

	//! assemble the element residual into the global residual
	void AssembleResidual(vector<int>& en, vector<int>& elm, vector<double>& fe, vector<double>& R);

	//! adjust the residual matrix for prescribed displacements
	void AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke);

public:
	FEStiffnessMatrix* GetStiffnessMatrix() { return m_pK; }

public:
	//{ --- NonLinearSystem overrides ---
		//! Perform an update
		void Update(vector<double>& ui);

		//! Evaluate system, i.e. calculate residual
		void Evaluate(vector<double>& R) { Residual(R); }
	//}

	//{ --- Solution functions ---

		//! prepares the data for the first QN iteration
		void PrepStep(double time);

		//! Performs a Newton-Raphson iteration
		bool Quasin(double time);

		//! Update Stresses
		void UpdateStresses();

		//! Update poroelastic data
		void UpdatePoro(vector<double>& ui);

		//! Update rigid body data
		void UpdateRigidBodies(vector<double>& ui);

		//! Update solute data
		void UpdateSolute(vector<double>& ui);

		//! Lagrangian augmentation
		bool Augment();
	//}

	//{ --- Stiffness matrix routines ---

		//! contact stiffness
		void ContactStiffness();

		//! calculates the global stiffness matrix
		bool StiffnessMatrix();

		//! reform the stiffness matrix
		bool ReformStiffness();

		//! calculate the rigid stiffnes matrices
		void RigidStiffness(vector<int>& en, vector<int>& elm, matrix& ke);

		//! calculates stiffness contributon of linear constraints
		void LinearConstraintStiffness();
	//}

	//{ --- Residual routines ---

		//! Calculates concentrated nodal forces
		void NodalForces(vector<double>& F);

		//! Calculate inertial forces for dynamic problems
		void InertialForces(vector<double>& R);

		//! Calculate the contact forces
		void ContactForces(vector<double>& R);

		//! Calculates residual
		bool Residual(vector<double>& R);

		//! Calculate linear constraint forces
		void LinearConstraintForces(vector<double>& R);
	//}

protected:
	void GetPressureData(vector<double>& pi, vector<double>& ui);

public:
	//! serialize data to/from dump file
	void Serialize(DumpFile& ar);

public:
	vector<double> m_Fn;	//!< concentrated nodal force vector
	vector<double> m_Fr;	//!< nodal reaction forces
	vector<double> m_Ui;	//!< Total displacement vector for iteration
	vector<double> m_Ut;	//!< Total dispalcement vector at time t (incl all previous timesteps)
	vector<double> m_Fd;	//!< residual correction due to prescribed displacements

	// poro data
	vector<double>	m_pi;	//!< pressure increment vector
	vector<double>	m_Pi;	//!< Total pressure vector for iteration

	vector<double>	m_ci;	//!< concentration increment vector
	vector<double>	m_Ci;	//!< Total concentration vector for iteration

	// convergence norms
	double		m_normRi;	//!< initial residual norm
	double		m_normEi;	//!< initial energy norm
	double		m_normEm;	//!< max energy norm
	double		m_normUi;	//!< initial displacement norm

	// poro data
	double		m_normPi;	//!< initial pressure norm
	double		m_normP;	//!< current pressure norm
	double		m_normp;	//!< incremement pressure norm

	// matrix reshape flag
	bool	m_breshape;		//!< Matrix reshape flag
};
