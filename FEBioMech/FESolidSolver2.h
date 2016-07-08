#pragma once

#include "FECore/FENewtonSolver.h"
#include "FECore/FETypes.h"
#include "FECore/FERigidBody.h"
#include "FECore/FEGlobalVector.h"
#include "FERigidSolver.h"

//-----------------------------------------------------------------------------
//! The FESolidSolver2 class solves large deformation solid mechanics problems
//! It can deal with quasi-static and dynamic problems
//! 
class FESolidSolver2 : public FENewtonSolver
{
public:
	//! constructor
	FESolidSolver2(FEModel* pfem);

	//! destructor
	virtual ~FESolidSolver2();

	//! serialize data to/from dump file
	void Serialize(DumpStream& ar);

	//! Initializes data structures
	bool Init();

	//! initialize the step
	bool InitStep(double time);

	//! Initialize linear equation system
	bool InitEquations();

public:
	//! assemble the element residual into the global residual
	//! \todo This was implemented for nodal forces
	void AssembleResidual(int node, int dof, double f, vector<double>& R);

	//! adjust the residual matrix for prescribed displacements
	void AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke);

	//! assemble global stiffness matrix \todo this is only used by rigid joints
	void AssembleStiffness(vector<int>& elm, matrix& ke);

	//! adjust the residual matrix for prescribed displacements
	void AssembleStiffness2(vector<int>& lmi, vector<int>& lmj, matrix& ke);


public:
	//{ --- evaluation and update ---
		//! Perform an update
		void Update(vector<double>& ui);

		//! Evaluate system, i.e. calculate residual
		void Evaluate(vector<double>& R) { Residual(R); }
	//}

	//{ --- Solution functions ---

		//! prepares the data for the first QN iteration
		virtual void PrepStep(const FETimeInfo& timeInfo);

		//! Performs a Newton-Raphson iteration
		virtual bool Quasin(double time);

		//! update nodal positions, velocities, accelerations, etc.
		virtual void UpdateKinematics(vector<double>& ui);

    //! update DOF increments
        virtual void UpdateIncrements(vector<double>& Ui, vector<double>& ui, bool emap);
    
		//! Update Stresses
		void UpdateStresses();

		//! update contact data
		virtual void UpdateContact();

		//! update constraint data
		virtual void UpdateConstraints();

		//! Lagrangian augmentation
		bool Augment();
	//}

	//{ --- Stiffness matrix routines ---

		//! calculates the global stiffness matrix
		virtual bool StiffnessMatrix(const FETimeInfo& tp);

		//! contact stiffness
		void ContactStiffness();

		//! calculate the rigid stiffnes matrices
		void RigidStiffness(vector<int>& en, vector<int>& elm, matrix& ke);

		//! calculates stiffness contributon of nonlinear constraints
		void NonLinearConstraintStiffness(const FETimeInfo& tp);
	//}

	//{ --- Residual routines ---

		//! Calculates concentrated nodal forces
		// NOTE: I made this function virtual so that derived class (i.e. the bi/multi-phasic solvers)
		//       can handle applied pressure and concentration "forces". But I really want to get rid 
		//       of this function eventually.
		virtual void NodalForces(vector<double>& F, const FETimeInfo& tp);

		//! Calculate inertial forces for dynamic problems
		void InertialForces(FEGlobalVector& R);

		//! Calculate the contact forces
		void ContactForces(FEGlobalVector& R);

		//! Calculates residual
		virtual bool Residual(vector<double>& R);

		//! Calculate nonlinear constraint forces
		void NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp);
	//}

public:
	// convergence tolerances
	double	m_Rtol;			//!< residual tolerance
	double	m_Dtol;			//!< displacement tolerance
	double	m_Etol;			//!< energy tolerance
	double	m_Rmin;			//!< min residual value

	// strategy parameters
	bool	m_bdivreform;	//!< reform when diverging
	bool	m_bdoreforms;	//!< do reformations

	// equation numbers
	int		m_nreq;			//!< start of rigid body equations

public:
	vector<double> m_Fn;	//!< concentrated nodal force vector
	vector<double> m_Fr;	//!< nodal reaction forces
	vector<double> m_Ui;	//!< Total displacement vector for iteration
	vector<double> m_Ut;	//!< Total dispalcement vector at time t (incl all previous timesteps)
	vector<double> m_Fd;	//!< residual correction due to prescribed displacements

    // Newmark parameters (for dynamic analyses)
	double	m_alpha;		//!< Newmark parameter alpha (force integration)
	double	m_beta;			//!< Newmark parameter beta (displacement integration)
	double	m_gamma;		//!< Newmark parameter gamme (velocity integration)
    
public:
	bool		m_baugment;		//!< augmentation flag

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

protected:
	FERigidSolverNew	m_rigidSolver;

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
