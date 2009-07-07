// FESolver.h: interface for the FESolver class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FESOLVER_H__EFCA10FB_7487_44A1_A81D_8B53BE3BECEA__INCLUDED_)
#define AFX_FESOLVER_H__EFCA10FB_7487_44A1_A81D_8B53BE3BECEA__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/LinearSolver.h"
#include "FEStiffnessMatrix.h"
#include "fem.h"
#include "Timer.h"
#include "FEException.h"
#include "Interrupt.h"

//-----------------------------------------------------------------------------
//! The solver class.

//! This class is responsible for solving a single time step in the FE
//! analysis, possibly using results from previous solutions.
//! In the future this class might become a base-class for different FE solvers.
//! For instances there could be a different solver for quasi-static, dynamic and
//! eigenvalue problems.
//! This class is derived from Interruptable, because we want to be able to capture
//! the CTRL+C interruptions


class FESolver : public Interruptable
{
public:
	//! constructor
	FESolver(FEM& fem);

	//! destructor
	virtual ~FESolver();

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

	//{ --- Solution functions ---

		//! prepares the data for the first QN iteration
		void PrepStep(double time);

		//! Performs a Newton-Raphson iteration
		bool Quasin(double time);

		//! Perform an update
		void Update(vector<double>& ui, double s);

		//! Update Stresses
		void UpdateStresses();

		//! Update poroelstic data
		void UpdatePoro(vector<double>& ui, double s);

		//! Update rigid body data
		void UpdateRigidBodies(vector<double>& ui, double s);

		//! Performs a linesearch
		double LineSearch();

		//! Do a BFGS stiffness update
		bool BFGSUpdate(double s);

		//! Lagrangian augmentation
		bool Augment();

		//! solve the system of equations
		void SolveEquations(vector<double>& x, vector<double>& b);

		//! solve the system of equations
		void SolveEquations(matrix& x, matrix& b);
	//}

	//{ --- Stiffness matrix routines ---

		//! calculates the solid element stiffness matrix
		void ElementStiffness(FESolidElement& el, matrix& ke);

		//! calculates the shell element stiffness matrix
		void ElementStiffness(FEShellElement& el, matrix& ke);

		//! calculates the solid element inertial stiffness matrix
		void ElementInertialStiffness(FESolidElement& el, matrix& ke);

		//! calculates the element biphasic stiffness matrix
		bool ElementPoroStiffness(FESolidElement& el, matrix& ke);

		//! material stiffness component
		void MaterialStiffness(FESolidElement& el, matrix& ke);

		//! material stiffness for UDG hex elements
		void UDGMaterialStiffness(FESolidElement& el, matrix& ke);

		//! geometrical stiffness (i.e. initial stress)
		void GeometricalStiffness(FESolidElement& el, matrix& ke);

		//! geometrical stiffness for UDG hex elements
		void UDGGeometricalStiffness(FESolidElement& el, matrix& ke);

		//! Dilatational stiffness component for nearly-incompressible materials
		void DilatationalStiffness(FESolidElement& elem, matrix& ke);

		//! dilatational stiffness for UDG hex elements
		void UDGDilatationalStiffness(FESolidElement& el, matrix& ke);

		//! hourglass stiffness for UDG hex elements
		void UDGHourglassStiffness(FESolidElement& el, matrix& ke);

		//! Dilatational stiffness component for nearly-incompressible materials
		void DilatationalStiffness(FEShellElement& elem, matrix& ke);

		//! Calculates pressure stiffness
		bool PressureStiffness(FESurfaceElement& el, matrix& ke);

		//! contact stiffness
		void ContactStiffness();

		//! calculates the global stiffness matrix
		bool StiffnessMatrix();

		//! reform the stiffness matrix
		bool ReformStiffness();

		//! recalculates the shape of the stiffness matrix
		bool CreateStiffness(bool breset);

		//! calculate the rigid stiffnes matrices
		void RigidStiffness(vector<int>& en, vector<int>& elm, matrix& ke);

		//! calculates the discrete element stiffness
		void DiscreteElementStiffness();

		//! calculates stiffness contributon of linear constraints
		void LinearConstraintStiffness();
	//}

	//{ --- Residual routines ---

		//! Calculatess external body forces for solid elements
		void BodyForces(FESolidElement& elem, vector<double>& fe);

		//! Calculate extenral body forces for shell elements
		void BodyForces(FEShellElement& el, vector<double>& fe);

		//! Calculates concentrated nodal forces
		void NodalForces(vector<double>& F);

		//! Calculates external pressure forces
		bool PressureForce(FESurfaceElement& el, vector<double>& fe);

		//! Calculates the linear external pressure forces (ie. non-follower forces)
		bool LinearPressureForce(FESurfaceElement& el, vector<double>& fe);

		//! Calculates the internal fluid forces
		bool InternalFluidWork(FESolidElement& elem, vector<double>& fe);

		//! Calculates the internal stress vector for solid elements
		void InternalForces(FESolidElement& el, vector<double>& fe);

		//! Calculates the internal stress vector for enhanced strain hex elements
		void UDGInternalForces(FESolidElement& el, vector<double>& fe);

		//! calculates hourglass forces for the UDG element
		void UDGHourglassForces(FESolidElement& el, vector<double>& fe);

		//! Calculates the internal stress vector for shell elements
		void InternalForces(FEShellElement& el, vector<double>& fe);

		//! Calculate inertial forces for dynamic problems
		void InertialForces(vector<double>& R);

		//! Calculate the contact forces
		void ContactForces(vector<double>& R);

		//! Calculates discrete element forces
		void DiscreteElementForces(vector<double>& R);

		//! Calculates residual
		bool Residual(vector<double>& R);

		//! Calculate linear constraint forces
		void LinearConstraintForces(vector<double>& R);
	//}

protected:
	double HexVolume(FESolidElement& el, int state = 0);
	void AvgCartDerivs(FESolidElement& el, double GX[8], double GY[8], double GZ[8], int state = 0);
	void AvgDefGrad(FESolidElement& el, mat3d& F, double GX[8], double GY[8], double GZ[8]);

public:
	//! serialize data to/from dump file
	void Serialize(Archive& ar);

public:
	Logfile&	m_log;	//!< reference to log file
	FEM&		m_fem;	//!< reference the FE data structure

	LinearSolver*		m_psolver;	//!< the linear solver
	FEStiffnessMatrix*	m_pK;		//!< global stiffness matrix
	FEStiffnessMatrix*	m_pM;		//!< mass matrix

	vector<double> m_Fn;	//!< concentrated nodal force vector
	vector<double> m_Fr;	//!< nodal reaction forces
	vector<double> m_ui;	//!< displacement increment vector
	vector<double> m_Ui;	//!< Total displacement vector for iteration
	vector<double> m_Ut;	//!< Total dispalcement vector at time t (incl all previous timesteps)
	vector<double> m_Fd;	//!< residual correction due to prescribed displacements
	vector<double> m_R0;	//!< residual at iteration i-1
	vector<double> m_R1;	//!< residual at iteration i

	// BFGS update vectors
	matrix			m_V;
	matrix			m_W;
	vector<double>	m_D, m_G, m_H;	//!< temp vectors for calculating BFGS update vectors
	double			m_cmax;			//!< maximum value for the condition number

	// convergence tolerances
	double	m_Rtol;			//!< residual tolerance
	double	m_Dtol;			//!< displacement tolerance
	double	m_Etol;			//!< energy tolerance
	double	m_LStol;		//!< line search tolerance
	double	m_LSmin;		//!< minimum line search step
	int		m_LSiter;		//!< max nr of line search iterations
	int		m_maxups;		//!< max nr of QN iters permitted between stiffness reformations
	int		m_maxref;		//!< max nr of reformations per time step

	// counters
	int		m_nrhs;			//!< nr of right hand side evalutations
	int		m_niter;		//!< nr of quasi-newton iterations
	int		m_nref;			//!< nr of stiffness retormations
	int		m_nups;			//!< nr of stiffness updates
	int		m_naug;			//!< nr of augmentations

	// convergence norms
	double		m_normRi;	//!< initial residual norm
	double		m_normR1;	//!< current residual norm
	double		m_normEi;	//!< initial energy norm
	double		m_normE1;	//!< current energy norm
	double		m_normEm;	//!< max energy norm
	double		m_normUi;	//!< initial displacement norm
	double		m_normU;	//!< current displacement norm
	double		m_normu;	//!< incremement displacement norm

	// timers
	Timer	m_SolverTime;	//!< tracks time spent in solver

	// matrix reshape flag
	bool	m_breshape;		//!< Matrix reshape flag
};

#endif // !defined(AFX_FESOLVER_H__EFCA10FB_7487_44A1_A81D_8B53BE3BECEA__INCLUDED_)
