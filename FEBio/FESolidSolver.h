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

		//! Lagrangian augmentation
		bool Augment();
	//}

	//{ --- Stiffness matrix routines ---

		//! calculates the solid element stiffness matrix
		void ElementStiffness(FESolidElement& el, matrix& ke);

		//! calculates the shell element stiffness matrix
		void ElementStiffness(FEShellElement& el, matrix& ke);

		//! calculates the truss element stiffness matrix
		void ElementStiffness(FETrussElement& el, matrix& ke);

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

		//! Calculates the internal stress vector for solid elements
		void InternalForces(FETrussElement& el, vector<double>& fe);

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

	void GetPressureData(vector<double>& pi, vector<double>& ui);

public:
	//! serialize data to/from dump file
	void Serialize(Archive& ar);

public:
	vector<double> m_Fn;	//!< concentrated nodal force vector
	vector<double> m_Fr;	//!< nodal reaction forces
	vector<double> m_ui;	//!< displacement increment vector
	vector<double> m_Ui;	//!< Total displacement vector for iteration
	vector<double> m_Ut;	//!< Total dispalcement vector at time t (incl all previous timesteps)
	vector<double> m_Fd;	//!< residual correction due to prescribed displacements
	vector<double> m_R0;	//!< residual at iteration i-1
	vector<double> m_R1;	//!< residual at iteration i

	// poro data
	vector<double>	m_pi;	//!< pressure increment vector
	vector<double>	m_Pi;	//!< Total pressure vector for iteration

	// convergence norms
	double		m_normRi;	//!< initial residual norm
	double		m_normR1;	//!< current residual norm
	double		m_normEi;	//!< initial energy norm
	double		m_normE1;	//!< current energy norm
	double		m_normEm;	//!< max energy norm
	double		m_normUi;	//!< initial displacement norm
	double		m_normU;	//!< current displacement norm
	double		m_normu;	//!< incremement displacement norm

	// poro data
	double		m_normPi;	//!< initial pressure norm
	double		m_normP;	//!< current pressure norm
	double		m_normp;	//!< incremement pressure norm

	// matrix reshape flag
	bool	m_breshape;		//!< Matrix reshape flag
};
