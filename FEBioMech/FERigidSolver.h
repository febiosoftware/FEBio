#pragma once
#include <FECore/FETypes.h>
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
class matrix;
class FEModel;
class SparseMatrix;
class FEGlobalVector;
class FERigidBody;
class FESolidSolver2;

//-----------------------------------------------------------------------------
//! This is a helper class that helps the solid deformables solvers update the 
//! state of the rigid system.
class FERigidSolver
{
public:
	FERigidSolver(FEModel* fem);

	// initialize the equation
	// neq is the number of equations that already have been assigned
	// returns the new total number of equations or -1 on error
	int InitEquations(int neq);

	// This is called at the start of each time step
	void PrepStep(const FETimeInfo& timeInfo, vector<double>& ui);

	// update rigid DOF increments
	void UpdateIncrements(vector<double>& Ui, vector<double>& ui, bool emap);

	// update rigid bodies
	void UpdateRigidBodies(vector<double> &Ui, vector<double>& ui);

	// correct stiffness matrix for rigid bodies
	void RigidStiffness(SparseMatrix& K, vector<double>& ui, vector<double>& F, vector<int>& en, vector<int>& elm, matrix& ke, double alpha);

	// adjust residual for rigid bodies
	void AssembleResidual(int node_id, int dof, double f, vector<double>& R);

	// this is called during residual evaluation
	// Currently, this is used for resetting rigid body forces
	void Residual();

	// evaluate inertial data
	void InertialForces(FEGlobalVector& R, const FETimeInfo& timeInfo, double beta, double gamma);

	// contribution from rigid bodies to stiffness matrix
	void StiffnessMatrix(FESolidSolver2* solver, const FETimeInfo& tp);

	// calculate contribution to mass matrix from a rigid body
	void RigidMassMatrix(FESolidSolver2* solver, FERigidBody& RB, const FETimeInfo& timeInfo);

private:
	FEModel*	m_fem;
	int			m_dofX, m_dofY, m_dofZ;
	int			m_dofVX, m_dofVY, m_dofVZ;
};
