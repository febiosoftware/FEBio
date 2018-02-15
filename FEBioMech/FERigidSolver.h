#pragma once
#include "FEBodyForce.h"
#include <FECore/FETypes.h>
#include <FECore/FESolver.h>
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
class matrix;
class FEModel;
class SparseMatrix;
class FEGlobalVector;
class FERigidBody;
class FESolidSolver2;
class DumpStream;

//-----------------------------------------------------------------------------
//! This is a helper class that helps the solid deformables solvers update the 
//! state of the rigid system.
class FERigidSolver
{
public:
	FERigidSolver(FEModel* fem);

	// destructor
	virtual ~FERigidSolver(){}

	// initialize the equation
	// neq is the number of equations that already have been assigned
	// returns the new total number of equations or -1 on error
	int InitEquations(int neq);

	// This is called at the start of each time step
	void PrepStep(const FETimeInfo& timeInfo, vector<double>& ui);

	// correct stiffness matrix for rigid bodies
	void RigidStiffness(SparseMatrix& K, vector<double>& ui, vector<double>& F, vector<int>& en, vector<int>& elm, matrix& ke, double alpha);

    // correct stiffness matrix for rigid bodies accounting for rigid-body-deformable-shell interfaces
    void RigidStiffnessSolid(SparseMatrix& K, vector<double>& ui, vector<double>& F, vector<int>& en, vector<int>& elm, matrix& ke, double alpha);
    
    // correct stiffness matrix for rigid bodies accounting for rigid-body-deformable-shell interfaces
    void RigidStiffnessShell(SparseMatrix& K, vector<double>& ui, vector<double>& F, vector<int>& en, vector<int>& elm, matrix& ke, double alpha);
    
	// adjust residual for rigid-deformable interface nodes
	void AssembleResidual(int node_id, int dof, double f, vector<double>& R);
    
	// this is called during residual evaluation
	// Currently, this is used for resetting rigid body forces
	void Residual();

	// contribution from rigid bodies to stiffness matrix
	void StiffnessMatrix(SparseMatrix& K, const FETimeInfo& tp);

	// calculate contribution to mass matrix from a rigid body
	void RigidMassMatrix(FESolver* solver, const FETimeInfo& timeInfo);

	//! Serialization
	void Serialize(DumpStream& ar);

public:
	void AllowMixedBCs(bool b) { m_bAllowMixedBCs = b; }

protected:
	FEModel*	m_fem;
	int			m_dofX, m_dofY, m_dofZ;
	int			m_dofVX, m_dofVY, m_dofVZ;
    int			m_dofU, m_dofV, m_dofW;
	bool		m_bAllowMixedBCs;
};

//-----------------------------------------------------------------------------
class FERigidSolverOld : public FERigidSolver
{
public:
	FERigidSolverOld(FEModel* fem) : FERigidSolver(fem) { AllowMixedBCs(true); }

	//! update rigid body kinematics for dynamic problems
	void UpdateRigidKinematics();

	//! Update rigid body data
	void UpdateRigidBodies(vector<double>& Ui, vector<double>& ui, bool bnewUpdate);
};

//-----------------------------------------------------------------------------
class FERigidSolverNew : public FERigidSolver
{
public:
	FERigidSolverNew(FEModel* fem) : FERigidSolver(fem){}

    // evaluate body forces
    void BodyForces(FEGlobalVector& R, const FETimeInfo& timeInfo, FEBodyForce& pbf);
    
    // evaluate inertial data
	void InertialForces(FEGlobalVector& R, const FETimeInfo& timeInfo);

	// update rigid DOF increments
	void UpdateIncrements(vector<double>& Ui, vector<double>& ui, bool emap);

	// update rigid bodies
	void UpdateRigidBodies(vector<double> &Ui, vector<double>& ui);
};
