/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include "FEBodyForce.h"
#include <FECore/FETimeInfo.h>
#include <FECore/FESolver.h>
#include <vector>

//-----------------------------------------------------------------------------
class matrix;
class FEModel;
class SparseMatrix;
class FEGlobalVector;
class FERigidBody;
class FESolidSolver2;
class DumpStream;
class FEElementMatrix;
class FEMechModel;

//-----------------------------------------------------------------------------
//! This is a helper class that helps the solid deformables solvers update the 
//! state of the rigid system.
class FEBIOMECH_API FERigidSolver
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
	void RigidStiffness(SparseMatrix& K, std::vector<double>& ui, std::vector<double>& F, const FEElementMatrix& ke, double alpha);

    // correct stiffness matrix for rigid bodies accounting for rigid-body-deformable-shell interfaces
    void RigidStiffnessSolid(SparseMatrix& K, std::vector<double>& ui, std::vector<double>& F, const std::vector<int>& en, const std::vector<int>& lmi, const std::vector<int>& lmj, const matrix& ke, double alpha);
    
    // correct stiffness matrix for rigid bodies accounting for rigid-body-deformable-shell interfaces
    void RigidStiffnessShell(SparseMatrix& K, std::vector<double>& ui, std::vector<double>& F, const std::vector<int>& en, const std::vector<int>& lmi, const std::vector<int>& lmj, const matrix& ke, double alpha);
    
	// adjust residual for rigid-deformable interface nodes
	void AssembleResidual(int node_id, int dof, double f, std::vector<double>& R);
    
	// this is called during residual evaluation
	// Currently, this is used for resetting rigid body forces
	void Residual();

	// contribution from rigid bodies to stiffness matrix
	void StiffnessMatrix(SparseMatrix& K, const FETimeInfo& tp);

	// calculate contribution to mass matrix from a rigid body
	void RigidMassMatrix(FELinearSystem& LS, const FETimeInfo& timeInfo);

	//! Serialization
	void Serialize(DumpStream& ar);

	//! return the FEModel
	FEModel* GetFEModel();

public:
	void AllowMixedBCs(bool b) { m_bAllowMixedBCs = b; }

protected:
	FEMechModel*	m_fem;
	int			m_dofX, m_dofY, m_dofZ;
	int			m_dofVX, m_dofVY, m_dofVZ;
    int			m_dofU, m_dofV, m_dofW;
    int         m_dofSX, m_dofSY, m_dofSZ;
    int         m_dofSVX, m_dofSVY, m_dofSVZ;
	bool		m_bAllowMixedBCs;
};

//-----------------------------------------------------------------------------
class FEBIOMECH_API FERigidSolverOld : public FERigidSolver
{
public:
	FERigidSolverOld(FEModel* fem) : FERigidSolver(fem) { AllowMixedBCs(true); }

	//! update rigid body kinematics for dynamic problems
	void UpdateRigidKinematics();

	//! Update rigid body data
	void UpdateRigidBodies(std::vector<double>& Ui, std::vector<double>& ui, bool bnewUpdate);
};

//-----------------------------------------------------------------------------
class FEBIOMECH_API FERigidSolverNew : public FERigidSolver
{
public:
	FERigidSolverNew(FEModel* fem) : FERigidSolver(fem){}

	// evaluate inertial data
	void InertialForces(FEGlobalVector& R, const FETimeInfo& timeInfo);

	// update rigid DOF increments
	void UpdateIncrements(std::vector<double>& Ui, std::vector<double>& ui, bool emap);

	// update rigid bodies
	void UpdateRigidBodies(vector<double> &Ui, vector<double>& ui);
};
