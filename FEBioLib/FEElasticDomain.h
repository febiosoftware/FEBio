#pragma once
#include "FECore/FENLSolver.h"
#include "FECore/FEModel.h"
#include "FECore/FEBodyForce.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
//! Abstract interface class for elastic domains.

//! An elastic domain is used by the structural mechanics solver.
//! This interface defines the functions that have to be implemented by an
//! elastic domain. There are basically two categories: residual functions
//! that contribute to the global residual vector. And stiffness matrix 
//! function that calculate contributions to the global stiffness matrix.
class FEElasticDomain
{
public:
	virtual ~FEElasticDomain(){}

	//! Updates the element stresses
	virtual void UpdateStresses(FEModel& fem) = 0;

	// --- R E S I D U A L ---

	// calculate total residual (TODO: remove this)
//	virtual void Residual(FENLSolver* psolver, vector<double>& R) = 0;

	//! calculate the internal forces
	virtual void InternalForces(FENLSolver* psolver, vector<double>& R) = 0;

	//! Calculate the body force vector
	virtual void BodyForce(FENLSolver* psolver, FEBodyForce& bf, vector<double>& R) = 0;

	//! calculate the interial forces (for dynamic problems)
	virtual void InertialForces(FENLSolver* psolver, vector<double>& R, vector<double>& F) = 0;

	// --- S T I F F N E S S   M A T R I X ---

	//! Calculate global stiffness matrix (only contribution from internal force derivative)
	// TODO: maybe I should rename this the InternalStiffness matrix?
	virtual void StiffnessMatrix   (FENLSolver* psolver) = 0;

	//! Calculate stiffness contribution of body forces
	virtual void BodyForceStiffness(FENLSolver* psolver) = 0;

	//! calculate the inertial stiffness matrix (for dynamic problems)
	virtual void InertialStiffness (FENLSolver* psolver) = 0;
};
