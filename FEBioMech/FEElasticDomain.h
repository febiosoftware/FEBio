#pragma once
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
class FEModel;
class FEGlobalVector;
class FEBodyForce;
class FESolver;

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
	FEElasticDomain(FEModel* pfem);
	virtual ~FEElasticDomain(){}

	// --- R E S I D U A L ---

	//! calculate the internal forces
	virtual void InternalForces(FEGlobalVector& R) = 0;

	//! Calculate the body force vector
	virtual void BodyForce(FEGlobalVector& R, FEBodyForce& bf) = 0;

	//! calculate the interial forces (for dynamic problems)
	virtual void InertialForces(FEGlobalVector& R, vector<double>& F) = 0;

	//! calculate the interial forces (used by FESolidSolver2)
	virtual void InertialForces2(FEGlobalVector& R, vector<double>& F) {}

	// --- S T I F F N E S S   M A T R I X ---

	//! Calculate global stiffness matrix (only contribution from internal force derivative)
	//! \todo maybe I should rename this the InternalStiffness matrix?
	virtual void StiffnessMatrix   (FESolver* psolver) = 0;

	//! Calculate stiffness contribution of body forces
	virtual void BodyForceStiffness(FESolver* psolver, FEBodyForce& bf) = 0;

	//! calculate the mass matrix (for dynamic problems)
	virtual void MassMatrix(FESolver* psolver, double scale) = 0;

public:
	FEModel* GetFEModel() { return m_pfem; }

protected:
	FEModel*	m_pfem;
	int					m_dofX;		//!< X-dof index
	int					m_dofY;		//!< Y-dof index
	int					m_dofZ;		//!< Z-dof index
    int					m_dofU;
    int					m_dofV;
    int					m_dofW;
    int                 m_dofUP;
    int                 m_dofVP;
    int                 m_dofWP;
    int					m_dofRU;
	int					m_dofRV;	
	int					m_dofRW;

    int                 m_dofVX;
    int                 m_dofVY;
    int                 m_dofVZ;
    
    int                 m_dofVU;
    int                 m_dofVV;
    int                 m_dofVW;
    int                 m_dofVUP;
    int                 m_dofVVP;
    int                 m_dofVWP;
    
    int                 m_dofAU;
    int                 m_dofAV;
    int                 m_dofAW;
    int                 m_dofAUP;
    int                 m_dofAVP;
    int                 m_dofAWP;
};
