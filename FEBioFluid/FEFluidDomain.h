//
//  FEFluidDomain.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 12/16/15.
//  Copyright © 2015 febio.org. All rights reserved.
//
#pragma once
#include "FEBioMech/FEBodyForce.h"
#include "FECore/FESolver.h"
#include "FECore/FEModel.h"
#include <vector>
using namespace std;

class FEModel;

//-----------------------------------------------------------------------------
//! Abstract interface class for fluid domains.

//! An fluid domain is used by the fluid mechanics solver.
//! This interface defines the functions that have to be implemented by a
//! fluid domain. There are basically two categories: residual functions
//! that contribute to the global residual vector. And stiffness matrix
//! function that calculate contributions to the global stiffness matrix.
class FEFluidDomain
{
public:
    FEFluidDomain(FEModel* pfem);
    virtual ~FEFluidDomain(){}
    
    //! Updates the element stresses
    virtual void UpdateStresses(FEModel& fem) = 0;
    
    // --- R E S I D U A L ---
    
    //! calculate the internal forces
    virtual void InternalForces(FEGlobalVector& R) = 0;
    
    //! Calculate the body force vector
    virtual void BodyForce(FEGlobalVector& R, FEBodyForce& bf) = 0;
    
    //! calculate the interial forces (for dynamic problems)
    virtual void InertialForces(FEGlobalVector& R) = 0;
    
    // --- S T I F F N E S S   M A T R I X ---
    
    //! Calculate global stiffness matrix (only contribution from internal force derivative)
    //! \todo maybe I should rename this the InternalStiffness matrix?
    virtual void StiffnessMatrix   (FESolver* psolver) = 0;
    
    //! Calculate stiffness contribution of body forces
    virtual void BodyForceStiffness(FESolver* psolver, FEBodyForce& bf) = 0;
    
    //! calculate the mass matrix (for dynamic problems)
    virtual void MassMatrix(FESolver* psolver) = 0;
    
    //! transient analysis
    void SetTransientAnalysis() { m_btrans = true; }
    void SetSteadyStateAnalysis() { m_btrans = false; }
    
public:
    FEModel* GetFEModel() { return m_pfem; }
    
protected:
    FEModel*	m_pfem;
    bool        m_btrans;   // flag for transient (true) or steady-state (false) analysis
};
