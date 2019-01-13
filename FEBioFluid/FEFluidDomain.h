//
//  FEFluidDomain.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 12/16/15.
//  Copyright © 2015 febio.org. All rights reserved.
//
#pragma once
#include <vector>
#include "febiofluid_api.h"
using namespace std;

class FEModel;
class FESolver;
class FEBodyForce;
class FEGlobalVector;
class FETimeInfo;

//-----------------------------------------------------------------------------
//! Abstract interface class for fluid domains.

//! An fluid domain is used by the fluid mechanics solver.
//! This interface defines the functions that have to be implemented by a
//! fluid domain. There are basically two categories: residual functions
//! that contribute to the global residual vector. And stiffness matrix
//! function that calculate contributions to the global stiffness matrix.
class FEBIOFLUID_API FEFluidDomain
{
public:
    FEFluidDomain(FEModel* pfem);
    virtual ~FEFluidDomain(){}
    
    // --- R E S I D U A L ---
    
    //! calculate the internal forces
    virtual void InternalForces(FEGlobalVector& R, const FETimeInfo& tp) = 0;
    
    //! Calculate the body force vector
    virtual void BodyForce(FEGlobalVector& R, const FETimeInfo& tp, FEBodyForce& bf) = 0;
    
    //! calculate the interial forces (for dynamic problems)
    virtual void InertialForces(FEGlobalVector& R, const FETimeInfo& tp) = 0;
    
    // --- S T I F F N E S S   M A T R I X ---
    
    //! Calculate global stiffness matrix (only contribution from internal force derivative)
    //! \todo maybe I should rename this the InternalStiffness matrix?
    virtual void StiffnessMatrix   (FESolver* psolver, const FETimeInfo& tp) = 0;
    
    //! Calculate stiffness contribution of body forces
    virtual void BodyForceStiffness(FESolver* psolver, const FETimeInfo& tp, FEBodyForce& bf) = 0;
    
    //! calculate the mass matrix (for dynamic problems)
    virtual void MassMatrix(FESolver* psolver, const FETimeInfo& tp) = 0;
    
    //! transient analysis
    void SetTransientAnalysis() { m_btrans = true; }
    void SetSteadyStateAnalysis() { m_btrans = false; }
    
protected:
    bool        m_btrans;   // flag for transient (true) or steady-state (false) analysis
};
