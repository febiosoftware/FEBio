//
//  FEBiphasicFSIDomain.hpp
//  FEBioFluid
//
//  Created by Jay Shim on 12/30/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#pragma once
#include "febiofluid_api.h"
#include <vector>
using namespace std;

class FEModel;
class FELinearSystem;
class FEBodyForce;
class FEGlobalVector;
class FETimeInfo;

//-----------------------------------------------------------------------------
//! Abstract interface class for fluid-FSI domains.

//! A fluid-FSI domain is used by the fluid-FSI solver.
//! This interface defines the functions that have to be implemented by a
//! fluid-FSI domain. There are basically two categories: residual functions
//! that contribute to the global residual vector. And stiffness matrix
//! function that calculate contributions to the global stiffness matrix.
class FEBIOFLUID_API FEBiphasicFSIDomain
{
public:
    FEBiphasicFSIDomain(FEModel* pfem);
    virtual ~FEBiphasicFSIDomain(){}
    
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
    virtual void StiffnessMatrix   (FELinearSystem& LS, const FETimeInfo& tp) = 0;
    
    //! Calculate stiffness contribution of body forces
    virtual void BodyForceStiffness(FELinearSystem& LS, const FETimeInfo& tp, FEBodyForce& bf) = 0;
    
    //! calculate the mass matrix (for dynamic problems)
    virtual void MassMatrix(FELinearSystem& LS, const FETimeInfo& tp) = 0;
    
    //! transient analysis
    void SetTransientAnalysis() { m_btrans = true; }
    void SetSteadyStateAnalysis() { m_btrans = false; }
    void SetQuasiAnalysis() { m_bquasi = true; }
    void UnsetQuasiStateAnalysis() { m_bquasi = false; }
    
protected:
    bool        m_btrans;   // flag for transient (true) or steady-state (false) analysis
    bool        m_bquasi;   // flag for quasi (true) or not (false)
};
