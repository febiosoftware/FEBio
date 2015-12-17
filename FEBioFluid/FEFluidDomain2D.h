//
//  FEFluidDomain2D.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 12/15/15.
//  Copyright © 2015 febio.org. All rights reserved.
//

#pragma once
#include "FECore/FEShellDomain.h"
#include "FECore/FEModel.h"
#include "FEFluidDomain.h"
#include "FEFluid.h"

//-----------------------------------------------------------------------------
//! domain described by 3D volumetric elements
//!
class FEFluidDomain2D : public FEShellDomain, public FEFluidDomain
{
public:
    //! constructor
    FEFluidDomain2D(FEModel* pfem);
    ~FEFluidDomain2D() {}
    
    //! assignment operator
    FEFluidDomain2D& operator = (FEFluidDomain2D& d);
    
    //! initialize class
    bool Initialize(FEModel& fem);
    
    //! activate
    void Activate();
    
    //! initialize elements
    virtual void InitElements();
    
    //! Unpack solid element data
    void UnpackLM(FEElement& el, vector<int>& lm);
    
public: // overrides from FEDomain
    
    //! get the material
    FEMaterial* GetMaterial() { return m_pMat; }
    
    //! set the material
    void SetMaterial(FEMaterial* pm);
    
    //! create a copy
    FEDomain* Copy();
    
public: // overrides from FEElasticDomain
    
    // update stresses
    void UpdateStresses(FEModel& fem);
    
    // update the element stress
    void UpdateElementStress(int iel, double dt);
    
    //! internal stress forces
    void InternalForces(FEGlobalVector& R);
    
    //! body forces
    void BodyForce(FEGlobalVector& R, FEBodyForce& BF);
    
    //! intertial forces for dynamic problems
    void InertialForces(FEGlobalVector& R);
    
    //! calculates the global stiffness matrix for this domain
    void StiffnessMatrix(FESolver* psolver);
    
    //! calculates inertial stiffness
    void MassMatrix(FESolver* psolver);
    
    //! body force stiffness
    void BodyForceStiffness(FESolver* psolver, FEBodyForce& bf);
    
public:
    // --- S T I F F N E S S ---
    
    //! calculates the solid element stiffness matrix
    void ElementStiffness(FEModel& fem, int iel, matrix& ke);
    
    //! material stiffness component
    void ElementMaterialStiffness(FEShellElement& el, matrix& ke);
    
    //! calculates the solid element mass matrix
    void ElementMassMatrix(FEShellElement& el, matrix& ke);
    
    //! calculates the stiffness matrix due to body forces
    void ElementBodyForceStiffness(FEBodyForce& bf, FEShellElement& el, matrix& ke);
    
    // --- R E S I D U A L ---
    
    //! Calculates the internal stress vector for solid elements
    void ElementInternalForce(FEShellElement& el, vector<double>& fe);
    
    //! Calculatess external body forces for solid elements
    void ElementBodyForce(FEBodyForce& BF, FEShellElement& elem, vector<double>& fe);
    
    //! Calculates the inertial force vector for solid elements
    void ElementInertialForce(FEShellElement& el, vector<double>& fe);
    
    // ---
    
protected:
    FEFluid*	m_pMat;

protected:
	int	m_dofVX, m_dofVY;
	int	m_dofE;
};
