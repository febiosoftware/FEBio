//
//  FEFluidDomain2D.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 12/15/15.
//  Copyright © 2015 febio.org. All rights reserved.
//

#pragma once
#include "FECore/FEDomain2D.h"
#include "FEFluidDomain.h"
#include "FEFluid.h"

//-----------------------------------------------------------------------------
//! domain described by 2D elements
//!
class FEFluidDomain2D : public FEDomain2D, public FEFluidDomain
{
public:
    //! constructor
    FEFluidDomain2D(FEModel* pfem);
    ~FEFluidDomain2D() {}
    
    //! assignment operator
    FEFluidDomain2D& operator = (FEFluidDomain2D& d);
    
    //! initialize class
	bool Init() override;
    
    //! initialize elements
    void PreSolveUpdate(const FETimeInfo& timeInfo);
    
public: // overrides from FEDomain
    
    //! get the material
    FEMaterial* GetMaterial() { return m_pMat; }
    
    //! set the material
    void SetMaterial(FEMaterial* pm);
    
public: // overrides from FEElasticDomain
    
    // update domain data
    void Update(const FETimeInfo& tp);
    
    // update the element stress
    void UpdateElementStress(int iel, const FETimeInfo& tp);
    
    //! internal stress forces
    void InternalForces(FEGlobalVector& R, const FETimeInfo& tp);
    
    //! body forces
    void BodyForce(FEGlobalVector& R, const FETimeInfo& tp, FEBodyForce& BF);
    
    //! intertial forces for dynamic problems
    void InertialForces(FEGlobalVector& R, const FETimeInfo& tp);
    
    //! calculates the global stiffness matrix for this domain
    void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp);
    
    //! calculates inertial stiffness
    void MassMatrix(FESolver* psolver, const FETimeInfo& tp);
    
    //! body force stiffness
    void BodyForceStiffness(FESolver* psolver, const FETimeInfo& tp, FEBodyForce& bf);
    
public:
    // --- S T I F F N E S S ---
    
    //! calculates the solid element stiffness matrix
    void ElementStiffness(int iel, matrix& ke);
    
    //! material stiffness component
    void ElementMaterialStiffness(FEElement2D& el, matrix& ke);
    
    //! calculates the solid element mass matrix
    void ElementMassMatrix(FEElement2D& el, matrix& ke);
    
    //! calculates the stiffness matrix due to body forces
    void ElementBodyForceStiffness(FEBodyForce& bf, FEElement2D& el, matrix& ke);
    
    // --- R E S I D U A L ---
    
    //! Calculates the internal stress vector for solid elements
    void ElementInternalForce(FEElement2D& el, vector<double>& fe);
    
    //! Calculatess external body forces for solid elements
    void ElementBodyForce(FEBodyForce& BF, FEElement2D& elem, vector<double>& fe);
    
    //! Calculates the inertial force vector for solid elements
    void ElementInertialForce(FEElement2D& el, vector<double>& fe);
    
    // ---
    
protected:
    FEFluid*	m_pMat;

protected:
    int	m_dofWX, m_dofWY, m_dofWZ;
    int	m_dofEF;
    int	m_dofWXP, m_dofWYP, m_dofWZP;
    int m_dofEFP;
    int	m_dofAWX, m_dofAWY, m_dofAWZ;
    int m_dofAEF;
    int	m_dofAWXP, m_dofAWYP, m_dofAWZP;
    int m_dofAEFP;
};
