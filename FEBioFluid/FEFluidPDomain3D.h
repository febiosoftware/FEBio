//
//  FEFluidPDomain3D.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 2/16/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#ifndef FEFluidPDomain3D_hpp
#define FEFluidPDomain3D_hpp

#include <FECore/FESolidDomain.h>
#include "FEFluidDomain.h"
#include "FEFluidP.h"

//-----------------------------------------------------------------------------
//! domain described by 3D volumetric elements
//!
class FEBIOFLUID_API FEFluidPDomain3D : public FESolidDomain, public FEFluidDomain
{
public:
    //! constructor
    FEFluidPDomain3D(FEModel* pfem);
    ~FEFluidPDomain3D() {}
    
    //! assignment operator
    FEFluidPDomain3D& operator = (FEFluidPDomain3D& d);
    
    //! initialize class
    bool Init() override;
    
    //! initialize elements
    void PreSolveUpdate(const FETimeInfo& timeInfo) override;
    
public: // overrides from FEDomain
    
    //! get the material
    FEMaterial* GetMaterial() override { return m_pMat; }
    
    //! set the material
    void SetMaterial(FEMaterial* pm) override;
    
public: // overrides from FEElasticDomain
    
    // update stresses
    void Update(const FETimeInfo& tp) override;
    
    // update the element stress
    void UpdateElementStress(int iel, const FETimeInfo& tp);
    
    //! internal stress forces
    void InternalForces(FEGlobalVector& R, const FETimeInfo& tp) override;
    
    //! body forces
    void BodyForce(FEGlobalVector& R, const FETimeInfo& tp, FEBodyForce& BF) override;
    
    //! intertial forces for dynamic problems
    void InertialForces(FEGlobalVector& R, const FETimeInfo& tp) override;
    
    //! calculates the global stiffness matrix for this domain
    void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) override;
    
    //! calculates inertial stiffness
    void MassMatrix(FESolver* psolver, const FETimeInfo& tp) override;
    
    //! body force stiffness
    void BodyForceStiffness(FESolver* psolver, const FETimeInfo& tp, FEBodyForce& bf) override;
    
public:
    // --- S T I F F N E S S ---
    
    //! calculates the solid element stiffness matrix
    void ElementStiffness(FESolidElement& el, matrix& ke, const FETimeInfo& tp);
    
    //! calculates the solid element mass matrix
    void ElementMassMatrix(FESolidElement& el, matrix& ke, const FETimeInfo& tp);
    
    //! calculates the stiffness matrix due to body forces
    void ElementBodyForceStiffness(FEBodyForce& bf, FESolidElement& el, matrix& ke, const FETimeInfo& tp);
    
    // --- R E S I D U A L ---
    
    //! Calculates the internal stress vector for solid elements
    void ElementInternalForce(FESolidElement& el, vector<double>& fe, const FETimeInfo& tp);
    
    //! Calculatess external body forces for solid elements
    void ElementBodyForce(FEBodyForce& BF, FESolidElement& elem, vector<double>& fe, const FETimeInfo& tp);
    
    //! Calculates the inertial force vector for solid elements
    void ElementInertialForce(FESolidElement& el, vector<double>& fe, const FETimeInfo& tp);
    
protected:
    FEFluidP*    m_pMat;
    
protected:
    int         m_dofWX, m_dofWY, m_dofWZ;
    int         m_dofEF;
    int         m_dofWXP, m_dofWYP, m_dofWZP;
    int         m_dofEFP;
    int         m_dofAWX, m_dofAWY, m_dofAWZ;
    int         m_dofAEF;
    int         m_dofAWXP, m_dofAWYP, m_dofAWZP;
    int         m_dofAEFP;
};

#endif /* FEFluidPDomain3D_hpp */
