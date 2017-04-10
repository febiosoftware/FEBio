#pragma once
#include "FECore/FESolidDomain.h"
#include "FEFluidDomain.h"
#include "FEFluid.h"

//-----------------------------------------------------------------------------
//! domain described by 3D volumetric elements
//!
class FEFluidDomain3D : public FESolidDomain, public FEFluidDomain
{
public:
    //! constructor
    FEFluidDomain3D(FEModel* pfem);
    ~FEFluidDomain3D() {}
    
    //! assignment operator
    FEFluidDomain3D& operator = (FEFluidDomain3D& d);
    
    //! initialize class
    bool Initialize();
    
    //! initialize elements
    void PreSolveUpdate(const FETimeInfo& timeInfo);
    
public: // overrides from FEDomain
    
    //! get the material
    FEMaterial* GetMaterial() { return m_pMat; }
    
    //! set the material
    void SetMaterial(FEMaterial* pm);
    
public: // overrides from FEElasticDomain
    
    // update stresses
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
    FEFluid*	m_pMat;
    
protected:
    int	m_dofVX, m_dofVY, m_dofVZ;
    int	m_dofE;
    int m_dofEP;
    int m_dofAE;
    int m_dofAEP;
};
