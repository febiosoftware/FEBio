#pragma once
#include "FECore/FESolidDomain.h"
#include "FECore/FEModel.h"
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
    bool Initialize(FEModel& fem);
    
    //! activate
    void Activate();
    
    //! initialize elements
    void InitElements();
    
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
    void Update();
    
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
    void ElementMaterialStiffness(FESolidElement& el, matrix& ke);
    
    //! calculates the solid element mass matrix
    void ElementMassMatrix(FESolidElement& el, matrix& ke);
    
    //! calculates the stiffness matrix due to body forces
    void ElementBodyForceStiffness(FEBodyForce& bf, FESolidElement& el, matrix& ke);
    
    // --- R E S I D U A L ---
    
    //! Calculates the internal stress vector for solid elements
    void ElementInternalForce(FESolidElement& el, vector<double>& fe);
    
    //! Calculatess external body forces for solid elements
    void ElementBodyForce(FEBodyForce& BF, FESolidElement& elem, vector<double>& fe);
    
    //! Calculates the inertial force vector for solid elements
    void ElementInertialForce(FESolidElement& el, vector<double>& fe);
    
    // ---
    
protected:
    FEFluid*	m_pMat;
    
protected:
    int	m_dofVX, m_dofVY, m_dofVZ;
    int	m_dofE;
};
