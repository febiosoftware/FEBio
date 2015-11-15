#pragma once
#include "FECore/FESolidDomain.h"
#include "FECore/FEModel.h"
#include "FEFluid.h"

//-----------------------------------------------------------------------------
//! domain described by 3D volumetric elements
//!
class FEFluidDomain : public FESolidDomain
{
public:
    //! constructor
    FEFluidDomain(FEModel* pfem);
    
    //! assignment operator
    FEFluidDomain& operator = (FEFluidDomain& d);
    
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
    
    void SetTransientAnalysis() { m_btrans = true; }
    void SetSteadyStateAnalysis() { m_btrans = false; }
    void Setform(bool form) { m_boldform = form; }
    
    // ---
    
protected:
    FEFluid*	m_pMat;
    bool        m_btrans;   // flag for transient (true) or steady-state (false) analysis
    bool        m_boldform; // flag for using old form of virtual work integral
};
