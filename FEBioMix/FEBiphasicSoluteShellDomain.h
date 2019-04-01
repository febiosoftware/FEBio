#pragma once
#include <FEBioMech/FESSIShellDomain.h>
#include "FEBiphasicSolute.h"
#include "FEBiphasicSoluteDomain.h"

//-----------------------------------------------------------------------------
//! Domain class for biphasic-solute 3D solid elements
//! Note that this class inherits from FEElasticSolidDomain since this domain
//! also needs to calculate elastic stiffness contributions.
//!
class FEBIOMIX_API FEBiphasicSoluteShellDomain : public FESSIShellDomain, public FEBiphasicSoluteDomain
{
public:
    //! constructor
    FEBiphasicSoluteShellDomain(FEModel* pfem);
    
    //! reset domain data
    void Reset() override;
    
    //! get the material (overridden from FEDomain)
    FEMaterial* GetMaterial() override { return m_pMat; }
    
    //! set the material
    void SetMaterial(FEMaterial* pmat) override;
    
    //! Unpack solid element data (overridden from FEDomain)
    void UnpackLM(FEElement& el, vector<int>& lm) override;
    
    //! Activate
    void Activate() override;
    
    //! initialize material points in the domain
    void InitMaterialPoints() override;
    
    //! initialize elements for this domain
    void PreSolveUpdate(const FETimeInfo& timeInfo) override;
    
    // update domain data
    void Update(const FETimeInfo& tp) override;
    
    // update element stress
    void UpdateElementStress(int iel);
    
public:
    // internal work (overridden from FEElasticDomain)
    void InternalForces(FEGlobalVector& R) override;
    
    // internal work (steady-state analyses)
    void InternalForcesSS(FEGlobalVector& R) override;
    
public:
    //! calculates the global stiffness matrix for this domain
    void StiffnessMatrix(FESolver* psolver, bool bsymm) override;
    
    //! calculates the global stiffness matrix for this domain (steady-state case)
    void StiffnessMatrixSS(FESolver* psolver, bool bsymm) override;
    
protected:
    //! element internal force vector
    void ElementInternalForce(FEShellElement& el, vector<double>& fe);
    
    //! element internal force vector (steady-state analyses)
    void ElementInternalForceSS(FEShellElement& el, vector<double>& fe);
    
    //! calculates the element solute-poroelastic stiffness matrix
    bool ElementBiphasicSoluteStiffness(FEShellElement& el, matrix& ke, bool bsymm);
    
    //! calculates the element solute-poroelastic stiffness matrix
    bool ElementBiphasicSoluteStiffnessSS(FEShellElement& el, matrix& ke, bool bsymm);
    
protected: // overridden from FEElasticDomain, but not implemented in this domain
    void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override {}
    void InertialForces(FEGlobalVector& R, vector<double>& F) override {}
    void StiffnessMatrix(FESolver* psolver) override {}
    void BodyForceStiffness(FESolver* psolver, FEBodyForce& bf) override {}
    void MassMatrix(FESolver* psolver, double scale) override {}
    
protected:
    int					m_dofSX;
    int					m_dofSY;
    int					m_dofSZ;
};
