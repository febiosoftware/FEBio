#pragma once
#include "FECore/FESolidDomain.h"
#include "FEMultiphasic.h"
#include "FEMultiphasicDomain.h"

//-----------------------------------------------------------------------------
//! Domain class for multiphasic 3D solid elements
//! Note that this class inherits from FEElasticSolidDomain since this domain
//! also needs to calculate elastic stiffness contributions.
//!
class FEBIOMIX_API FEMultiphasicSolidDomain : public FESolidDomain, public FEMultiphasicDomain
{
public:
    //! constructor
    FEMultiphasicSolidDomain(FEModel* pfem);
    
    //! Reset data
    void Reset() override;
    
    //! get the material (overridden from FEDomain)
    FEMaterial* GetMaterial() override { return m_pMat; }
    
    //! set the material
    void SetMaterial(FEMaterial* pmat) override;
    
    //! Unpack solid element data (overridden from FEDomain)
    void UnpackLM(FEElement& el, vector<int>& lm) override;
    
    //! initialize elements for this domain
    void PreSolveUpdate(const FETimeInfo& timeInfo) override;
    
    //! calculates the global stiffness matrix for this domain
    void StiffnessMatrix(FESolver* psolver, bool bsymm) override;
    
    //! calculates the global stiffness matrix for this domain (steady-state case)
    void StiffnessMatrixSS(FESolver* psolver, bool bsymm) override;
    
    //! initialize class
	bool Init() override;
    
    //! activate
    void Activate() override;
    
    //! initialize material points in the domain
    void InitMaterialPoints() override;
    
    // update domain data
    void Update(const FETimeInfo& tp) override;
    
    // update element state data
    void UpdateElementStress(int iel, double dt);
    
public:
    
    // internal work (overridden from FEElasticDomain)
    void InternalForces(FEGlobalVector& R) override;
    
    // internal work (steady-state case)
    void InternalForcesSS(FEGlobalVector& R) override;
    
public:
    //! element internal force vector
    void ElementInternalForce(FESolidElement& el, vector<double>& fe);
    
    //! element internal force vector (steady-state case)
    void ElementInternalForceSS(FESolidElement& el, vector<double>& fe);
    
    //! calculates the element triphasic stiffness matrix
    bool ElementMultiphasicStiffness(FESolidElement& el, matrix& ke, bool bsymm);
    
    //! calculates the element triphasic stiffness matrix
    bool ElementMultiphasicStiffnessSS(FESolidElement& el, matrix& ke, bool bsymm);
    
protected: // overridden from FEElasticDomain, but not implemented in this domain
    void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override {}
    void InertialForces(FEGlobalVector& R, vector<double>& F) override {}
    void StiffnessMatrix(FESolver* psolver) override {}
    void BodyForceStiffness(FESolver* psolver, FEBodyForce& bf) override {}
    void MassMatrix(FESolver* psolver, double scale) override {}
};
