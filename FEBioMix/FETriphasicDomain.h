#pragma once
#include "FECore/FESolidDomain.h"
#include "FECore/FETypes.h"
#include "FETriphasic.h"
#include "FEBioMech/FEElasticDomain.h"

//-----------------------------------------------------------------------------
//! Domain class for triphasic 3D solid elements
//! Note that this class inherits from FEElasticSolidDomain since this domain
//! also needs to calculate elastic stiffness contributions.
//!
class FETriphasicDomain : public FESolidDomain, public FEElasticDomain
{
public:
	//! constructor
	FETriphasicDomain(FEModel* pfem);
	
	//! reset domain data
	void Reset() override;

	//! get material (overridden from FEDomain)
	FEMaterial* GetMaterial() override { return m_pMat; }

	//! set the material
	void SetMaterial(FEMaterial* pmat) override;

	//! Unpack solid element data (overridden from FEDomain)
	void UnpackLM(FEElement& el, vector<int>& lm) override;

	//! initialize class
	bool Init() override;

	//! activate
	void Activate() override;

	//! initialize elements for this domain
	void PreSolveUpdate(const FETimeInfo& timeInfo) override;
	
	// update domain data
	void Update(const FETimeInfo& tp) override;

	//! update element state data
	void UpdateElementStress(int iel);

public:
	// internal work (overridden from FEElasticDomain)
	void InternalForces(FEGlobalVector& R) override;

    // internal work (steady-state case)
    void InternalForcesSS(FEGlobalVector& R);
    
	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolver* psolver, bool bsymm);

	//! calculates the global stiffness matrix for this domain (steady-state case)
	void StiffnessMatrixSS(FESolver* psolver, bool bsymm);

protected:
	//! element internal force vector
	void ElementInternalForce(FESolidElement& el, vector<double>& fe);

    //! element internal force vector (steady-state case)
    void ElementInternalForceSS(FESolidElement& el, vector<double>& fe);
    
	//! calculates the element triphasic stiffness matrix
	bool ElementTriphasicStiffness(FESolidElement& el, matrix& ke, bool bsymm);
	
	//! calculates the element triphasic stiffness matrix
	bool ElementTriphasicStiffnessSS(FESolidElement& el, matrix& ke, bool bsymm);
	
protected: // overridden from FEElasticDomain, but not implemented in this domain
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override {}
	void InertialForces(FEGlobalVector& R, vector<double>& F) override {}
	void StiffnessMatrix(FESolver* psolver) override {}
	void BodyForceStiffness(FESolver* psolver, FEBodyForce& bf) override {}
	void MassMatrix(FESolver* psolver, double scale) override {}
	
protected:
	FETriphasic*	m_pMat;
	int				m_dofP;		//!< pressure dof index
	int				m_dofC;		//!< concentration dof index
};
