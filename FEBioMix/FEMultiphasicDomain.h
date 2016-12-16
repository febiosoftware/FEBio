#pragma once
#include "FECore/FESolidDomain.h"
#include "FEMultiphasic.h"
#include "FEBioMech/FEElasticDomain.h"
#include "FECore/FETypes.h"

//-----------------------------------------------------------------------------
//! Domain class for multiphasic 3D solid elements
//! Note that this class inherits from FEElasticSolidDomain since this domain
//! also needs to calculate elastic stiffness contributions.
//!
class FEMultiphasicDomain : public FESolidDomain, public FEElasticDomain
{
public:
	//! constructor
	FEMultiphasicDomain(FEModel* pfem);
	
	//! Reset data
	void Reset();

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial() { return m_pMat; }

	//! set the material
	void SetMaterial(FEMaterial* pmat);

	//! Unpack solid element data (overridden from FEDomain)
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! initialize elements for this domain
	void PreSolveUpdate(const FETimeInfo& timeInfo);
	
	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolver* psolver, bool bsymm);

	//! calculates the global stiffness matrix for this domain (steady-state case)
	void StiffnessMatrixSS(FESolver* psolver, bool bsymm);
	
	//! initialize class
	bool Initialize();

	//! activate
	void Activate();
	
	// update domain data
	void Update(const FETimeInfo& tp);

	// update element state data
	void UpdateElementStress(int iel, double dt);

public:

	// internal work (overridden from FEElasticDomain)
	void InternalForces(FEGlobalVector& R);

    // internal work (steady-state case)
    void InternalForcesSS(FEGlobalVector& R);
    
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
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf) {}
	void InertialForces(FEGlobalVector& R, vector<double>& F) {}
	void StiffnessMatrix(FESolver* psolver) {}
	void BodyForceStiffness(FESolver* psolver, FEBodyForce& bf) {}
	void MassMatrix(FESolver* psolver, double scale) {}

protected:
	FEMultiphasic*	m_pMat;
	int				m_dofP;		//!< pressure dof index
	int				m_dofC;		//!< concentration dof index
};
