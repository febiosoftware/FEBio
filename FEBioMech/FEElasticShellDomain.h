#pragma once
#include "FESSIShellDomain.h"
#include "FEElasticDomain.h"
#include "FESolidMaterial.h"

//-----------------------------------------------------------------------------
//! Domain described by 3D shell elements
class FEElasticShellDomain : public FESSIShellDomain, public FEElasticDomain
{
public:
	FEElasticShellDomain(FEModel* pfem);

	//! \todo do I really need this?
	FEElasticShellDomain& operator = (FEElasticShellDomain& d);

	//! Initialize domain
	bool Init() override;

	//! Activate the domain
	void Activate();

    //! initialize elements
    void PreSolveUpdate(const FETimeInfo& timeInfo);
    
	//! Unpack shell element data
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial() { return m_pMat; }

	//! set the material
	void SetMaterial(FEMaterial* pmat);

public: // overrides from FEElasticDomain

	//! calculates the residual
//	void Residual(FESolver* psolver, vector<double>& R);

	//! internal stress forces
	void InternalForces(FEGlobalVector& R);

	//! Calculates inertial forces for dynamic problems
    void InertialForces(FEGlobalVector& R, vector<double>& F);

	//! calculate body force
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf);

	// update stresses
	void Update(const FETimeInfo& tp);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolver* psolver);

	// inertial stiffness
    void MassMatrix(FESolver* psolver, double scale);

	// body force stiffness
    void BodyForceStiffness  (FESolver* psolver, FEBodyForce& bf);

public:

	// --- S T I F F N E S S --- 

	//! calculates the shell element stiffness matrix
	void ElementStiffness(int iel, matrix& ke);

    //! calculates the solid element mass matrix
    void ElementMassMatrix(FEShellElement& el, matrix& ke, double a);
    
    //! calculates the stiffness matrix due to body forces
    void ElementBodyForceStiffness(FEBodyForce& bf, FEShellElement& el, matrix& ke);
    
	// --- R E S I D U A L ---

	//! Calculates the internal stress vector for shell elements
	void ElementInternalForce(FEShellElement& el, vector<double>& fe);

	//! Calculate extenral body forces for shell elements
	void ElementBodyForce(FEModel& fem, FEShellElement& el, vector<double>& fe);

	//! Calculate extenral body forces for shell elements
	void ElementBodyForce(FEBodyForce& BF, FEShellElement& el, vector<double>& fe);
    
    //! Calculates the inertial force for shell elements
    void ElementInertialForce(FEShellElement& el, vector<double>& fe);
    
protected:
	FESolidMaterial*	m_pMat;
    double              m_alphaf;
    double              m_alpham;
    double              m_beta;
};
