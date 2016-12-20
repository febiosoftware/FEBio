#pragma once
#include "FECore/FEShellDomain.h"
#include "FEElasticDomain.h"
#include "FESolidMaterial.h"

//-----------------------------------------------------------------------------
//! Domain described by 3D shell elements
class FEElasticShellDomainOld : public FEShellDomain, public FEElasticDomain
{
public:
	FEElasticShellDomainOld(FEModel* pfem);

	//! \todo do I really need this?
	FEElasticShellDomainOld& operator = (FEElasticShellDomainOld& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }

	//! Initialize domain
	bool Initialize(FEModel& fem);

	//! Activate the domain
	void Activate();

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

	//! Calculates inertial forces for dynamic problems | todo implement this (removed assert DSR)
	void InertialForces(FEGlobalVector& R, vector<double>& F) { }

	//! calculate body force
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf);

	// update stresses
	void Update(const FETimeInfo& tp);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolver* psolver);

	// inertial stiffness \todo implement this (removed assert DSR)
	void MassMatrix(FESolver* psolver, double scale) { }

	// body force stiffness \todo implement this (removed assert DSR)
	void BodyForceStiffness  (FESolver* psolver, FEBodyForce& bf) { }

public:
	//! calculates covariant basis vectors at an integration point
	void CoBaseVectors0(FEShellElement& el, int n, vec3d g[3]);

	//! calculates contravariant basis vectors at an integration point
	void ContraBaseVectors0(FEShellElement& el, int n, vec3d g[3]);

	// inverse jacobian with respect to reference frame
	double invjac0(FEShellElement& el, double J[3][3], int n);

	// jacobian with respect to reference frame
	double detJ0(FEShellElement& el, int n);

    //! calculates covariant basis vectors at an integration point
    void CoBaseVectors(FEShellElement& el, int n, vec3d g[3]);
    
    //! calculates contravariant basis vectors at an integration point
    void ContraBaseVectors(FEShellElement& el, int n, vec3d g[3]);
    
    // jacobian with respect to current frame
    double detJ(FEShellElement& el, int n);
    
	// calculate deformation gradient
	double defgrad(FEShellElement& el, mat3d& F, int n);

	// inverse jacobian with respect to current frame
	double invjact(FEShellElement& el, double J[3][3], int n);
 
public:

	// --- S T I F F N E S S --- 

	//! calculates the shell element stiffness matrix
	void ElementStiffness(int iel, matrix& ke);

	// --- R E S I D U A L ---

	//! Calculates the internal stress vector for shell elements
	void ElementInternalForce(FEShellElement& el, vector<double>& fe);

	//! Calculate extenral body forces for shell elements
	void ElementBodyForce(FEModel& fem, FEShellElement& el, vector<double>& fe);

	//! Calculate extenral body forces for shell elements
	void ElementBodyForce(FEBodyForce& BF, FEShellElement& el, vector<double>& fe);

protected:
	FESolidMaterial*	m_pMat;
	int					m_dofU;
	int					m_dofV;
	int					m_dofW;
};
