#pragma once
#include "FECore/FEShellDomain.h"
#include "FEElasticDomain.h"
#include "FESolidMaterial.h"

//-----------------------------------------------------------------------------
//! Domain described by 3D shell elements
class FEElasticShellDomain : public FEShellDomain, public FEElasticDomain
{
public:
	FEElasticShellDomain(FEModel* pfem);

	//! \todo do I really need this?
	FEElasticShellDomain& operator = (FEElasticShellDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }

	//! Initialize domain
	bool Initialize();

	//! Activate the domain
	void Activate();

	//! update prior to solve
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

    //! evaluate a vector function over the shell
    vec3d evaluate(FEShellElement& el, vec3d* vn, vec3d* dvn, int n);
    
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
    
    //! calculates the solid element mass matrix
    void ElementMassMatrix(FEShellElement& el, matrix& ke, double a);
    
    //! calculates the stiffness matrix due to body forces
    void ElementBodyForceStiffness(FEBodyForce& bf, FEShellElement& el, matrix& ke);

public:
	//! Find interfaces between solid element faces and shell elements
	void FindSSI();

protected:
	FESolidMaterial*	m_pMat;
	bool                    m_binit;    //!< initialization flag
};
