//
//  FEElasticANSShellDomain.hpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 12/20/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#ifndef FEElasticANSShellDomain_hpp
#define FEElasticANSShellDomain_hpp

#include "FESSIShellDomain.h"
#include "FEElasticDomain.h"
#include "FESolidMaterial.h"

//-----------------------------------------------------------------------------
//! Domain described by 3D shell elements
class FEElasticANSShellDomain : public FESSIShellDomain, public FEElasticDomain
{
public:
    FEElasticANSShellDomain(FEModel* pfem);
    
    //! \todo do I really need this?
    FEElasticANSShellDomain& operator = (FEElasticANSShellDomain& d);
    
    //! Initialize domain
	bool Init() override;
    
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
    //    void Residual(FESolver* psolver, vector<double>& R);
    
    //! internal stress forces
    void InternalForces(FEGlobalVector& R);
    
    //! Calculates inertial forces for dynamic problems
    void InertialForces(FEGlobalVector& R, vector<double>& F);
    
    //! calculate body force
    void BodyForce(FEGlobalVector& R, FEBodyForce& bf);
    
    // update stresses
    void Update(const FETimeInfo& tp);
    
    //! initialize elements for this domain
    void PreSolveUpdate(const FETimeInfo& timeInfo);
    
    //! calculates the global stiffness matrix for this domain
    void StiffnessMatrix(FESolver* psolver);
    
    // inertial stiffness
    void MassMatrix(FESolver* psolver, double scale);
    
    // body force stiffness
    void BodyForceStiffness  (FESolver* psolver, FEBodyForce& bf);
    
    // evaluate strain E and matrix hu and hw
	void EvaluateEh(FEShellElementNew& el, const int n, const vec3d* Gcnt, mat3ds& E, vector<matrix>& hu, vector<matrix>& hw, vector<vec3d>& Nu, vector<vec3d>& Nw);
    
public:
    
    // --- S T I F F N E S S ---
    
    //! calculates the shell element stiffness matrix
    void ElementStiffness(int iel, matrix& ke);
    
    // --- R E S I D U A L ---
    
    //! Calculates the internal stress vector for shell elements
	void ElementInternalForce(FEShellElementNew& el, vector<double>& fe);
    
    //! Calculate extenral body forces for shell elements
	void ElementBodyForce(FEModel& fem, FEShellElementNew& el, vector<double>& fe);
    
    //! Calculate extenral body forces for shell elements
	void ElementBodyForce(FEBodyForce& BF, FEShellElementNew& el, vector<double>& fe);
    
    //! calculates the solid element mass matrix
	void ElementMassMatrix(FEShellElementNew& el, matrix& ke, double a);
    
    //! calculates the stiffness matrix due to body forces
	void ElementBodyForceStiffness(FEBodyForce& bf, FEShellElementNew& el, matrix& ke);
    
public:
    
    // --- A N S  M E T H O D ---
    
    // Evaluate contravariant components of mat3ds tensor
    void mat3dsCntMat61(const mat3ds s, const vec3d* Gcnt, matrix& S);
    
    // Evaluate contravariant components of tens4ds tensor
    void tens4dsCntMat66(const tens4ds c, const vec3d* Gcnt, matrix& C);
    
    // Evaluate the strain using the ANS method
	void CollocationStrainsANS(FEShellElementNew& el, vector<double>& E, vector< vector<vec3d>>& HU, vector< vector<vec3d>>& HW, matrix& NS, matrix& NN);
    
	void EvaluateANS(FEShellElementNew& el, const int n, const vec3d* Gcnt, mat3ds& Ec, vector<matrix>& hu, vector<matrix>& hw, vector<double>& E, vector< vector<vec3d>>& HU, vector< vector<vec3d>>& HW);
    
protected:
    FESolidMaterial*    m_pMat;
};

#endif /* FEElasticANSShellDomain_hpp */
