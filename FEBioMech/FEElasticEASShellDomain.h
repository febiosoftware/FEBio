//
//  FEElasticEASShellDomain.hpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 12/6/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#ifndef FEElasticEASShellDomain_hpp
#define FEElasticEASShellDomain_hpp

#include "FESSIShellDomain.h"
#include "FEElasticDomain.h"
#include "FESolidMaterial.h"

//-----------------------------------------------------------------------------
//! Domain described by 3D shell elements
class FEElasticEASShellDomain : public FESSIShellDomain, public FEElasticDomain
{
public:
    FEElasticEASShellDomain(FEModel* pfem);
    
    //! \todo do I really need this?
    FEElasticEASShellDomain& operator = (FEElasticEASShellDomain& d);
    
    //! Initialize domain
    bool Initialize();
    
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
    void EvaluateEh(FEShellElement& el, const int n, const vec3d* Gcnt, mat3ds& E, vector<matrix>& hu, vector<matrix>& hw, vector<vec3d>& Nu, vector<vec3d>& Nw);
    
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
    
    // --- E A S  M E T H O D ---
    
    // Generate the G matrix for the EAS method
    void GenerateGMatrix(FEShellElement& el, const int n, const double Jeta, matrix& G);
    
    // Evaluate contravariant components of mat3ds tensor
    void mat3dsCntMat61(const mat3ds s, const vec3d* Gcnt, matrix& S);
    
    // Evaluate contravariant components of tens4ds tensor
    void tens4dsCntMat66(const tens4ds c, const vec3d* Gcnt, matrix& C);
    
    // Evaluate the matrices and vectors relevant to the EAS method
    void EvaluateEAS(FEShellElement& el, vector<double>& E, vector< vector<vec3d>>& HU, vector< vector<vec3d>>& HW, vector<mat3ds>& S, vector<tens4ds>& c);
    
    // Evaluate the strain using the ANS method
    void CollocationStrainsANS(FEShellElement& el, vector<double>& E, vector< vector<vec3d>>& HU, vector< vector<vec3d>>& HW, matrix& NS, matrix& NN);
    
    void EvaluateANS(FEShellElement& el, const int n, const vec3d* Gcnt, mat3ds& Ec, vector<matrix>& hu, vector<matrix>& hw, vector<double>& E, vector< vector<vec3d>>& HU, vector< vector<vec3d>>& HW);
    
    // Update alpha in EAS method
    void UpdateEAS(vector<double>& ui);
    void UpdateIncrementsEAS(vector<double>& ui, const bool binc);
    
protected:
    FESolidMaterial*    m_pMat;
    int                 m_nEAS;
};

#endif /* FEElasticEASShellDomain_hpp */
