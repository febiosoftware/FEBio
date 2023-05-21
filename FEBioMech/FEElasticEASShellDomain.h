/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include "FESSIShellDomain.h"
#include "FEElasticDomain.h"
#include "FESolidMaterial.h"

//-----------------------------------------------------------------------------
//! Domain described by 3D shell elements
class FEBIOMECH_API FEElasticEASShellDomain : public FESSIShellDomain, public FEElasticDomain
{
public:
    FEElasticEASShellDomain(FEModel* pfem);
    
    //! \todo do I really need this?
    FEElasticEASShellDomain& operator = (FEElasticEASShellDomain& d);
    
    //! Initialize domain
	bool Init() override;
    
    //! Activate the domain
    void Activate() override;
    
    //! Unpack shell element data
    void UnpackLM(FEElement& el, vector<int>& lm) override;
    
    //! Set flag for update for dynamic quantities
    void SetDynamicUpdateFlag(bool b);
    
    //! serialization
    void Serialize(DumpStream& ar) override;
    
    //! get the material (overridden from FEDomain)
    FEMaterial* GetMaterial() override { return m_pMat; }
    
    //! set the material
    void SetMaterial(FEMaterial* pmat) override;

	// get the total dof list
	const FEDofList& GetDOFList() const override;
    
public: // overrides from FEElasticDomain
    
    //! calculates the residual
    //    void Residual(FESolver* psolver, vector<double>& R);
    
    //! internal stress forces
    void InternalForces(FEGlobalVector& R) override;
    
    //! Calculates inertial forces for dynamic problems
    void InertialForces(FEGlobalVector& R, vector<double>& F) override;
    
    //! calculate body force
    void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override;
    
    // update stresses
    void Update(const FETimeInfo& tp) override;
    
    // update the element stress
    void UpdateElementStress(int iel, const FETimeInfo& tp);
    
    //! initialize elements for this domain
    void PreSolveUpdate(const FETimeInfo& timeInfo) override;
    
    //! calculates the global stiffness matrix for this domain
    void StiffnessMatrix(FELinearSystem& LS) override;
    
    // inertial stiffness
    void MassMatrix(FELinearSystem& LS, double scale) override;
    
    // body force stiffness
    void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) override;
    
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
    
    //! Calculates the inertial force for shell elements
    void ElementInertialForce(FEShellElementNew& el, vector<double>& fe);

    //! calculates the solid element mass matrix
	void ElementMassMatrix(FEShellElementNew& el, matrix& ke, double a);
    
    //! calculates the stiffness matrix due to body forces
	void ElementBodyForceStiffness(FEBodyForce& bf, FEShellElementNew& el, matrix& ke);
    
public:
    
    // --- E A S  M E T H O D ---
    
    // Generate the G matrix for the EAS method
	void GenerateGMatrix(FEShellElementNew& el, const int n, const double Jeta, matrix& G);
    
    // Evaluate contravariant components of mat3ds tensor
    void mat3dsCntMat61(const mat3ds s, const vec3d* Gcnt, matrix& S);
    
    // Evaluate contravariant components of tens4ds tensor
    void tens4dsCntMat66(const tens4ds c, const vec3d* Gcnt, matrix& C);
    void tens4dmmCntMat66(const tens4dmm c, const vec3d* Gcnt, matrix& C);

    // Evaluate the matrices and vectors relevant to the EAS method
	void EvaluateEAS(FEShellElementNew& el, vector<double>& E, vector< vector<vec3d>>& HU, vector< vector<vec3d>>& HW, vector<mat3ds>& S, vector<tens4dmm>& c);
    
    // Evaluate the strain using the ANS method
	void CollocationStrainsANS(FEShellElementNew& el, vector<double>& E, vector< vector<vec3d>>& HU, vector< vector<vec3d>>& HW, matrix& NS, matrix& NN);
    
	void EvaluateANS(FEShellElementNew& el, const int n, const vec3d* Gcnt, mat3ds& Ec, vector<matrix>& hu, vector<matrix>& hw, vector<double>& E, vector< vector<vec3d>>& HU, vector< vector<vec3d>>& HW);
    
    // Update alpha in EAS method
    void UpdateEAS(vector<double>& ui) override;
    void UpdateIncrementsEAS(vector<double>& ui, const bool binc) override;
    
protected:
    FESolidMaterial*    m_pMat;
    int                 m_nEAS;
    bool                m_update_dynamic;    //!< flag for updating quantities only used in dynamic analysis
    
    bool    m_secant_stress;    //!< use secant approximation to stress
    bool    m_secant_tangent;   //!< flag for using secant tangent

    FEDofList   m_dofV;
    FEDofList   m_dofSV;
	FEDofList	m_dofSA;
	FEDofList	m_dofR;
	FEDofList	m_dof;
};
