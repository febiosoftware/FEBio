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
class FEElasticShellDomain : public FESSIShellDomain, public FEElasticDomain
{
public:
	FEElasticShellDomain(FEModel* pfem);

	//! \todo do I really need this?
	FEElasticShellDomain& operator = (FEElasticShellDomain& d);

	//! Activate the domain
	void Activate() override;

    //! initialize elements
    void PreSolveUpdate(const FETimeInfo& timeInfo) override;
    
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

	//! get the total dofs
	const FEDofList& GetDOFList() const override;

public: // overrides from FEElasticDomain

	//! calculates the residual
//	void Residual(FESolver* psolver, vector<double>& R);

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

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FELinearSystem& LS) override;

	// inertial stiffness
    void MassMatrix(FELinearSystem& LS, double scale) override;

	// body force stiffness
    void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) override;

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
    bool                m_update_dynamic;    //!< flag for updating quantities only used in dynamic analysis

	bool	m_secant_stress;	//!< use secant approximation to stress
	bool	m_secant_tangent;   //!< flag for using secant tangent

protected:
	FEDofList	m_dofV;
	FEDofList	m_dofSV;
	FEDofList	m_dofSA;
	FEDofList	m_dofR;
	FEDofList	m_dof;
};
