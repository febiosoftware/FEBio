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
#include <FECore/FEShellDomain.h>
#include <FECore/FEDofList.h>
#include "FEElasticDomain.h"
#include "FESolidMaterial.h"

//-----------------------------------------------------------------------------
//! Domain described by 3D shell elements
class FEElasticShellDomainOld : public FEShellDomainOld, public FEElasticDomain
{
public:
	FEElasticShellDomainOld(FEModel* pfem);

	//! \todo do I really need this?
	FEElasticShellDomainOld& operator = (FEElasticShellDomainOld& d);

	//! Initialize domain
	bool Init() override;

	//! Activate the domain
	void Activate() override;

	//! Unpack shell element data
	void UnpackLM(FEElement& el, vector<int>& lm) override;

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial() override { return m_pMat; }

	//! set the material
	void SetMaterial(FEMaterial* pmat) override;

	//! get total dof list
	const FEDofList& GetDOFList() const override;

public: // overrides from FEElasticDomain

	//! calculates the residual
//	void Residual(FESolver* psolver, vector<double>& R);

	//! internal stress forces
	void InternalForces(FEGlobalVector& R) override;

	//! Calculates inertial forces for dynamic problems | todo implement this (removed assert DSR)
	void InertialForces(FEGlobalVector& R, vector<double>& F) override { }

	//! calculate body force
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override;

	// update stresses
	void Update(const FETimeInfo& tp) override;

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FELinearSystem& LS) override;

	// inertial stiffness \todo implement this (removed assert DSR)
	void MassMatrix(FELinearSystem& LS, double scale) override { }

	// body force stiffness \todo implement this (removed assert DSR)
	void BodyForceStiffness  (FELinearSystem& LS, FEBodyForce& bf) override { }

public:
	//! calculates covariant basis vectors at an integration point
	void CoBaseVectors0(FEShellElementOld& el, int n, vec3d g[3]);

	//! calculates contravariant basis vectors at an integration point
	void ContraBaseVectors0(FEShellElementOld& el, int n, vec3d g[3]);

	// inverse jacobian with respect to reference frame
	double invjac0(FEShellElementOld& el, double J[3][3], int n);

	// jacobian with respect to reference frame
	double detJ0(FEShellElementOld& el, int n);

    //! calculates covariant basis vectors at an integration point
	void CoBaseVectors(FEShellElementOld& el, int n, vec3d g[3]);
    
    //! calculates contravariant basis vectors at an integration point
	void ContraBaseVectors(FEShellElementOld& el, int n, vec3d g[3]);
    
    // jacobian with respect to current frame
	double detJ(FEShellElementOld& el, int n);
    
	// calculate deformation gradient
	double defgrad(FEShellElementOld& el, mat3d& F, int n);

	// inverse jacobian with respect to current frame
	double invjact(FEShellElementOld& el, double J[3][3], int n);
 
public:

	// --- S T I F F N E S S --- 

	//! calculates the shell element stiffness matrix
	void ElementStiffness(int iel, matrix& ke);

	// --- R E S I D U A L ---

	//! Calculates the internal stress vector for shell elements
	void ElementInternalForce(FEShellElementOld& el, vector<double>& fe);

	//! Calculate extenral body forces for shell elements
	void ElementBodyForce(FEModel& fem, FEShellElementOld& el, vector<double>& fe);

	//! Calculate extenral body forces for shell elements
	void ElementBodyForce(FEBodyForce& BF, FEShellElementOld& el, vector<double>& fe);

protected:
	FESolidMaterial*	m_pMat;
	FEDofList			m_dofSU;	// shell displacement dofs
	FEDofList			m_dofSR;	// shell rotation dofs
	FEDofList			m_dofR;		// rigid rotation dofs
	FEDofList			m_dof;
};
