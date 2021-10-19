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
#include "FEDiscreteElasticMaterial.h"
#include <FECore/FEDiscreteDomain.h>
#include "FEElasticDomain.h"

class FEDiscreteElasticDomain : public FEDiscreteDomain, public FEElasticDomain
{
public:
	FEDiscreteElasticDomain(FEModel* fem);

	//! Unpack LM data
	void UnpackLM(FEElement& el, vector<int>& lm) override;

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial() override;

	//! set the material
	void SetMaterial(FEMaterial* pmat) override;

	void Activate() override;

	// get the total dofs
	const FEDofList& GetDOFList() const override;

	void PreSolveUpdate(const FETimeInfo& timeInfo) override;

public: // overridden from FEElasticDomain

		//! calculate stiffness matrix
	void StiffnessMatrix(FELinearSystem& LS) override;
	void MassMatrix(FELinearSystem& LS, double scale) override {}
	void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) override {}

	//! Calculates inertial forces for dynamic problems | todo implement (removed assert DSR)
	void InertialForces(FEGlobalVector& R, vector<double>& F) override { }

	//! update domain data
	void Update(const FETimeInfo& tp) override;

	//! internal stress forces
	void InternalForces(FEGlobalVector& R) override;

	//! calculate bodyforces (not used since springs are considered mass-less)
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override {}

protected:
	FEDiscreteElasticMaterial*	m_pMat;
	FEDofList	m_dofU, m_dofR, m_dof;
};
