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
#include <FECore/FEBeamDomain.h>
#include "FEElasticDomain.h"
#include "febiomech_api.h"

class FEElementMatrix;
class FEElasticBeamMaterial;

class FEBIOMECH_API FEElasticBeamDomain : public FEBeamDomain, public FEElasticDomain
{
public:
	FEElasticBeamDomain(FEModel* fem);

public: // from FEDomain
	// create function
	bool Create(int elements, FE_Element_Spec espec) override;

	//! Get the list of dofs on this domain
	const FEDofList& GetDOFList() const override;

	void SetMaterial(FEMaterial* pm) override;

public: // from FEMeshPartition

	// return number of beam elements
	int Elements() const override;

	//! return a reference to an element
	FEElement& ElementRef(int i) override;
	const FEElement& ElementRef(int i) const override;

	// update internal domain variables
	void Update(const FETimeInfo& tp) override;

public:
	//! return a beam element
	FEBeamElement& Element(int i);

public: // from FEElasticDomain
	//! calculate the internal forces
	void InternalForces(FEGlobalVector& R) override;

	//! Calculate the body force vector
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override;

	//! calculate the interial forces (for dynamic problems)
	void InertialForces(FEGlobalVector& R, std::vector<double>& F) override;

	//! Calculate global stiffness matrix (only contribution from internal force derivative)
	//! \todo maybe I should rename this the InternalStiffness matrix?
	void StiffnessMatrix(FELinearSystem& LS) override;

	//! Calculate stiffness contribution of body forces
	void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) override;

	//! calculate the mass matrix (for dynamic problems)
	void MassMatrix(FELinearSystem& LS, double scale) override;

private:
	void ElementInternalForces(FEBeamElement& el, std::vector<double>& fe);

	void ElementStiffnessMatrix(FEBeamElement& el, FEElementMatrix& ke);

	void UpdateElement(FEBeamElement& el);

private:
	FEDofList	m_dofs;
	FEElasticBeamMaterial* m_mat;
	std::vector<FEBeamElement>	m_Elem;
};
