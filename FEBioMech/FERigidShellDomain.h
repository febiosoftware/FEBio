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
#include "FEElasticShellDomain.h"
#include "FEElasticShellDomainOld.h"

class FEBodyForce;
class FERigidMaterial;

//-----------------------------------------------------------------------------
//! domain class for 3D rigid shell elements
//!
class FERigidShellDomainOld : public FEElasticShellDomainOld
{
public:
	//! constructor
	FERigidShellDomainOld(FEModel* pfem);

	//! Initialize
	bool Init() override;

	//! reset data
	void Reset() override;

public:
	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FELinearSystem& LS) override;

	//! calculates the internal forces (nothing to do)
	void InternalForces(FEGlobalVector& R) override;

	//! calculates mass matrix (nothing to do)
	void MassMatrix(FELinearSystem& LS, double scale) override;

	//! calculates the inertial forces (nothing to do)
	void InertialForces(FEGlobalVector& R, std::vector<double>& F) override;

	// update domain data
	void Update(const FETimeInfo& tp) override;
};

//-----------------------------------------------------------------------------
//! domain class for 3D rigid shell elements
//!
class FERigidShellDomain : public FEElasticShellDomain
{
public:
	//! constructor
	FERigidShellDomain(FEModel* pfem);

	//! Initialize
	bool Init() override;

	//! reset data
	void Reset() override;

	//! initialize shells 
//	void InitShells() override;

public:
	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FELinearSystem& LS) override;

	//! calculates the internal forces (nothing to do)
	void InternalForces(FEGlobalVector& R) override;

	//! calculates mass matrix (nothing to do)
	void MassMatrix(FELinearSystem& LS, double scale) override;

	//! calculates the inertial forces (nothing to do)
	void InertialForces(FEGlobalVector& R, std::vector<double>& F) override;

	// update domain data
	void Update(const FETimeInfo& tp) override;
};

//-----------------------------------------------------------------------------
// Implements a rigid shell domain. We need to inherit from FEShellDomain and FEElasticDomain.
// The latter is needed because most of the solid solver classes assume that domains inherit this class.
class FERigidShellDomainNew : public FEShellDomain, public FEElasticDomain
{
public:
	FERigidShellDomainNew(FEModel* fem);

public: // from FEMeshPartition
	//! return number of elements
	int Elements() const override { return (int)m_Elem.size(); };

	//! return a reference to an element \todo this is not the preferred interface but I've added it for now
	FEElement& ElementRef(int i) override { return m_Elem[i]; };
	const FEElement& ElementRef(int i) const override { return m_Elem[i]; };

public: // from FEDomain
	// create function
	bool Create(int elements, FE_Element_Spec espec) override;

	const FEDofList& GetDOFList() const override { return m_dof; }

	void Update(const FETimeInfo& tp) override;

	void Reset() override;

	void SetMaterial(FEMaterial* pm) override;

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial() override;

public: // from FEShellDomain
	// get a shell element
	FEShellElement& Element(int i) override { return m_Elem[i]; }

	void AssignDefaultShellThickness() override;

public: // from FEElasticDomain
	void InternalForces(FEGlobalVector& R) override;

	void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override;

	void InertialForces(FEGlobalVector& R, std::vector<double>& F) override;

	void StiffnessMatrix(FELinearSystem& LS) override;

	void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) override;

	void MassMatrix(FELinearSystem& LS, double scale) override;

private:
	//! Calculate extenral body forces for shell elements
	void ElementBodyForce(FEBodyForce& BF, FEShellElement& el, vector<double>& fe);

	void ElementBodyForceStiffness(FEBodyForce& BF, FEShellElement& el, matrix& ke);

public:
	double detJ0(FEShellElement& el, int n);

protected:
	double	m_h0; // TODO: move to base class?

protected:
	std::vector<FEShellElement>	m_Elem;
	FEDofList	m_dof;
	FERigidMaterial* m_pMat;

	DECLARE_FECORE_CLASS();
};
