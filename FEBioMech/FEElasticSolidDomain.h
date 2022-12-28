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
#include <FECore/FESolidDomain.h>
#include "FEElasticDomain.h"
#include "FESolidMaterial.h"
#include <FECore/FEDofList.h>

//-----------------------------------------------------------------------------
//! domain described by Lagrange-type 3D volumetric elements
//!
class FEBIOMECH_API FEElasticSolidDomain : public FESolidDomain, public FEElasticDomain
{
public:
	//! constructor
	FEElasticSolidDomain(FEModel* pfem);

	//! assignment operator
	FEElasticSolidDomain& operator = (FEElasticSolidDomain& d);

	//! activate
	void Activate() override;

	//! initialize elements
	void PreSolveUpdate(const FETimeInfo& timeInfo) override;

	//! Unpack solid element data
	void UnpackLM(FEElement& el, vector<int>& lm) override;

	//! Set flag for update for dynamic quantities
	void SetDynamicUpdateFlag(bool b);

	//! serialization
	void Serialize(DumpStream& ar) override;

	// get the total dof list
	const FEDofList& GetDOFList() const override;

public: // overrides from FEDomain

	//! get the material
	FEMaterial* GetMaterial() override { return m_pMat; }

	//! set the material
	void SetMaterial(FEMaterial* pm) override;

public: // overrides from FEElasticDomain

	// update stresses
	void Update(const FETimeInfo& tp) override;

	// update the element stress
	virtual void UpdateElementStress(int iel, const FETimeInfo& tp);

	//! intertial forces for dynamic problems
	void InertialForces(FEGlobalVector& R, vector<double>& F) override;

	//! internal stress forces
	void InternalForces(FEGlobalVector& R) override;

	//! body forces
	void BodyForce(FEGlobalVector& R, FEBodyForce& BF) override;

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FELinearSystem& LS) override;

	//! calculates inertial stiffness
	void MassMatrix(FELinearSystem& LS, double scale) override;

	//! body force stiffness
	void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) override;

public:
	// --- S T I F F N E S S ---

	//! calculates the solid element stiffness matrix
	virtual void ElementStiffness(const FETimeInfo& tp, int iel, matrix& ke);

	//! geometrical stiffness (i.e. initial stress)
	virtual void ElementGeometricalStiffness(FESolidElement& el, matrix& ke);

	//! material stiffness component
	virtual void ElementMaterialStiffness(FESolidElement& el, matrix& ke);

	// --- R E S I D U A L ---

	//! Calculates the internal stress vector for solid elements
	void ElementInternalForce(FESolidElement& el, vector<double>& fe);

    //! Calculates the inertial force vector for solid elements
    void ElementInertialForce(FESolidElement& el, vector<double>& fe);
    
protected:
    double              m_alphaf;
    double              m_alpham;
    double              m_beta;
	bool				m_update_dynamic;	//!< flag for updating quantities only used in dynamic analysis

	bool	m_secant_stress;	//!< use secant approximation to stress
	bool	m_secant_tangent;   //!< flag for using secant tangent

protected:
	FEDofList	m_dofU;		// displacement dofs
	FEDofList	m_dofR;		// rigid rotation rofs
	FEDofList	m_dofSU;	// shell displacement dofs
	FEDofList	m_dofV;		// velocity dofs
	FEDofList	m_dofSV;	// shell velocity dofs
	FEDofList	m_dofSA;	// shell acceleration dofs
	FEDofList	m_dof;		// total dof list

	FESolidMaterial*	m_pMat;
};

class FEStandardElasticSolidDomain : public FEElasticSolidDomain
{
public:
	FEStandardElasticSolidDomain(FEModel* fem);

private:
	std::string		m_elemType;

	DECLARE_FECORE_CLASS();
};
