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
#include "FESolutesMaterial.h"
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! domain described by 3D volumetric elements
//!
class FEBIOFLUID_API FESolutesDomain : public virtual FESolidDomain
{
public:
	//! constructor
	FESolutesDomain(FEModel* pfem);
	~FESolutesDomain() {}

	//! initialize elements
	void PreSolveUpdate(const FETimeInfo& timeInfo) override;

	void SetSteadyStateAnalysis() { m_btrans = false; }
	void SetTransientAnalysis() { m_btrans = true; }

public: // overrides from FEDomain

	//! get the material
	FEMaterial* GetMaterial() override { return m_pMat; }

	//! set the material
	void SetMaterial(FEMaterial* pm) override;

	//! get the total dofs
	const FEDofList& GetDOFList() const override;

public: // overrides from FEElasticDomain

		//! initialize class
	bool Init() override;

    //! serialize data to archive
    void Serialize(DumpStream& ar) override;

	//! Reset data
	void Reset() override;

	//! activate
	void Activate() override;

	//! initialize material points in the domain
	void InitMaterialPoints() override;

	// update stresses
	void Update(const FETimeInfo& tp) override;

	// update the element stress
	void UpdateElementStress(int iel, const FETimeInfo& tp);

	//! internal stress forces
	void InternalForces(FEGlobalVector& R);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FELinearSystem& LS);

public:
	// --- S T I F F N E S S ---

	//! calculates the solid element stiffness matrix
	void ElementStiffness(FESolidElement& el, matrix& ke);

	// --- R E S I D U A L ---

	//! Calculates the internal stress vector for solid elements
	void ElementInternalForce(FESolidElement& el, vector<double>& fe);

protected:
	bool	m_btrans;
	int     m_dofC;
	int     m_dofAC;
	FEDofList	m_dof;

	FESolutesMaterial*     m_pMat;
};
