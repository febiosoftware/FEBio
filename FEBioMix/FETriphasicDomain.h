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
#include "FECore/FESolidDomain.h"
#include "FETriphasic.h"
#include "FEBioMech/FEElasticDomain.h"
#include <FECore/FEDofList.h>

//-----------------------------------------------------------------------------
//! Domain class for triphasic 3D solid elements
//! Note that this class inherits from FEElasticSolidDomain since this domain
//! also needs to calculate elastic stiffness contributions.
//!
class FEBIOMIX_API FETriphasicDomain : public FESolidDomain, public FEElasticDomain
{
public:
	//! constructor
	FETriphasicDomain(FEModel* pfem);
	
	//! reset domain data
	void Reset() override;

	//! get material (overridden from FEDomain)
	FEMaterial* GetMaterial() override;

	//! get the total dof
	const FEDofList& GetDOFList() const override;

	//! set the material
	void SetMaterial(FEMaterial* pmat) override;

	//! Unpack solid element data (overridden from FEDomain)
	void UnpackLM(FEElement& el, vector<int>& lm) override;

	//! activate
	void Activate() override;

	//! initialize elements for this domain
	void PreSolveUpdate(const FETimeInfo& timeInfo) override;
	
	// update domain data
	void Update(const FETimeInfo& tp) override;

	//! update element state data
	void UpdateElementStress(int iel);

public:
	// internal work (overridden from FEElasticDomain)
	void InternalForces(FEGlobalVector& R) override;

    // internal work (steady-state case)
    void InternalForcesSS(FEGlobalVector& R);
    
	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FELinearSystem& LS, bool bsymm);

	//! calculates the global stiffness matrix for this domain (steady-state case)
	void StiffnessMatrixSS(FELinearSystem& LS, bool bsymm);

protected:
	//! element internal force vector
	void ElementInternalForce(FESolidElement& el, vector<double>& fe);

    //! element internal force vector (steady-state case)
    void ElementInternalForceSS(FESolidElement& el, vector<double>& fe);
    
	//! calculates the element triphasic stiffness matrix
	bool ElementTriphasicStiffness(FESolidElement& el, matrix& ke, bool bsymm);
	
	//! calculates the element triphasic stiffness matrix
	bool ElementTriphasicStiffnessSS(FESolidElement& el, matrix& ke, bool bsymm);
	
protected: // overridden from FEElasticDomain, but not implemented in this domain
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override {}
	void InertialForces(FEGlobalVector& R, vector<double>& F) override {}
	void StiffnessMatrix(FELinearSystem& LS) override {}
	void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) override {}
	void MassMatrix(FELinearSystem& LS, double scale) override {}
	
protected:
	FETriphasic*	m_pMat;
	FEDofList		m_dofU;
	FEDofList		m_dofR;
	FEDofList		m_dof;
	int				m_dofP;		//!< pressure dof index
	int				m_dofC;		//!< concentration dof index
};
