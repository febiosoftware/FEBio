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
#include <FECore/FEDofList.h>
#include "FEBiphasicDomain.h"
#include "FEBiphasic.h"

//-----------------------------------------------------------------------------
//! Domain class for biphasic 3D solid elements
//! Note that this class inherits from FEElasticSolidDomain since the biphasic domain
//! also needs to calculate elastic stiffness contributions.
//!
class FEBIOMIX_API FEBiphasicSolidDomain : public FESolidDomain, public FEBiphasicDomain
{
public:
	//! constructor
	FEBiphasicSolidDomain(FEModel* pfem);

	//! initialize class
	bool Init() override;

	//! activate
	void Activate() override;

	//! reset domain data
	void Reset() override;

	//! intitialize element data
	void PreSolveUpdate(const FETimeInfo& timeInfo) override;

	//! Unpack solid element data  (overridden from FEDomain)
	void UnpackLM(FEElement& el, vector<int>& lm) override;

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial() override { return m_pMat; }

	//! set the material
	void SetMaterial(FEMaterial* pmat) override;

	//! get the total dof list
	const FEDofList& GetDOFList() const override;

	//! serialization
	void Serialize(DumpStream& ar) override;

public:
	// update domain data
	void Update(const FETimeInfo& tp) override;

	// update element stress
	void UpdateElementStress(int iel);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FELinearSystem& LS, bool bsymm) override;

	//! calculates the global stiffness matrix (steady-state case)
	void StiffnessMatrixSS(FELinearSystem& LS, bool bsymm) override;
	
public:
	// internal work (overridden from FEElasticDomain)
	void InternalForces(FEGlobalVector& R) override;

    // internal work (steady-state case)
    void InternalForcesSS(FEGlobalVector& R) override;
    
public:
	//! element internal force vector
	void ElementInternalForce(FESolidElement& el, vector<double>& fe);
	
    //! element internal force vector (steady-state case)
    void ElementInternalForceSS(FESolidElement& el, vector<double>& fe);
    
	//! calculates the element biphasic stiffness matrix
	bool ElementBiphasicStiffness(FESolidElement& el, matrix& ke, bool bsymm);
	
	//! calculates the element biphasic stiffness matrix for steady-state response
	bool ElementBiphasicStiffnessSS(FESolidElement& el, matrix& ke, bool bsymm);
	
public: // overridden from FEElasticDomain, but not all implemented in this domain
    void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override;
	void InertialForces(FEGlobalVector& R, vector<double>& F) override {}
	void StiffnessMatrix(FELinearSystem& LS) override {}
    void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) override;
    void ElementBodyForceStiffness(FEBodyForce& BF, FESolidElement &el, matrix &ke);
	void MassMatrix(FELinearSystem& LS, double scale) override {}

public: // biphasic domain "properties"
	// NOTE: I'm thinking about defining properties for domain classes. These would be similar to material
	// properties (and may require material properties to be evaluated), but are different in that they are
	// not meant to be customized. For example, the biphasic solver assumes Darcy's law in the evaluation 
	// of the fluid flux. Although one can see the fluid flux as a material property, since the code makes explicit
	// use of this constitutive equation (apparent from the fact that the biphasic material needs to define the permeability and its
	// strain derivate) it is not a true material property: i.e. it is not meant to be changed and is an inherent
	// assumption in this implementation. Consequently, the fluid flux would be a good example of a domain property.
	// That is why I've taken this calculation out of the FEBiphasic class and placed it here. 
	vec3d FluidFlux(FEMaterialPoint& mp) override;

	// Evaluate the nodal pressures
	// Note that the data vector stores the values for all of the nodes of the mesh, not just the domain nodes.
	// The values will be set to zero for nodes that don't belong to this domain.
	void GetNodalPressures(vector<double>& data);

private:
	// NOTE: This is a temporary construction. Just trying something out here.
	// Idea is here to construct a data export for the nodal pressures (the ones from the integration points, not the nodal dofs).
	vector<double>	m_nodePressure;	//!< nodal pressures projected from the integration points

	// This function updates the m_nodePressure variable
	void UpdateNodalPressures();

protected:
	int			m_varU, m_varP;	// displacement, pressure field indices

	FEDofList	m_dofU;		// displacement dofs
	FEDofList	m_dofSU;	// shell displacement dofs
	FEDofList	m_dofR;		// rigid rotation
	FEDofList	m_dof;
};
