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
#include <FEBioMech/FESSIShellDomain.h>
#include "FEBiphasicDomain.h"
#include "FEBiphasic.h"

//-----------------------------------------------------------------------------
//! Domain class for biphasic 3D shell elements
//!
class FEBIOMIX_API FEBiphasicShellDomain : public FESSIShellDomain, public FEBiphasicDomain
{
public:
    //! constructor
    FEBiphasicShellDomain(FEModel* pfem);
    
    //! activate
    void Activate() override;
    
    //! reset domain data
    void Reset() override;
    
    //! intitialize element data
    void PreSolveUpdate(const FETimeInfo& timeInfo) override;
    
    //! Unpack shell element data  (overridden from FEDomain)
    void UnpackLM(FEElement& el, vector<int>& lm) override;
    
	//! get material (overridden from FEDomain)
	FEMaterial* GetMaterial() override;

	//! get the total dof
	const FEDofList& GetDOFList() const override;

    //! set the material
    void SetMaterial(FEMaterial* pmat) override;
    
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
    void ElementInternalForce(FEShellElement& el, vector<double>& fe);
    
    //! element internal force vector (stead-state case)
    void ElementInternalForceSS(FEShellElement& el, vector<double>& fe);
    
    //! calculates the element biphasic stiffness matrix
    bool ElementBiphasicStiffness(FEShellElement& el, matrix& ke, bool bsymm);
    
    //! calculates the element biphasic stiffness matrix for steady-state response
    bool ElementBiphasicStiffnessSS(FEShellElement& el, matrix& ke, bool bsymm);

public: // overridden from FEElasticDomain, but not all implemented in this domain
    void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override;
    void ElementBodyForce(FEBodyForce& BF, FEShellElement& el, vector<double>& fe);
    void InertialForces(FEGlobalVector& R, vector<double>& F) override {}
    void StiffnessMatrix(FELinearSystem& LS) override {}
    void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) override;
    void ElementBodyForceStiffness(FEBodyForce& BF, FEShellElement &el, matrix &ke);
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

protected:
    int					m_dofSX;
    int					m_dofSY;
    int					m_dofSZ;
	FEDofList			m_dof;
};
