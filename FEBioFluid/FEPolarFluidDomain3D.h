/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2022 University of Utah, The Trustees of Columbia University in
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
#include "FEPolarFluidDomain.h"
#include "FEPolarFluid.h"

//-----------------------------------------------------------------------------
//! domain described by 3D volumetric elements
//!
class FEBIOFLUID_API FEPolarFluidDomain3D : public virtual FESolidDomain, public FEPolarFluidDomain
{
public:
    //! constructor
    FEPolarFluidDomain3D(FEModel* pfem);
    ~FEPolarFluidDomain3D() {}
    
    //! assignment operator
    FEPolarFluidDomain3D& operator = (FEPolarFluidDomain3D& d);
    
    //! serialize data to archive
    void Serialize(DumpStream& ar) override;

    //! initialize elements
    void PreSolveUpdate(const FETimeInfo& timeInfo) override;
    
public: // overrides from FEDomain
    
    //! get the material
    FEMaterial* GetMaterial() override { return m_pMat; }
    
    //! set the material
    void SetMaterial(FEMaterial* pm) override;
    
    // get total dof list
    const FEDofList& GetDOFList() const override;
    
public: // overrides from FEElasticDomain
    
    // update stresses
    void Update(const FETimeInfo& tp) override;
    
    // update the element stress
    void UpdateElementStress(int iel, const FETimeInfo& tp);
    
    //! internal stress forces
    void InternalForces(FEGlobalVector& R) override;
    
    //! body forces
    void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override;
    
    //! body moments
    void BodyMoment(FEGlobalVector& R, FEBodyMoment& bm) override;
    
    //! intertial forces for dynamic problems
    void InertialForces(FEGlobalVector& R) override;
    
    //! calculates the global stiffness matrix for this domain
    void StiffnessMatrix(FELinearSystem& LS) override;
    
    //! calculates inertial stiffness
    void MassMatrix(FELinearSystem& LS) override;
    
    //! body force stiffness
    void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) override;
    
    //! body momdent stiffness
    void BodyMomentStiffness(FELinearSystem& LS, FEBodyMoment& bm) override;
    
public:
    // --- S T I F F N E S S ---
    
    //! calculates the solid element stiffness matrix
    void ElementStiffness(FESolidElement& el, matrix& ke);
    
    //! calculates the solid element mass matrix
    void ElementMassMatrix(FESolidElement& el, matrix& ke);
    
    //! calculates the stiffness matrix due to body forces
    void ElementBodyForceStiffness(FEBodyForce& bf, FESolidElement& el, matrix& ke);
    
    //! calculates the stiffness matrix due to body moments
    void ElementBodyMomentStiffness(FEBodyMoment& bm, FESolidElement& el, matrix& ke);
    
    // --- R E S I D U A L ---
    
    //! Calculates the internal stress vector for solid elements
    void ElementInternalForce(FESolidElement& el, vector<double>& fe);
    
    //! Calculatess external body forces for solid elements
    void ElementBodyForce(FEBodyForce& BF, FESolidElement& elem, vector<double>& fe);
    
    //! Calculatess external body moments for solid elements
    void ElementBodyMoment(FEBodyMoment& bm, FESolidElement& elem, vector<double>& fe);
    
    //! Calculates the inertial force vector for solid elements
    void ElementInertialForce(FESolidElement& el, vector<double>& fe);
    
protected:
    FEPolarFluid*   m_pMat;
    
protected:
    FEDofList   m_dofW;
    FEDofList   m_dofAW;
    FEDofList   m_dofG;
    FEDofList   m_dofAG;
    FEDofList   m_dof;
    int         m_dofEF;
    int         m_dofAEF;
};
