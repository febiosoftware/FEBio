/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FECore/FEDomain2D.h"
#include "FEFluidDomain.h"
#include "FEFluid.h"

//-----------------------------------------------------------------------------
//! domain described by 2D elements
//!
class FEBIOFLUID_API FEFluidDomain2D : public FEDomain2D, public FEFluidDomain
{
public:
    //! constructor
    FEFluidDomain2D(FEModel* pfem);
    ~FEFluidDomain2D() {}
    
    //! assignment operator
    FEFluidDomain2D& operator = (FEFluidDomain2D& d);
    
    //! initialize class
	bool Init() override;
    
    //! initialize elements
    void PreSolveUpdate(const FETimeInfo& timeInfo) override;
    
public: // overrides from FEDomain
    
    //! get the material
    FEMaterial* GetMaterial() override { return m_pMat; }
    
    //! set the material
    void SetMaterial(FEMaterial* pm) override;
    
public: // overrides from FEElasticDomain
    
    // update domain data
    void Update(const FETimeInfo& tp) override;
    
    // update the element stress
    void UpdateElementStress(int iel, const FETimeInfo& tp);
    
    //! internal stress forces
    void InternalForces(FEGlobalVector& R, const FETimeInfo& tp) override;
    
    //! body forces
    void BodyForce(FEGlobalVector& R, const FETimeInfo& tp, FEBodyForce& BF) override;
    
    //! intertial forces for dynamic problems
    void InertialForces(FEGlobalVector& R, const FETimeInfo& tp) override;
    
    //! calculates the global stiffness matrix for this domain
    void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) override;
    
    //! calculates inertial stiffness
    void MassMatrix(FESolver* psolver, const FETimeInfo& tp) override;
    
    //! body force stiffness
    void BodyForceStiffness(FESolver* psolver, const FETimeInfo& tp, FEBodyForce& bf) override;
    
public:
    // --- S T I F F N E S S ---
    
    //! calculates the solid element stiffness matrix
    void ElementStiffness(int iel, matrix& ke);
    
    //! material stiffness component
    void ElementMaterialStiffness(FEElement2D& el, matrix& ke);
    
    //! calculates the solid element mass matrix
    void ElementMassMatrix(FEElement2D& el, matrix& ke);
    
    //! calculates the stiffness matrix due to body forces
    void ElementBodyForceStiffness(FEBodyForce& bf, FEElement2D& el, matrix& ke);
    
    // --- R E S I D U A L ---
    
    //! Calculates the internal stress vector for solid elements
    void ElementInternalForce(FEElement2D& el, vector<double>& fe);
    
    //! Calculatess external body forces for solid elements
    void ElementBodyForce(FEBodyForce& BF, FEElement2D& elem, vector<double>& fe);
    
    //! Calculates the inertial force vector for solid elements
    void ElementInertialForce(FEElement2D& el, vector<double>& fe);
    
    // ---
    
protected:
    FEFluid*	m_pMat;

protected:
    int	m_dofWX, m_dofWY, m_dofWZ;
    int	m_dofEF;
    int	m_dofWXP, m_dofWYP, m_dofWZP;
    int m_dofEFP;
    int	m_dofAWX, m_dofAWY, m_dofAWZ;
    int m_dofAEF;
    int	m_dofAWXP, m_dofAWYP, m_dofAWZP;
    int m_dofAEFP;
};
