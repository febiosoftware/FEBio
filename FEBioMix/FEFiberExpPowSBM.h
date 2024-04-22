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
#include <FEBioMech/FEElasticMaterial.h>
#include <FEBioMech/FERemodelingElasticMaterial.h>
#include <FECore/FEModelParam.h>
#include "febiomix_api.h"

//-----------------------------------------------------------------------------
//! Material class for single fiber, tension only
//! Exponential-power law
//! Fiber modulus depends on SBM content

class FEBIOMIX_API FEFiberExpPowSBM : public FEElasticMaterial, public FERemodelingInterface
{
public:
    FEFiberExpPowSBM(FEModel* pfem) : FEElasticMaterial(pfem) { m_sbm = -1; m_fiber = nullptr; }
    
    //! Initialization
    bool Init() override;
    
    //! Cauchy stress
    mat3ds Stress(FEMaterialPoint& mp) override;
    
    // Spatial tangent
    tens4ds Tangent(FEMaterialPoint& mp) override;
    
    //! Strain energy density
    double StrainEnergyDensity(FEMaterialPoint& mp) override;
    
    //! evaluate referential mass density
    double Density(FEMaterialPoint& pt) override;
    
    //! return fiber modulus
    double FiberModulus(FEMaterialPoint& pt, double rhor) { return m_ksi0(pt)*pow(rhor/m_rho0(pt), m_g(pt));}
    
    //! Create material point data
    FEMaterialPointData* CreateMaterialPointData() override;
    
    //! update specialize material point data
    void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp) override;
    
public: // --- remodeling interface ---

    //! calculate strain energy density at material point
    double StrainEnergy(FEMaterialPoint& pt) override;

    //! calculate tangent of strain energy density with mass density
    double Tangent_SE_Density(FEMaterialPoint& pt) override;
    
    //! calculate tangent of stress with mass density
    mat3ds Tangent_Stress_Density(FEMaterialPoint& pt) override;
    
public:
    FEParamDouble   m_alpha;	// coefficient of (In-1) in exponential
    FEParamDouble   m_beta;		// power of (In-1) in exponential
    FEParamDouble   m_ksi0;		// fiber modulus ksi = ksi0*(rhor/rho0)^gamma
    FEParamDouble   m_rho0;     // rho0
    FEParamDouble   m_g;        // gamma
    int             m_lsbm;     //!< local id of solid-bound molecule

public:
    FEVec3dValuator*    m_fiber;    //!< fiber orientation

    // declare the parameter list
    DECLARE_FECORE_CLASS();
};
