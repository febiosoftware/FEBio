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
#include "FEElasticMaterial.h"
#include "FEDamageCriterion.h"
#include "FEDamageCDF.h"
#include "FEReactiveFatigueMaterialPoint.h"
#include <FECore/FEMaterial.h>

//-----------------------------------------------------------------------------
// This material models fatigue and damage in any hyper-elastic materials.

class FEReactiveFatigue : public FEElasticMaterial
{
public:
    FEReactiveFatigue(FEModel* pfem);
    
public:
    //! calculate stress at material point
    mat3ds Stress(FEMaterialPoint& pt) override;
    
    //! calculate tangent stiffness at material point
    tens4ds Tangent(FEMaterialPoint& pt) override;
    
    //! calculate strain energy density at material point
    double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
    //! damage
    double Damage(FEMaterialPoint& pt);
    
    //! data initialization and checking
    bool Init() override;
    
    // returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override;
    
    // get the elastic material
    FEElasticMaterial* GetElasticMaterial() override { return m_pBase; }
    
    // update fatigue material point at each iteration
    void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp) override;
    
public:
    
    FEElasticMaterial*  m_pBase;    // base elastic material
    FEDamageCDF*        m_pIdmg;    // damage model for intact bonds
    FEDamageCDF*        m_pFdmg;    // damage model for fatigued bonds
    FEDamageCriterion*  m_pIcrt;    // damage criterion
    FEDamageCriterion*  m_pFcrt;    // fatigue criterion
    
public:
    FEParamDouble       m_k0;       // reaction rate for fatigue reaction
    FEParamDouble       m_beta;     // power exponent for fatigue reaction
    
    DECLARE_FECORE_CLASS();
};
