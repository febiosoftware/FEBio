/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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

//-----------------------------------------------------------------------------
// Define a material point that stores the fatigue and damage variables.
class FEFatigueMaterialPoint : public FEMaterialPoint
{
public:
    FEFatigueMaterialPoint(FEMaterialPoint *pt) : FEMaterialPoint(pt) {}
    
    FEMaterialPoint* Copy();
    
    void Init();
    void Update(const FETimeInfo& timeInfo);
    
    void Serialize(DumpStream& ar);
    
public:
    double      m_D;            //!< damage (0 = no damage, 1 = complete damage)
    
    double      m_wi;           //!< intact bond fraction at intermediate time
    double      m_wip;          //!< intact bond fraction at previous time
    double      m_wit;          //!< intact bond fraction at current time

    double      m_awi;          //!< time deriv of intact bond fraction at intermediate time
    double      m_awip;         //!< time deriv of intact bond fraction at previous time
    double      m_awit;         //!< time deriv of intact bond fraction at current time

    double      m_wf;           //!< fatigued bond fraction at intermediate time
    double      m_wfp;          //!< fatigued bond fraction at previous time
    double      m_wft;          //!< fatigued bond fraction at current time

    double      m_awf;          //!< time deriv of fatigued bond fraction at intermediate time
    double      m_awfp;         //!< time deriv of fatigued bond fraction at previous time
    double      m_awft;         //!< time deriv of fatigued bond fraction at current time
    
    double      m_Xd;           //!< damage criterion at intermediate time
    double      m_Xdp;          //!< damage criterion at previous time
    double      m_Xdt;          //!< damage criterion at current time
    
    double      m_aXd;          //!< time deriv of damage criterion at intermediate time
    double      m_aXdp;         //!< time deriv of damage criterion at previous time
    double      m_aXdt;         //!< time deriv of damage criterion at current time
    
    double      m_Xf;           //!< fatigue criterion at intermediate time
    double      m_Xfp;          //!< fatigue criterion at previous time
    double      m_Xft;          //!< fatigue criterion at current time
    
    double      m_aXf;          //!< time deriv of fatigue criterion at intermediate time
    double      m_aXfp;         //!< time deriv of fatigue criterion at previous time
    double      m_aXft;         //!< time deriv of fatigue criterion at current time
    
};

//-----------------------------------------------------------------------------
// This material models fatigue and damage in any hyper-elastic materials.

class FEFatigueMaterial : public FEElasticMaterial
{
public:
    FEFatigueMaterial(FEModel* pfem);
    
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
    FEMaterialPoint* CreateMaterialPointData() override
    {
        return new FEFatigueMaterialPoint(m_pBase->CreateMaterialPointData());
    }
    
    // get the elastic material
    FEElasticMaterial* GetElasticMaterial() { return m_pBase; }
    
    // update fatigue material point at each iteration
    void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp) override;
    
public:
    
    FEElasticMaterial*  m_pBase;    // base elastic material
    FEDamageCDF*        m_pIdmg;    // damage model for intact bonds
    FEDamageCDF*        m_pFdmg;    // damage model for fatigued bonds
    FEDamageCriterion*  m_pDcrt;    // damage criterion
    FEDamageCriterion*  m_pFcrt;    // fatigue criterion
    
public:
    double      m_k0;       // reaction rate for fatigue reaction
    double      m_beta;     // power exponent for fatigue reaction
    
    DECLARE_FECORE_CLASS();
};
