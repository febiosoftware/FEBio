//
//  FEFatigueMaterial.hpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 9/29/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#ifndef FEFatigueMaterial_hpp
#define FEFatigueMaterial_hpp

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
    FEElasticMaterial* GetElasticMaterial() override { return m_pBase; }
    
    // update fatigue material point at each iteration
    void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp) override;
    
public:
    
    //! Set the local coordinate system for a material point (overridden from FEMaterial)
    void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp) override;
    
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

#endif /* FEFatigueMaterial_hpp */
