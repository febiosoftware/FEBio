//
//  FEDamageCriterion.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 9/18/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEDamageCriterion.h"
#include "FEDamageMaterial.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
// Simo's damage criterion uses sqrt(2*strain energy density)
double FEDamageCriterionSimo::DamageCriterion(FEMaterialPoint& pt)
{
    FEElasticMaterial* m_pMat = dynamic_cast<FEElasticMaterial*>(GetParent());
    FEElasticMaterial* m_pBase = m_pMat;
    double sed = m_pBase->StrainEnergyDensity(pt);
    
    // clean up round-off errors
    if (sed < 0) sed = 0;
    
    return sqrt(2*sed);
}

//-----------------------------------------------------------------------------
// Strain energy density damage criterion
double FEDamageCriterionSED::DamageCriterion(FEMaterialPoint& pt)
{
    FEElasticMaterial* m_pMat = dynamic_cast<FEElasticMaterial*>(GetParent());
    FEElasticMaterial* m_pBase = m_pMat;
    double sed = m_pBase->StrainEnergyDensity(pt);
    
    return sed;
}

//-----------------------------------------------------------------------------
// Specific strain energy damage criterion
double FEDamageCriterionSSE::DamageCriterion(FEMaterialPoint& pt)
{
    FEElasticMaterial* m_pMat = dynamic_cast<FEElasticMaterial*>(GetParent());
    FEElasticMaterial* m_pBase = m_pMat;
    double sed = m_pBase->StrainEnergyDensity(pt);
    FEParamDouble& rho = m_pBase->Density();
    
    return sed/rho(pt);
}

//-----------------------------------------------------------------------------
// von Mises stress damage criterion
double FEDamageCriterionVMS::DamageCriterion(FEMaterialPoint& pt)
{
    FEElasticMaterial* m_pMat = dynamic_cast<FEElasticMaterial*>(GetParent());
    FEElasticMaterial* m_pBase = m_pMat;
    mat3ds s = m_pBase->Stress(pt);
    
    double vms = sqrt((SQR(s.xx()-s.yy()) + SQR(s.yy()-s.zz()) + SQR(s.zz()-s.xx())
                       + 6*(SQR(s.xy()) + SQR(s.yz()) + SQR(s.xz())))/2);
    
    return vms;
}

//-----------------------------------------------------------------------------
// max shear stress damage criterion
double FEDamageCriterionMSS::DamageCriterion(FEMaterialPoint& pt)
{
    // evaluate stress tensor
    FEElasticMaterial* m_pMat = dynamic_cast<FEElasticMaterial*>(GetParent());
    FEElasticMaterial* m_pBase = m_pMat;
    mat3ds s = m_pBase->Stress(pt);
    
    // evaluate principal normal stresses
    double ps[3], ms[3];
    s.eigen2(ps);
    ms[0] = fabs(ps[1] - ps[2])/2;
    ms[1] = fabs(ps[2] - ps[0])/2;
    ms[2] = fabs(ps[0] - ps[1])/2;
    
    // evaluate max shear stresses and use largest of three
    double mss = max(max(ms[0],ms[1]),ms[2]);
    
    return mss;
}

//-----------------------------------------------------------------------------
// max normal stress damage criterion
double FEDamageCriterionMNS::DamageCriterion(FEMaterialPoint& pt)
{
    // evaluate stress tensor
    FEElasticMaterial* m_pMat = dynamic_cast<FEElasticMaterial*>(GetParent());
    FEElasticMaterial* m_pBase = m_pMat;
    mat3ds s = m_pBase->Stress(pt);
    
    // evaluate principal normal stresses
    double ps[3];
    s.eigen2(ps);
    
    // evaluate max normal stress
    double mns = max(max(ps[0],ps[1]),ps[2]);
    
    return mns;
}

//-----------------------------------------------------------------------------
// max normal Lagrange strain damage criterion
double FEDamageCriterionMNLS::DamageCriterion(FEMaterialPoint& mp)
{
    // evaluate strain tensor
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    mat3ds E = pt.Strain();
    
    // evaluate principal normal strains
    double ps[3];
    E.eigen2(ps);
    
    // evaluate max normal Lagrange strain
    double mnls = max(max(ps[0],ps[1]),ps[2]);
    
    return mnls;
}
