//
//  FEDamageCriterionUC.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 9/19/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEDamageCriterionUC.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
// Simo's damage criterion uses sqrt(2*strain energy density)
double FEDamageCriterionUCSimo::DamageCriterion(FEMaterialPoint& pt)
{
    FEUncoupledMaterial* m_pMat = dynamic_cast<FEUncoupledMaterial*>(GetParent());
    FEUncoupledMaterial* m_pBase = dynamic_cast<FEUncoupledMaterial*>(m_pMat->GetElasticMaterial());
    double sed = m_pBase->DevStrainEnergyDensity(pt);
    
    // clean up round-off errors
    if (sed < 0) sed = 0;
    
    return sqrt(2*sed);
}

//-----------------------------------------------------------------------------
// Strain energy density damage criterion
double FEDamageCriterionUCSED::DamageCriterion(FEMaterialPoint& pt)
{
    FEUncoupledMaterial* m_pMat = dynamic_cast<FEUncoupledMaterial*>(GetParent());
    FEUncoupledMaterial* m_pBase = dynamic_cast<FEUncoupledMaterial*>(m_pMat->GetElasticMaterial());
    double sed = m_pBase->DevStrainEnergyDensity(pt);
    
    return sed;
}

//-----------------------------------------------------------------------------
// von Mises stress damage criterion
double FEDamageCriterionUCVMS::DamageCriterion(FEMaterialPoint& pt)
{
    FEUncoupledMaterial* m_pMat = dynamic_cast<FEUncoupledMaterial*>(GetParent());
    FEUncoupledMaterial* m_pBase = dynamic_cast<FEUncoupledMaterial*>(m_pMat->GetElasticMaterial());
    mat3ds s = m_pBase->DevStress(pt);
    
    double vms = sqrt((SQR(s.xx()-s.yy()) + SQR(s.yy()-s.zz()) + SQR(s.zz()-s.xx())
                       + 6*(SQR(s.xy()) + SQR(s.yz()) + SQR(s.xz())))/2);
    
    return vms;
}

//-----------------------------------------------------------------------------
// max shear stress damage criterion
double FEDamageCriterionUCMSS::DamageCriterion(FEMaterialPoint& pt)
{
    // evaluate stress tensor
    FEUncoupledMaterial* m_pMat = dynamic_cast<FEUncoupledMaterial*>(GetParent());
    FEUncoupledMaterial* m_pBase = dynamic_cast<FEUncoupledMaterial*>(m_pMat->GetElasticMaterial());
    mat3ds s = m_pBase->DevStress(pt);
    
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
double FEDamageCriterionUCMNS::DamageCriterion(FEMaterialPoint& pt)
{
    // evaluate stress tensor
    FEUncoupledMaterial* m_pMat = dynamic_cast<FEUncoupledMaterial*>(GetParent());
    FEUncoupledMaterial* m_pBase = dynamic_cast<FEUncoupledMaterial*>(m_pMat->GetElasticMaterial());
    mat3ds s = m_pBase->DevStress(pt);
    
    // evaluate principal normal stresses
    double ps[3];
    s.eigen2(ps);
    
    // evaluate max normal stress
    double mns = max(max(ps[0],ps[1]),ps[2]);
    
    return mns;
}

//-----------------------------------------------------------------------------
// max normal Lagrange strain damage criterion
double FEDamageCriterionUCMNLS::DamageCriterion(FEMaterialPoint& mp)
{
    // evaluate strain tensor
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    mat3ds E = (pt.DevRightCauchyGreen() - mat3dd(1))/2;
    
    // evaluate principal normal strains
    double ps[3];
    E.eigen2(ps);
    
    // evaluate max normal Lagrange strain
    double mnls = max(max(ps[0],ps[1]),ps[2]);
    
    return mnls;
}

//-----------------------------------------------------------------------------
// max principal stretch ratio damage criterion
double FEDamageCriterionUCMPSR::DamageCriterion(FEMaterialPoint& mp)
{
    // evaluate strain tensor
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    mat3ds C = pt.DevRightCauchyGreen();
    
    // evaluate principal square stretches
    double ps[3];
    C.eigen2(ps);
    
    // evaluate max normal Lagrange strain
    double mpsr = sqrt(max(max(ps[0],ps[1]),ps[2]));
    
    return mpsr;
}
