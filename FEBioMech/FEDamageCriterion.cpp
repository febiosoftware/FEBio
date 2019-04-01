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

#include "stdafx.h"
#include "FEDamageCriterion.h"
#include "FEDamageMaterial.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
// Simo's damage criterion uses sqrt(2*strain energy density)
double FEDamageCriterionSimo::DamageCriterion(FEMaterialPoint& pt)
{
    FEDamageMaterial* m_pMat = dynamic_cast<FEDamageMaterial*>(GetParent());
	FEElasticMaterial* m_pBase = m_pMat->GetElasticMaterial();
    double sed = m_pBase->StrainEnergyDensity(pt);
    
    // clean up round-off errors
    if (sed < 0) sed = 0;
    
    return sqrt(2*sed);
}

//-----------------------------------------------------------------------------
// Strain energy density damage criterion
double FEDamageCriterionSED::DamageCriterion(FEMaterialPoint& pt)
{
	FEDamageMaterial* m_pMat = dynamic_cast<FEDamageMaterial*>(GetParent());
	FEElasticMaterial* m_pBase = m_pMat->GetElasticMaterial();
	double sed = m_pBase->StrainEnergyDensity(pt);
    
    return sed;
}

//-----------------------------------------------------------------------------
// Specific strain energy damage criterion
double FEDamageCriterionSSE::DamageCriterion(FEMaterialPoint& pt)
{
	FEDamageMaterial* m_pMat = dynamic_cast<FEDamageMaterial*>(GetParent());
	FEElasticMaterial* m_pBase = m_pMat->GetElasticMaterial();
	double sed = m_pBase->StrainEnergyDensity(pt);
    FEParamDouble& rho = m_pBase->Density();
    
    return sed/rho(pt);
}

//-----------------------------------------------------------------------------
// von Mises stress damage criterion
double FEDamageCriterionVMS::DamageCriterion(FEMaterialPoint& pt)
{
	FEDamageMaterial* m_pMat = dynamic_cast<FEDamageMaterial*>(GetParent());
	FEElasticMaterial* m_pBase = m_pMat->GetElasticMaterial();
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
	FEDamageMaterial* m_pMat = dynamic_cast<FEDamageMaterial*>(GetParent());
	FEElasticMaterial* m_pBase = m_pMat->GetElasticMaterial();
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
	FEDamageMaterial* m_pMat = dynamic_cast<FEDamageMaterial*>(GetParent());
	FEElasticMaterial* m_pBase = m_pMat->GetElasticMaterial();
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
