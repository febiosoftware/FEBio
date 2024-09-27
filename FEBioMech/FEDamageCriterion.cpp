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



#include "stdafx.h"
#include "FEDamageCriterion.h"
#include "FEDamageMaterial.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <algorithm>

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEDamageCriterionSimo, FEDamageCriterion)
    ADD_PROPERTY(m_nla, "nonlocal_averaging",FEProperty::Optional)->SetLongName("nonlocal averaging");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// Simo's damage criterion uses sqrt(2*strain energy density)
double FEDamageCriterionSimo::DCpt(FEMaterialPoint& pt)
{
    FEElasticMaterial* m_pMat = dynamic_cast<FEElasticMaterial*>(GetParent());
	FEElasticMaterial* m_pBase = m_pMat->GetElasticMaterial();
    double sed = m_pBase->StrainEnergyDensity(pt);
    
    // clean up round-off errors
    if (sed < 0) sed = 0;
    
    return sqrt(2*sed);
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEDamageCriterionSED, FEDamageCriterion)
ADD_PROPERTY(m_nla, "nonlocal_averaging",FEProperty::Optional)->SetLongName("nonlocal averaging");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// Strain energy density damage criterion
double FEDamageCriterionSED::DCpt(FEMaterialPoint& pt)
{
	FEElasticMaterial* m_pMat = dynamic_cast<FEElasticMaterial*>(GetParent());
	FEElasticMaterial* m_pBase = m_pMat->GetElasticMaterial();
	double sed = m_pBase->StrainEnergyDensity(pt);
    
    return sed;
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEDamageCriterionSSE, FEDamageCriterion)
ADD_PROPERTY(m_nla, "nonlocal_averaging",FEProperty::Optional)->SetLongName("nonlocal averaging");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// Specific strain energy damage criterion
double FEDamageCriterionSSE::DCpt(FEMaterialPoint& pt)
{
	FEElasticMaterial* m_pMat = dynamic_cast<FEElasticMaterial*>(GetParent());
	FEElasticMaterial* m_pBase = m_pMat->GetElasticMaterial();
	double sed = m_pBase->StrainEnergyDensity(pt);
    
    return sed/ m_pBase->Density(pt);
}

//-----------------------------------------------------------------------------
// von Mises stress damage criterion
double FEDamageCriterionVMS::DCpt(FEMaterialPoint& pt)
{
	FEElasticMaterial* m_pMat = dynamic_cast<FEElasticMaterial*>(GetParent());
	FEElasticMaterial* m_pBase = m_pMat->GetElasticMaterial();
	mat3ds s = m_pBase->Stress(pt);
    
    double vms = sqrt((SQR(s.xx()-s.yy()) + SQR(s.yy()-s.zz()) + SQR(s.zz()-s.xx())
                       + 6*(SQR(s.xy()) + SQR(s.yz()) + SQR(s.xz())))/2);
    
    return vms;
}


//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEDamageCriterionVMS, FEDamageCriterion)
    ADD_PROPERTY(m_nla, "nonlocal_averaging",FEProperty::Optional)->SetLongName("nonlocal averaging");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! criterion tangent with respect to stress
mat3ds FEDamageCriterionVMS::CSTpt(FEMaterialPoint& pt)
{
    FEElasticMaterial* m_pMat = dynamic_cast<FEElasticMaterial*>(GetParent());
    FEElasticMaterial* m_pBase = m_pMat->GetElasticMaterial();
    mat3ds s = m_pBase->Stress(pt);
    double vms = DamageCriterion(pt);
    return s.dev()*(1.5/vms);
}

//-----------------------------------------------------------------------------
// Drucker yield criterion
BEGIN_FECORE_CLASS(FEDamageCriterionDrucker, FEDamageCriterion)
    ADD_PARAMETER(m_c, FE_RANGE_CLOSED(-27.0/8.0,9.0/4.0), "c");
    ADD_PROPERTY(m_nla, "nonlocal_averaging",FEProperty::Optional)->SetLongName("nonlocal averaging");
END_FECORE_CLASS();

double FEDamageCriterionDrucker::DCpt(FEMaterialPoint& pt)
{
    FEElasticMaterial* m_pMat = dynamic_cast<FEElasticMaterial*>(GetParent());
    FEElasticMaterial* m_pBase = m_pMat->GetElasticMaterial();
    mat3ds s = m_pBase->Stress(pt);
    mat3ds sdev = s.dev();
    double J2 = (sdev*sdev).trace()/2;
    double J3 = sdev.det();
    double c = m_c(pt);
    double tau = pow(pow(J2,3) - c*pow(J3,2),1./6.);
    
    return tau;
}

//-----------------------------------------------------------------------------
//! criterion tangent with respect to stress
mat3ds FEDamageCriterionDrucker::CSTpt(FEMaterialPoint& pt)
{
    FEElasticMaterial* m_pMat = dynamic_cast<FEElasticMaterial*>(GetParent());
    FEElasticMaterial* m_pBase = m_pMat->GetElasticMaterial();
    mat3ds s = m_pBase->Stress(pt);
    mat3ds sdev = s.dev();
    double J2 = (sdev*sdev).trace()/2;
    double J3 = sdev.det();
    double c = m_c(pt);
    double tau = pow(pow(J2,3) - c*pow(J3,2),1./6.);
    return (s.dev()*(pow(J2,2)/2) - (sdev.inverse()).dev()*(pow(J3, 2)*c/3))/pow(tau, 5);
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEDamageCriterionMSS, FEDamageCriterion)
    ADD_PROPERTY(m_nla, "nonlocal_averaging",FEProperty::Optional)->SetLongName("nonlocal averaging");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// max shear stress damage criterion
double FEDamageCriterionMSS::DCpt(FEMaterialPoint& pt)
{
    // evaluate stress tensor
	FEElasticMaterial* m_pMat = dynamic_cast<FEElasticMaterial*>(GetParent());
	FEElasticMaterial* m_pBase = m_pMat->GetElasticMaterial();
	mat3ds s = m_pBase->Stress(pt);
    
    // evaluate principal normal stresses
    double ps[3];
    s.eigen2(ps);   // sorted in ascending order
    
    // evaluate max shear stresses
    double mss = (ps[2] - ps[0])/2;
    
    return mss;
}

//-----------------------------------------------------------------------------
//! criterion tangent with respect to stress
mat3ds FEDamageCriterionMSS::CSTpt(FEMaterialPoint& pt)
{
    FEElasticMaterial* m_pMat = dynamic_cast<FEElasticMaterial*>(GetParent());
    FEElasticMaterial* m_pBase = m_pMat->GetElasticMaterial();
    mat3ds s = m_pBase->Stress(pt);
    // evaluate principal normal stresses
    double ps[3];
    vec3d v[3];
    s.eigen2(ps,v); // sorted in ascending order
    // clean up small differences
    const double eps = 1e-6;
    if (fabs(ps[2] - ps[0]) <= eps*fabs(ps[2])) ps[0] = ps[2];
    if (fabs(ps[2] - ps[1]) <= eps*fabs(ps[2])) ps[1] = ps[2];
    
    // return tangent
    mat3ds N;
    if (ps[2] > ps[1]) {
        if (ps[1] > ps[0])
            N = (dyad(v[2])-dyad(v[0]))/2;
        else
            N = (dyad(v[2])-(dyad(v[1]) + dyad(v[0]))/2)/2;
    }
    else if (ps[1] > ps[0])
        N = ((dyad(v[2]) + dyad(v[1]))/2 - dyad(v[0]))/2;
    else N.zero();
    
    return N;
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEDamageCriterionMNS, FEDamageCriterion)
    ADD_PROPERTY(m_nla, "nonlocal_averaging",FEProperty::Optional)->SetLongName("nonlocal averaging");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// max normal stress damage criterion
double FEDamageCriterionMNS::DCpt(FEMaterialPoint& pt)
{
    // evaluate stress tensor
	FEElasticMaterial* m_pMat = dynamic_cast<FEElasticMaterial*>(GetParent());
	FEElasticMaterial* m_pBase = m_pMat->GetElasticMaterial();
	mat3ds s = m_pBase->Stress(pt);
    
    // evaluate principal normal stresses
    double ps[3];
    s.eigen2(ps);   // sorted in ascending order
    
    return ps[2];
}

//-----------------------------------------------------------------------------
//! criterion tangent with respect to stress
mat3ds FEDamageCriterionMNS::CSTpt(FEMaterialPoint& pt)
{
    FEElasticMaterial* m_pMat = dynamic_cast<FEElasticMaterial*>(GetParent());
    FEElasticMaterial* m_pBase = m_pMat->GetElasticMaterial();
    mat3ds s = m_pBase->Stress(pt);
    // evaluate principal normal stresses
    double ps[3];
    vec3d v[3];
    s.eigen2(ps,v); // sorted in ascending order
    
    // clean up small differences
    const double eps = 1e-6;
    if (fabs(ps[2] - ps[0]) <= eps*fabs(ps[2])) ps[0] = ps[2];
    if (fabs(ps[2] - ps[1]) <= eps*fabs(ps[2])) ps[1] = ps[2];
    
    // return tangent
    mat3ds N;
    if (ps[2] > ps[1]) N = dyad(v[2]);
    else if (ps[2] > ps[0]) N = dyad(v[2]) + dyad(v[1]);
    else N = mat3dd(1);
    
    return N;
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEDamageCriterionMNLS, FEDamageCriterion)
    ADD_PROPERTY(m_nla, "nonlocal_averaging",FEProperty::Optional)->SetLongName("nonlocal averaging");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// max normal Lagrange strain damage criterion
double FEDamageCriterionMNLS::DCpt(FEMaterialPoint& mp)
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

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEDamageCriterionOSS, FEDamageCriterion)
    ADD_PROPERTY(m_nla, "nonlocal_averaging",FEProperty::Optional)->SetLongName("nonlocal averaging");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// octahedral shear strain damage criterion
double FEDamageCriterionOSS::DCpt(FEMaterialPoint& mp)
{
    // evaluate strain tensor
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    mat3ds E = pt.Strain();
    
    mat3ds devE = E.dev();
    double oss = sqrt(devE.dotdot(devE)*(2./3.));
    
    return oss;
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEDamageCriterionONS, FEDamageCriterion)
    ADD_PROPERTY(m_nla, "nonlocal_averaging",FEProperty::Optional)->SetLongName("nonlocal averaging");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// octahedral natural shear strain damage criterion
double FEDamageCriterionONS::DCpt(FEMaterialPoint& mp)
{
    // evaluate strain tensor
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    mat3ds h = pt.LeftHencky();
    
    mat3ds devh = h.dev();
    double ons = sqrt(devh.dotdot(devh)*(2./3.));
    
    return ons;
}
