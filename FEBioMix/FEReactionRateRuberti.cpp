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
#include "FEReactionRateRuberti.h"
#include "FEBiphasic.h"
#include "FEMultiphasic.h"
#include "FESoluteInterface.h"
#include <FEBioMech/FERemodelingElasticMaterial.h>
#include <FEBioMech/FEElasticMixture.h>

// Material parameters for the FEMultiphasic material
BEGIN_FECORE_CLASS(FEReactionRateRuberti, FEReactionRate)
	ADD_PARAMETER(m_kd, "kd");
	ADD_PARAMETER(m_sig, FE_RANGE_GREATER_OR_EQUAL(0.0), "sigma");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEReactionRateRuberti::FEReactionRateRuberti(FEModel* pfem) : FEReactionRate(pfem)
{
    m_kd = 0;
    m_sig = 0;
    m_comp = -1;
    m_fiber = nullptr;
    m_M = 0;
    m_lsbm = -1;
}

//-----------------------------------------------------------------------------
//! initialization
bool FEReactionRateRuberti::Init()
{
    FEMultiphasic* pbm = dynamic_cast<FEMultiphasic*>(GetAncestor());
    if (pbm == nullptr) return true;    // in case material is not multiphasic

    FESoluteInterface* psi = dynamic_cast<FESoluteInterface*>(GetAncestor());
    
    FEChemicalReaction* pcm = dynamic_cast<FEChemicalReaction*>(m_pReact);
    int sbm = pcm->m_vRtmp[0]->m_speciesID;
    m_lsbm = pbm->FindLocalSBMID(sbm);
    m_M = pbm->SBMMolarMass(m_lsbm);
    
    FEElasticMaterial* pem = pbm->GetSolid();
    FEElasticMixture* psm = dynamic_cast<FEElasticMixture*>(pem);
    if (psm == nullptr) return true;    // in case material is not a solid mixture

    for (int i=0; i<psm->Materials(); ++i) {
        pem = psm->GetMaterial(i);
        // check the types of materials that can employ this reaction rate
        FERemodelingInterface* pri = dynamic_cast<FERemodelingInterface*>(pem);
        if (pri) {
            if (sbm == pri->m_sbm) m_comp = pri->m_comp;
        }
    }
    if (m_comp == -1) return false;
    
    // get the fiber property (this remodeling rule only works with fibrous materials)
    m_fiber = dynamic_cast<FEVec3dValuator*>(psm->GetMaterial(m_comp)->GetProperty("fiber"));
    if (m_fiber == nullptr) return false;
    
    return true;
}

//-----------------------------------------------------------------------------
//! reaction rate at material point
double FEReactionRateRuberti::ReactionRate(FEMaterialPoint& pt)
{
    FEMultiphasic* pbm = dynamic_cast<FEMultiphasic*>(GetAncestor());
    FEBiphasicInterface* pbi = dynamic_cast<FEBiphasicInterface*>(GetAncestor());
    FESoluteInterface* psi = dynamic_cast<FESoluteInterface*>(GetAncestor());
    double phir = pbm->SolidReferentialVolumeFraction(pt);
    double c = psi->SBMConcentration(pt, m_lsbm);

    mat3d Q = pbm->GetLocalCS(pt);
    vec3d a0 = m_fiber->unitVector(pt);
    vec3d ar = Q * a0;
    FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    double lam = (et.m_F*ar).norm();
    double k = (lam > 1.0) ? m_kd(pt)*exp(-0.5*pow((lam-1)/m_sig(pt),2)) : m_kd(pt);
    
	return k;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
mat3ds FEReactionRateRuberti::Tangent_ReactionRate_Strain(FEMaterialPoint& pt)
{
    FEMultiphasic* pbm = dynamic_cast<FEMultiphasic*>(GetAncestor());
    FEBiphasicInterface* pbi = dynamic_cast<FEBiphasicInterface*>(GetAncestor());
    FESoluteInterface* psi = dynamic_cast<FESoluteInterface*>(GetAncestor());
    double phir = pbm->SolidReferentialVolumeFraction(pt);
    double c = psi->SBMConcentration(pt, m_lsbm);
    
    mat3d Q = pbm->GetLocalCS(pt);
    vec3d a0 = m_fiber->unitVector(pt);
    vec3d ar = Q * a0;
    FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    double lam = (et.m_F*ar).norm();
    double kd = m_kd(pt);
    double sig = m_sig(pt);
    double k = (lam > 1.0) ? kd*exp(-0.5*pow(lam-1,2)/sig) : kd;
    double dkdlam = (lam > 1.0) ? -(lam-1)/sig*kd*exp(-0.5*pow((lam-1)/sig,2)) : 0;

    mat3ds dzhatde = dyad(et.m_F*ar)*(dkdlam/lam/et.m_J);
	return dzhatde;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with effective fluid pressure at material point
double FEReactionRateRuberti::Tangent_ReactionRate_Pressure(FEMaterialPoint& pt)
{
	return 0;
}

