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
#include "FEReactionRateHuiskesLazyZone.h"
#include "FEBiphasic.h"
#include "FEMultiphasic.h"
#include <FEBioMech/FERemodelingElasticMaterial.h>
#include <FEBioMech/FEElasticMixture.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>

// Material parameters for the FEMultiphasic material
BEGIN_FECORE_CLASS(FEReactionRateHuiskesLazyZone, FEReactionRate)
ADD_PARAMETER(m_B, "B");
ADD_PARAMETER(m_psi0lb, FE_RANGE_GREATER_OR_EQUAL(0.0), "psi0lb")->setUnits(UNIT_SPECIFIC_ENERGY)->setLongName("lower bound threshold");
ADD_PARAMETER(m_psi0ub, FE_RANGE_GREATER_OR_EQUAL(0.0), "psi0ub")->setUnits(UNIT_SPECIFIC_ENERGY)->setLongName("upper bound threshold");
ADD_PARAMETER(m_D, FE_RANGE_GREATER_OR_EQUAL(0.0), "D")->setUnits(UNIT_LENGTH)->setLongName("sensor distance");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEReactionRateHuiskesLazyZone::FEReactionRateHuiskesLazyZone(FEModel* pfem) : FEReactionRate(pfem)
{
    m_B = 0;
    m_psi0lb = 0;
    m_psi0ub = 0;
    m_D = 0;
    m_comp = -1;
    m_binit = false;
    m_M = 0;
    m_lsbm = -1;
}

//-----------------------------------------------------------------------------
//! initialization
bool FEReactionRateHuiskesLazyZone::Init()
{
    if (m_binit) return true;
    FEMultiphasic* pbm = dynamic_cast<FEMultiphasic*>(GetAncestor());
    if (pbm == nullptr) return true;    // in case material is not multiphasic
    
    FEChemicalReaction* pcm = dynamic_cast<FEChemicalReaction*>(m_pReact);
    int sbm = pcm->m_vPtmp[0]->m_speciesID;
    m_lsbm = pbm->FindLocalSBMID(sbm);
    m_M = pbm->SBMMolarMass(m_lsbm);
    
    // get neighboring elements for given proximity
    if (m_D > 0) {
        double mult = 4;    //! multiplier of characteristic distance, such that exp(-mult) << 1
        FEMesh& mesh = GetFEModel()->GetMesh();
        if (m_topo.Create(&mesh) == false)
        {
            feLogError("Failed building mesh topo.");
            return false;
        }
        feLogInfo("Evaluating element proximity...");
        m_EPL.assign(mesh.Elements(), std::vector<int>());
        for (int i=0; i< mesh.Elements(); ++i) {
            std::vector<int> epl = m_topo.ElementProximityList(i, m_D*mult, false);
            m_EPL[i] = epl;
        }
        feLogInfo("Done.");
    }
    m_binit = true;
    
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
    
    return true;
}

//-----------------------------------------------------------------------------
//! reaction rate at material point
double FEReactionRateHuiskesLazyZone::ReactionRate(FEMaterialPoint& pt)
{
    FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());
    double rhor = pbm->SolidReferentialApparentDensity(pt);
    double phir = pbm->SolidReferentialVolumeFraction(pt);
    
    FERemodelingMaterialPoint* rpt = nullptr;
    double J = 1;
    double sed = 0;
    if (m_comp == -1) {
        rpt = pt.ExtractData<FERemodelingMaterialPoint>();
        FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
        J = et.m_J;
    }
    else {
        FEElasticMixtureMaterialPoint* mp = pt.ExtractData<FEElasticMixtureMaterialPoint>();
        FEMaterialPoint& mpi = *mp->GetPointData(m_comp);
        rpt = mpi.ExtractData<FERemodelingMaterialPoint>();
        FEElasticMaterialPoint& et = *mpi.ExtractData<FEElasticMaterialPoint>();
        J = et.m_J;
    }
    sed = rpt->m_sed;
    double B = m_B(pt);
    double psi0lb = m_psi0lb(pt);
    double psi0ub = m_psi0ub(pt);
    double psi = sed/rhor;
    double dpsi = 0;
    if (psi < psi0lb) dpsi = psi - psi0lb;
    else if (psi > psi0ub) dpsi = psi - psi0ub;

    double zhat = B*dpsi/(J-phir)/m_M;
    // nonlocal contributions are weighted by exp(-d/D) and the weighted sum is
    // normalized by the sum of weights W, so that D sets the averaging distance
    // of the stimulus without changing the magnitude of the response
    double W = 1;
    if (m_D > 0) {
        FEMesh& mesh = GetFEModel()->GetMesh();
        // m_EPL is indexed by the mesh-wide element index (not the domain-local ID)
        int ie = m_topo.GetElementIndexFromID(pt.m_elem->GetID());
        int NEPL = (int)m_EPL[ie].size();
        // serial loop: zhat, W and the temporaries below are shared accumulators
        for (int i=0; i<NEPL; ++i) {
            int je = m_EPL[ie][i];
            if (je > -1) {
                FEElement* el = mesh.Element(je);
                for (int k=0; k<el->GaussPoints(); ++k) {
                    FEMaterialPoint& mp = *(el->GetMaterialPoint(k));
                    if (&mp == &pt) continue;   // this point's own contribution is added above
                    double w = exp(-(pt.m_rt - mp.m_rt).norm()/m_D);
                    W += w;
                    if (m_comp == -1) {
                        rpt = mp.ExtractData<FERemodelingMaterialPoint>();
                    }
                    else {
                        FEElasticMixtureMaterialPoint* mmp = mp.ExtractData<FEElasticMixtureMaterialPoint>();
                        FEMaterialPoint& mpi = *mmp->GetPointData(m_comp);
                        rpt = mpi.ExtractData<FERemodelingMaterialPoint>();
                    }
                    // the neighbor's stimulus contributes to the mass supply at this point,
                    // so it is normalized by the fluid volume (J-phir) of this point
                    double psin = rpt->m_sed/pbm->SolidReferentialApparentDensity(mp);
                    double lbn = m_psi0lb(mp);
                    double ubn = m_psi0ub(mp);
                    double dpsin = 0;
                    if (psin < lbn) dpsin = psin - lbn;
                    else if (psin > ubn) dpsin = psin - ubn;
                    zhat += w*m_B(mp)*dpsin/(J-phir)/m_M;
                }
            }
        }
    }
    zhat /= W;
    return zhat;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
mat3ds FEReactionRateHuiskesLazyZone::Tangent_ReactionRate_Strain(FEMaterialPoint& pt)
{
    FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());
    
    FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    
    double rhor = pbm->SolidReferentialApparentDensity(pt);
    double phir = pbm->SolidReferentialVolumeFraction(pt);
    double p = pbm->GetActualFluidPressure(pt);
    
    double J = 1;
    if (m_comp == -1) {
        FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
        J = et.m_J;
    }
    else {
        FEElasticMixtureMaterialPoint* mmp = pt.ExtractData<FEElasticMixtureMaterialPoint>();
        FEMaterialPoint& mpi = *mmp->GetPointData(m_comp);
        FEElasticMaterialPoint& et = *mpi.ExtractData<FEElasticMaterialPoint>();
        J = et.m_J;
    }
    double zhat = ReactionRate(pt);
    mat3dd I(1);
    double B = m_B(pt);
    mat3ds dzhatde = mat3dd(-zhat/(J-phir));
    // the stimulus at this point only contributes outside the lazy zone
    FERemodelingMaterialPoint* rpt = (m_comp == -1) ? pt.ExtractData<FERemodelingMaterialPoint>() :
        pt.ExtractData<FEElasticMixtureMaterialPoint>()->GetPointData(m_comp)->ExtractData<FERemodelingMaterialPoint>();
    double psi = rpt->m_sed/rhor;
    if ((psi < m_psi0lb(pt)) || (psi > m_psi0ub(pt)))
        dzhatde += (et.m_s+I*p)*(B/rhor/m_M/(J-phir)/ProximityWeight(pt));
    // Neighboring points contribute to the reaction rate (see ReactionRate), but their
    // contributions depend on the strains at those points, not at this point, so they
    // do not contribute to this tangent (this cross-element coupling is neglected).
    return dzhatde;
}

//-----------------------------------------------------------------------------
//! sum of the proximity weights exp(-d/D) of all integration points that contribute
//! to the reaction rate at this point (including this point, whose weight is 1)
double FEReactionRateHuiskesLazyZone::ProximityWeight(FEMaterialPoint& pt)
{
    double W = 1;
    if (m_D > 0) {
        FEMesh& mesh = GetFEModel()->GetMesh();
        int ie = m_topo.GetElementIndexFromID(pt.m_elem->GetID());
        for (int je : m_EPL[ie]) {
            if (je < 0) continue;
            FEElement* el = mesh.Element(je);
            for (int k=0; k<el->GaussPoints(); ++k) {
                FEMaterialPoint& mp = *(el->GetMaterialPoint(k));
                if (&mp == &pt) continue;
                W += exp(-(pt.m_rt - mp.m_rt).norm()/m_D);
            }
        }
    }
    return W;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with effective fluid pressure at material point
double FEReactionRateHuiskesLazyZone::Tangent_ReactionRate_Pressure(FEMaterialPoint& pt)
{
    return 0;
}

