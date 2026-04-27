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
#include "FEReactionRateHuiskes.h"
#include "FEBiphasic.h"
#include "FEMultiphasic.h"
#include <FEBioMech/FERemodelingElasticMaterial.h>
#include <FEBioMech/FEElasticMixture.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>

// Material parameters for the FEMultiphasic material
BEGIN_FECORE_CLASS(FEReactionRateHuiskes, FEReactionRate)
	ADD_PARAMETER(m_B, "B");
	ADD_PARAMETER(m_psi0, FE_RANGE_GREATER_OR_EQUAL(0.0), "psi0")->setUnits(UNIT_SPECIFIC_ENERGY);
    ADD_PARAMETER(m_D, FE_RANGE_GREATER_OR_EQUAL(0.0), "D")->setUnits(UNIT_LENGTH)->setLongName("sensor distance");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEReactionRateHuiskes::FEReactionRateHuiskes(FEModel* pfem) : FEReactionRate(pfem) 
{ 
    m_B = 0;
    m_psi0 = 0;
    m_D = 0;
    m_comp = -1;
    m_binit = false;
    m_M = 0;
    m_lsbm = -1;
}

//-----------------------------------------------------------------------------
//! initialization
bool FEReactionRateHuiskes::Init()
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
            std::vector<int> epl = m_topo.ElementProximityList(i, m_D*mult);
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
double FEReactionRateHuiskes::ReactionRate(FEMaterialPoint& pt)
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
    double psi0 = m_psi0(pt);
	double zhat = B*(sed/rhor - psi0)/(J-phir)/m_M;
    if (m_D > 0) {
        FEMesh& mesh = GetFEModel()->GetMesh();
        int ie = pt.m_elem->GetLocalID();
        int NEPL = (int)m_EPL[ie].size();
#pragma omp parallel for shared (NEPL)
        for (int i=0; i<NEPL; ++i) {
            int je = m_EPL[ie][i];
            if (je > -1) {
                FEElement* el = mesh.Element(je);
                for (int k=0; k<el->GaussPoints(); ++k) {
                    FEMaterialPoint& mp = *(el->GetMaterialPoint(k));
                    double d = (pt.m_rt - mp.m_rt).unit();
                    if (m_comp == -1) {
                        rpt = mp.ExtractData<FERemodelingMaterialPoint>();
                        FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
                        J = et.m_J;
                    }
                    else {
                        FEElasticMixtureMaterialPoint* mmp = mp.ExtractData<FEElasticMixtureMaterialPoint>();
                        FEMaterialPoint& mpi = *mmp->GetPointData(m_comp);
                        rpt = mpi.ExtractData<FERemodelingMaterialPoint>();
                    }
                    sed = rpt->m_sed;
                    rhor = pbm->SolidReferentialApparentDensity(mp);
                    phir = pbm->SolidReferentialVolumeFraction(mp);
                    B = m_B(mp);
                    psi0 = m_psi0(mp);
                    zhat += exp(-d/m_D)*B*(sed/rhor - psi0)/(J-phir)/m_M;
                }
            }
        }
    }
	return zhat;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
mat3ds FEReactionRateHuiskes::Tangent_ReactionRate_Strain(FEMaterialPoint& pt)
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
    mat3ds dzhatde = (I*(-zhat) + (et.m_s+I*p)*(B/rhor/m_M))/(J-phir);
    if (m_D > 0) {
        mat3ds s(0);
        FEMesh& mesh = GetFEModel()->GetMesh();
        int ie = pt.m_elem->GetLocalID();
        int NEPL = (int)m_EPL[ie].size();
#pragma omp parallel for shared (NEPL)
        for (int i=0; i<NEPL; ++i) {
            int je = m_EPL[ie][i];
            if (je > -1) {
                FEElement* el = mesh.Element(je);
                for (int k=0; k<el->GaussPoints(); ++k) {
                    FEMaterialPoint& mp = *(el->GetMaterialPoint(k));
                    double d = (pt.m_rt - mp.m_rt).unit();
                    if (m_comp == -1) {
                        FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
                        // we just want the effective stress, so subtract the pressure term
                        p = pbm->GetActualFluidPressure(mp);
                        s = et.m_s+I*p;
                    }
                    else {
                        FEElasticMixtureMaterialPoint* mmp = mp.ExtractData<FEElasticMixtureMaterialPoint>();
                        FEMaterialPoint& mpi = *mmp->GetPointData(m_comp);
                        FEElasticMaterialPoint& et = *mpi.ExtractData<FEElasticMaterialPoint>();
                        // the stress stored in elastic mixture material points is an effective stress
                        s = et.m_s;
                    }
                    rhor = pbm->SolidReferentialApparentDensity(mp);
                    B = m_B(mp);
                    dzhatde += exp(-d/m_D)*s*(B/rhor/m_M)/(J-phir);
                }
            }
        }
    }
	return dzhatde;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with effective fluid pressure at material point
double FEReactionRateHuiskes::Tangent_ReactionRate_Pressure(FEMaterialPoint& pt)
{
	return 0;
}

