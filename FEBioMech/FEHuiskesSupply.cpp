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
#include "FEHuiskesSupply.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEHuiskesSupply, FESolidSupply)
	ADD_PARAMETER(m_B, "B");
	ADD_PARAMETER(m_k, "k")->setUnits(UNIT_SPECIFIC_ENERGY);
    ADD_PARAMETER(m_D, FE_RANGE_GREATER_OR_EQUAL(0.0), "D")->setUnits(UNIT_LENGTH)->setLongName("sensor distance");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEHuiskesSupply::FEHuiskesSupply(FEModel* pfem) : FESolidSupply(pfem)
{
	m_B = m_k = m_D = 0;
}

//-----------------------------------------------------------------------------
//! Initialization
bool FEHuiskesSupply::Init()
{
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
    
    return true;
}

//-----------------------------------------------------------------------------
//! Solid supply
double FEHuiskesSupply::Supply(FEMaterialPoint& pt)
{
	FERemodelingMaterialPoint* rpt = pt.ExtractData<FERemodelingMaterialPoint>();
	double rhor = rpt->m_rhor;
	double sed = rpt->m_sed;
	double rhorhat = m_B*(sed/rhor - m_k);
    
    // nonlocal contributions are weighted by exp(-d/D) and the weighted sum is
    // normalized by the sum of weights W, so that D sets the averaging distance
    // of the stimulus without changing the magnitude of the response
    double W = 1;
    if (m_D > 0) {
        FEMesh& mesh = GetFEModel()->GetMesh();
        // m_EPL is indexed by the mesh-wide element index (not the domain-local ID)
        int ie = m_topo.GetElementIndexFromID(pt.m_elem->GetID());
        int NEPL = (int)m_EPL[ie].size();
        // serial loop: rhorhat and W are shared accumulators
        for (int i=0; i<NEPL; ++i) {
            int je = m_EPL[ie][i];
            if (je > -1) {
                FEElement* el = mesh.Element(je);
                for (int k=0; k<el->GaussPoints(); ++k) {
                    FEMaterialPoint& mp = *(el->GetMaterialPoint(k));
                    if (&mp == &pt) continue;   // this point's own contribution is added above
                    FERemodelingMaterialPoint* rptn = mp.ExtractData<FERemodelingMaterialPoint>();
                    if (rptn == nullptr) continue;
                    double w = exp(-(pt.m_rt - mp.m_rt).norm()/m_D);
                    W += w;
                    rhorhat += w*m_B*(rptn->m_sed/rptn->m_rhor - m_k);
                }
            }
        }
    }

    return rhorhat/W;
}

//-----------------------------------------------------------------------------
//! Tangent of solid supply with respect to strain
mat3ds FEHuiskesSupply::Tangent_Supply_Strain(FEMaterialPoint &pt)
{
    FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
	FERemodelingMaterialPoint* rpt = pt.ExtractData<FERemodelingMaterialPoint>();
    mat3ds ruhat = et.m_s*(m_B/rpt->m_rhor/ProximityWeight(pt));

    // Neighboring points contribute to the supply (see Supply), but their contributions
    // depend on the strains at those points, not at this point, so they do not
    // contribute to this tangent (this cross-element coupling is neglected).

    return ruhat;
}

//-----------------------------------------------------------------------------
//! Tangent of solid supply with respect to referential density
double FEHuiskesSupply::Tangent_Supply_Density(FEMaterialPoint &mp)
{
	FERemodelingMaterialPoint& rpt = *mp.ExtractData<FERemodelingMaterialPoint>();
    double rhor = rpt.m_rhor;
    double sed = rpt.m_sed;
    double dsed = rpt.m_dsed;
	return (dsed - sed/rhor)*m_B/rhor/ProximityWeight(mp);
}

//-----------------------------------------------------------------------------
//! sum of the proximity weights exp(-d/D) of all integration points that contribute
//! to the supply at this point (including this point, whose weight is 1)
double FEHuiskesSupply::ProximityWeight(FEMaterialPoint& pt)
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
                if (mp.ExtractData<FERemodelingMaterialPoint>() == nullptr) continue;
                W += exp(-(pt.m_rt - mp.m_rt).norm()/m_D);
            }
        }
    }
    return W;
}
