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
            std::vector<int> epl = m_topo.ElementProximityList(i, m_D*mult);
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
                    rpt = mp.ExtractData<FERemodelingMaterialPoint>();
                    FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
                    rhorhat += exp(-d/m_D)*m_B*(rpt->m_sed/rpt->m_rhor - m_k);
                }
            }
        }
    }

    return rhorhat;
}

//-----------------------------------------------------------------------------
//! Tangent of solid supply with respect to strain
mat3ds FEHuiskesSupply::Tangent_Supply_Strain(FEMaterialPoint &pt)
{
    FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
	FERemodelingMaterialPoint* rpt = pt.ExtractData<FERemodelingMaterialPoint>();
    mat3ds ruhat = et.m_s*(m_B/rpt->m_rhor);

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
                    rpt = mp.ExtractData<FERemodelingMaterialPoint>();
                    FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
                    ruhat += et.m_s*(exp(-d/m_D)*m_B/rpt->m_rhor);
                }
            }
        }
    }

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
	return (dsed - sed/rhor)*m_B/rhor;
}

