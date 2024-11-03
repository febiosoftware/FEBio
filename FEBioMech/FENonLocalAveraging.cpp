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
#include "FENonLocalAveraging.h"
#include "FEDamageCriterion.h"
#include "FEReactivePlasticityMaterialPoint.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <algorithm>

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif
#ifdef max
#undef max
#endif

//-----------------------------------------------------------------------------
bool FENonLocalAveraging::Init()
{
    if (m_binit) return true;
    
    if (m_krnl == nullptr) return false;
    
    // get neighboring elements for given proximity
    double R = m_krnl->m_R;
    if (R > 0) {
        // initialize the kernel integral to infinity
        m_krnl->m_c = m_krnl->KernelIntegralInfinity();

        double mult = m_krnl->m_mult;
        FEMesh& mesh = GetFEModel()->GetMesh();
        feLogInfo("Evaluating mesh topology...");
        if (m_topo.Create(&mesh) == false)
        {
            feLogError("Failed building mesh topo.");
            return false;
        }
        feLogInfo("Evaluating element proximity...");
        m_EPL.assign(mesh.Elements(), std::vector<int>());
        double eplmin = 1e9;
        double eplmax = 0;
        double eplavg = 0;
        int NE = mesh.Elements();
        for (int i=0; i< mesh.Elements(); ++i) {
            std::vector<int> epl = m_topo.ElementProximityList(i, R*mult, false);
            m_EPL[i] = epl;
            double epls = epl.size();
            eplmin = std::min(eplmin, epls);
            eplmax = std::max(eplmax, epls);
            eplavg += epls;
        }
        eplavg /= NE;
        feLogInfo("Done.");
        feLog("Number of neighboring elements:\n");
        feLog("-------------------------------\n");
        feLog("Min. = %g, Max. = %g, Avg. = %g\n", eplmin, eplmax, eplavg);
    }
    m_binit = true;
    return true;
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FENLABazant, FENonLocalAveraging)
    ADD_PROPERTY(m_krnl, "kernel");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! criterion average
double FENLABazant::DamageCriterionAverage(FEMaterialPoint& pt)
{
    double ca = 0;
    
    double R = m_krnl->m_R;
    if (R > 0) {
        FEDamageCriterion* pmat = dynamic_cast<FEDamageCriterion*>(GetParent());
        assert(pmat);
        FEMesh& mesh = GetFEModel()->GetMesh();
        int ie = pt.m_elem->GetLocalID();
        int NEPL = (int)m_EPL[ie].size();
        double V = 0;
        //#pragma omp parallel for shared (NEPL)
        for (int i=0; i<NEPL; ++i) {
            int je = m_EPL[ie][i];
            if (je > -1) {
                FEElement* el = mesh.Element(je);

                // skip inactive elements. Needed for this to work with erosion
                if(!el->isActive()) continue;

                for (int k=0; k<el->GaussPoints(); ++k) {
                    FEMaterialPoint& mp = *(el->GetMaterialPoint(k));
                    double krnl = m_krnl->Kernel(pt, mp);
                    ca += pmat->DCpt(mp)*krnl*mp.m_V0;
                    V += krnl*mp.m_V0;
                }
            }
        }
        if (V > 0) ca /= V;
    }
    
    return ca;

}

//-----------------------------------------------------------------------------
//! criterion stress tangent average
mat3ds FENLABazant::DamageCriterionTangentAverage(FEMaterialPoint& pt)
{
    mat3ds cst(0);
    
    double R = m_krnl->m_R;
    if (R > 0) {
        FEDamageCriterion* pmat = dynamic_cast<FEDamageCriterion*>(GetParent());
        assert(pmat);
        FEMesh& mesh = GetFEModel()->GetMesh();
        int ie = pt.m_elem->GetLocalID();
        int NEPL = (int)m_EPL[ie].size();
        double V = 0;
        //#pragma omp parallel for shared (NEPL)
        for (int i=0; i<NEPL; ++i) {
            int je = m_EPL[ie][i];
            if (je > -1) {
                FEElement* el = mesh.Element(je);

                // skip inactive elements. Needed for this to work with erosion
                if(!el->isActive()) continue;

                for (int k=0; k<el->GaussPoints(); ++k) {
                    FEMaterialPoint& mp = *(el->GetMaterialPoint(k));
                    double krnl = m_krnl->Kernel(pt, mp);
                    cst += pmat->CSTpt(mp)*krnl*mp.m_V0;
                    V += krnl*mp.m_V0;
                }
            }
        }
        if (V > 0) cst /= V;
    }
    
    return cst;
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FENLABorino, FENonLocalAveraging)
    ADD_PROPERTY(m_krnl, "kernel");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! criterion average
double FENLABorino::DamageCriterionAverage(FEMaterialPoint& pt)
{
    double ca = 0;
    
    double R = m_krnl->m_R;
    if (R > 0) {
        FEDamageCriterion* pmat = dynamic_cast<FEDamageCriterion*>(GetParent());
        assert(pmat);
        FEMesh& mesh = GetFEModel()->GetMesh();
        int ie = pt.m_elem->GetLocalID();
        int NEPL = (int)m_EPL[ie].size();
        double V = 0;
        //#pragma omp parallel for shared (NEPL)
        for (int i=0; i<NEPL; ++i) {
            int je = m_EPL[ie][i];
            if (je > -1) {
                FEElement* el = mesh.Element(je);

                // skip inactive elements. Needed for this to work with erosion
                if(!el->isActive()) continue;

                for (int k=0; k<el->GaussPoints(); ++k) {
                    FEMaterialPoint& mp = *(el->GetMaterialPoint(k));
                    double krnl = m_krnl->Kernel(pt, mp);
                    ca += pmat->DCpt(mp)*krnl*mp.m_V0;
                    V += krnl*mp.m_V0;
                }
            }
        }
        if (V > 0) {
            ca /= V;
            if (V < 1) ca += (1 - V)*pmat->DCpt(pt);
        }
    }
    
    return ca;
    
}

//-----------------------------------------------------------------------------
//! criterion stress tangent average
mat3ds FENLABorino::DamageCriterionTangentAverage(FEMaterialPoint& pt)
{
    mat3ds cst(0);
    
    double R = m_krnl->m_R;
    if (R > 0) {
        FEDamageCriterion* pmat = dynamic_cast<FEDamageCriterion*>(GetParent());
        assert(pmat);
        FEMesh& mesh = GetFEModel()->GetMesh();
        int ie = pt.m_elem->GetLocalID();
        int NEPL = (int)m_EPL[ie].size();
        double V = 0;
        //#pragma omp parallel for shared (NEPL)
        for (int i=0; i<NEPL; ++i) {
            int je = m_EPL[ie][i];
            if (je > -1) {
                FEElement* el = mesh.Element(je);

                // skip inactive elements. Needed for this to work with erosion
                if(!el->isActive()) continue;

                for (int k=0; k<el->GaussPoints(); ++k) {
                    FEMaterialPoint& mp = *(el->GetMaterialPoint(k));
                    double krnl = m_krnl->Kernel(pt, mp);
                    cst += pmat->CSTpt(mp)*krnl*mp.m_V0;
                    V += krnl*mp.m_V0;
                }
            }
        }
        if (V > 0) {
            cst /= V;
            if (V < 1) cst += pmat->CSTpt(pt)*(1 - V);
        }
    }
    
    return cst;
}

//-----------------------------------------------------------------------------
bool FENLAElement::Init()
{
    if (m_binit) return true;
    
    FEMesh& mesh = GetFEModel()->GetMesh();
    feLogInfo("Evaluating mesh topology...");
    if (m_topo.Create(&mesh) == false)
    {
        feLogError("Failed building mesh topo.");
        return false;
    }
    m_EPL.assign(mesh.Elements(), std::vector<int>());
    int NE = mesh.Elements();
    for (int i=0; i< mesh.Elements(); ++i) {
        std::vector<int> epl = std::vector<int>(1,i);
        m_EPL[i] = epl;
        double epls = epl.size();
    }
    m_binit = true;
    return true;
}

//-----------------------------------------------------------------------------
//! criterion average
double FENLAElement::DamageCriterionAverage(FEMaterialPoint& pt)
{
    double ca = 0;
    
    FEDamageCriterion* pmat = dynamic_cast<FEDamageCriterion*>(GetParent());
    assert(pmat);
    FEMesh& mesh = GetFEModel()->GetMesh();
    int ie = pt.m_elem->GetLocalID();
    FEElement* el = mesh.Element(ie);
    double V = 0;
    for (int k=0; k<el->GaussPoints(); ++k) {
        FEMaterialPoint& mp = *(el->GetMaterialPoint(k));
        ca += pmat->DCpt(mp)*mp.m_V0;
        V += mp.m_V0;
    }
    if (V > 0) ca /= V;
    
    return ca;
    
}

//-----------------------------------------------------------------------------
//! criterion stress tangent average
mat3ds FENLAElement::DamageCriterionTangentAverage(FEMaterialPoint& pt)
{
    mat3ds cst(0);
    
    FEDamageCriterion* pmat = dynamic_cast<FEDamageCriterion*>(GetParent());
    assert(pmat);
    FEMesh& mesh = GetFEModel()->GetMesh();
    int ie = pt.m_elem->GetLocalID();
    FEElement* el = mesh.Element(ie);
    double V = 0;
    for (int k=0; k<el->GaussPoints(); ++k) {
        FEMaterialPoint& mp = *(el->GetMaterialPoint(k));
        cst += pmat->CSTpt(mp)*mp.m_V0;
        V += mp.m_V0;
    }
    if (V > 0) cst /= V;
    
    return cst;
}

//-----------------------------------------------------------------------------
//! criterion average
void FENLAElement::PlasticityCriterionAverage(FEMaterialPoint& pt, FEMaterialPoint& rp)
{
    // get the reactive plasticity material
    FECoreBase* pm = GetParent();
    FEReactivePlasticity* pmat = dynamic_cast<FEReactivePlasticity*>(pm->GetParent());
    assert(pmat);
    FEReactivePlasticityMaterialPoint& ca = *rp.ExtractData<FEReactivePlasticityMaterialPoint>();
    int n = (int)pmat->m_pFlow->BondFamilies(pt);
    for (int l=0; l<n; ++l) {
        ca.m_w[l] = 0;
        ca.m_Fusi[l] = mat3d(0,0,0,0,0,0,0,0,0);
        ca.m_Fvsi[l] = mat3d(0,0,0,0,0,0,0,0,0);
        ca.m_Ku[l] = 0;
        ca.m_Kv[l] = 0;
        ca.m_gp[l] = 0;
        ca.m_gpp[l] = 0;
        ca.m_gc[l] = 0;
        ca.m_byld[l] = false;
    }
    ca.m_Fp = mat3d(0,0,0,0,0,0,0,0,0);
    ca.m_Rhat = 0;

    FEMesh& mesh = GetFEModel()->GetMesh();
    int ie = pt.m_elem->GetLocalID();
    FEElement* el = mesh.Element(ie);
    if (el->GetMatID() != pt.m_elem->GetMatID()) return;
    double V = 0;
    for (int k=0; k<el->GaussPoints(); ++k) {
        FEMaterialPoint& mp = *(el->GetMaterialPoint(k));
        pmat->ElasticDeformationGradient(mp);
        // extract plastic material point
        FEReactivePlasticityMaterialPoint& pp = *mp.ExtractData<FEReactivePlasticityMaterialPoint>();
        for (int l=0; l<n; ++l) {
            ca.m_w[l] += pp.m_w[l]*mp.m_V0;
            ca.m_Fusi[l] += pp.m_Fusi[l]*mp.m_V0;
            ca.m_Fvsi[l] += pp.m_Fvsi[l]*mp.m_V0;
            ca.m_Ku[l] += pp.m_Ku[l]*mp.m_V0;
            ca.m_Kv[l] += pp.m_Kv[l]*mp.m_V0;
            ca.m_gp[l] += pp.m_gp[l]*mp.m_V0;
            ca.m_gpp[l] += pp.m_gpp[l]*mp.m_V0;
            ca.m_gc[l] += pp.m_gc[l]*mp.m_V0;
            ca.m_byld[l] = ca.m_byld[l] & pp.m_byld[l];
        }
        ca.m_Fp += pp.m_Fp*mp.m_V0;
        ca.m_Rhat += pp.m_Rhat*mp.m_V0;
        V += mp.m_V0;
    }
    if (V > 0) {
        for (int l=0; l<n; ++l) {
            ca.m_w[l] /= V;
            ca.m_Fusi[l] /= V;
            ca.m_Fvsi[l] /= V;
            ca.m_Ku[l] /= V;
            ca.m_Kv[l] /= V;
            ca.m_gp[l] /= V;
            ca.m_gpp[l] /= V;
            ca.m_gc[l] /= V;
        }
        ca.m_Fp /= V;
        ca.m_Rhat /= V;
        pmat->OctahedralPlasticStrain(rp);
        pmat->ReactiveHeatSupplyDensity(rp);
    }
}
