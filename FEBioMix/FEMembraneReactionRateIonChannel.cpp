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
#include "FEMembraneReactionRateIonChannel.h"
#include "FESoluteInterface.h"
#include "FESolutesMaterialPoint.h"
#include "FESolute.h"
#include <FEBioMech/FEElasticMaterialPoint.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>

// Material parameters for the FEMembraneReactionRateConst material
BEGIN_FECORE_CLASS(FEMembraneReactionRateIonChannel, FEMembraneReactionRate)
	ADD_PARAMETER(m_g, FE_RANGE_GREATER_OR_EQUAL(0.0), "g");
	ADD_PARAMETER(m_sol, "sol");
    ADD_PARAMETER(m_sbm, "sbm");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMembraneReactionRateIonChannel::FEMembraneReactionRateIonChannel(FEModel* pfem) : FEMembraneReactionRate(pfem)
{
    m_sol = -1;
    m_sbm = -1;
    m_lid = -1;
    m_g = 0;
    m_z = 0;
}

//-----------------------------------------------------------------------------
bool FEMembraneReactionRateIonChannel::Init()
{
    if (FEMembraneReactionRate::Init() == false) return false;
    
    // do only once
    if (m_lid == -1) {
        // get number of DOFS
        DOFS& fedofs = GetFEModel()->GetDOFS();
        int MAX_CDOFS = fedofs.GetVariableSize("concentration");
        // check validity of sol
        if (m_sol < 1 || m_sol > MAX_CDOFS) {
            feLogError("sol value outside of valid range for solutes");
            return false;
        }
        
        FEModel& fem = *GetFEModel();
        int N = GetFEModel()->GlobalDataItems();
        for (int i=0; i<N; ++i)
        {
            FESoluteData* psd = dynamic_cast<FESoluteData*>(fem.GetGlobalData(i));
            if (psd && (psd->GetID() == m_sol)) {
                m_lid = m_sol - 1;
                m_z = psd->m_z;
                break;
            }
        }
    }
    
    return true;
}

//-----------------------------------------------------------------------------
double FEMembraneReactionRateIonChannel::ReactionRate(FEMaterialPoint& pt)
{
    FESolutesMaterialPoint& ps = *pt.ExtractData<FESolutesMaterialPoint>();
    double ci = ps.m_ci[m_lid];
    double ce = ps.m_ce[m_lid];
    double R = GetGlobalConstant("R");
    double T = GetGlobalConstant("T");
    double Fc = GetGlobalConstant("Fc");
    FESoluteInterface* psi = m_pReact->m_psm; assert(psi);
    double ksi = (m_sbm > -1) ? psi->SBMArealConcentration(pt, m_sbm - 1): 1.0;

    double k = 0;
    if ((ci > 0) && (ce > 0))
        k = (ci != ce) ? R*T*m_g/ksi/pow(Fc*m_z,2)*log(ci/ce)/(ci-ce) : R*T*m_g/pow(Fc*m_z,2)/ce;
    
    return k;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
double FEMembraneReactionRateIonChannel::Tangent_ReactionRate_Strain(FEMaterialPoint& pt)
{
    // get the areal strain
    FEElasticMaterialPoint& pe = *(pt.ExtractData<FEElasticMaterialPoint>());
    FEShellElement*sel = dynamic_cast<FEShellElement*>(pt.m_elem);
    assert(sel);
    double Jg = pe.m_J*sel->Evaluate(sel->m_h0, pt.m_index)/sel->Evaluate(sel->m_ht, pt.m_index);
    return ReactionRate(pt)/Jg;
}

//-----------------------------------------------------------------------------
double FEMembraneReactionRateIonChannel::Tangent_ReactionRate_Ci(FEMaterialPoint& pt, const int isol)
{
    if (isol != m_lid)  return 0;
    
    FESolutesMaterialPoint& ps = *(pt.ExtractData<FESolutesMaterialPoint>());
    double ci = ps.m_ci[m_lid];
    double ce = ps.m_ce[m_lid];
    double R = GetGlobalConstant("R");
    double T = GetGlobalConstant("T");
    double Fc = GetGlobalConstant("Fc");
    FESoluteInterface* psi = m_pReact->m_psm; assert(psi);
    double ksi = (m_sbm > -1) ? psi->SBMArealConcentration(pt, m_sbm - 1) : 1.0;

    double dkdc = 0;
    if ((ci > 0) && (ce > 0))
        dkdc = (ci != ce) ? R*T/pow(Fc*m_z,2)*m_g/ksi*(ci*(1-log(ci/ce))-ce)/pow(ci-ce,2)/ci : -R*T*m_g/pow(ci*Fc*m_z,2)/2;
    
    return dkdc;
}

//-----------------------------------------------------------------------------
double FEMembraneReactionRateIonChannel::Tangent_ReactionRate_Ce(FEMaterialPoint& pt, const int isol)
{
    if (isol != m_lid)  return 0;
    
    FESolutesMaterialPoint& ps = *(pt.ExtractData<FESolutesMaterialPoint>());
    double ci = ps.m_ci[m_lid];
    double ce = ps.m_ce[m_lid];
    double R = GetGlobalConstant("R");
    double T = GetGlobalConstant("T");
    double Fc = GetGlobalConstant("Fc");
    FESoluteInterface* psi = m_pReact->m_psm; assert(psi);
    double ksi = (m_sbm > -1) ? psi->SBMArealConcentration(pt, m_sbm - 1) : 1.0;

    double dkdc = 0;
    if ((ci > 0) && (ce > 0))
        dkdc = (ci != ce) ? R*T*m_g/ksi/pow(Fc*m_z,2)*(ce*(1+log(ci/ce))-ci)/pow(ci-ce,2)/ce : -R*T*m_g/pow(ci*Fc*m_z,2)/2;
    
    return dkdc;
}
