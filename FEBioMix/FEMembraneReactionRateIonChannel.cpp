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

// Material parameters for the FEMembraneReactionRateConst material
BEGIN_FECORE_CLASS(FEMembraneReactionRateIonChannel, FEMembraneReactionRate)
	ADD_PARAMETER(m_g, FE_RANGE_GREATER_OR_EQUAL(0.0), "g");
	ADD_PARAMETER(m_sol, "sol");
    ADD_PARAMETER(m_sbm, "sbm");
END_FECORE_CLASS();

bool FEMembraneReactionRateIonChannel::Init()
{
    // reset m_sol and m_sbm to be zero-based
    m_sol -= 1;
    m_sbm -= 1;
    
    // membrane reaction rate is child of membrane reaction
    FEMembraneReaction* m_MRp = dynamic_cast<FEMembraneReaction*>(GetParent());
    if (m_MRp == nullptr) return false;
    m_z = (m_MRp->FindSoluteData(m_sol))->m_z;
    return true;
}

double FEMembraneReactionRateIonChannel::ReactionRate(FEMaterialPoint& pt)
{
    FESolutesMaterialPoint& ps = *(pt.ExtractData<FESolutesMaterialPoint>());
    double ci = ps.m_ci[m_sol];
    double ce = ps.m_ce[m_sol];
    double R = GetFEModel()->GetGlobalConstant("R");
    double T = GetFEModel()->GetGlobalConstant("T");
    double Fc = GetFEModel()->GetGlobalConstant("Fc");
    FEMultiphasic* pMP = m_pReact->m_pMP;
    assert(pMP);
    double ksi = pMP->SBMArealConcentration(pt, m_sbm);
    
    double k = 0;
    if ((ci > 0) && (ce > 0))
        k = (ci != ce) ? R*T*m_g/ksi/pow(Fc*m_z,2)*log(ci/ce)/(ci-ce) : R*T*m_g/pow(Fc*m_z,2)/ce;
    
    return k;
}

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

double FEMembraneReactionRateIonChannel::Tangent_ReactionRate_Ci(FEMaterialPoint& pt, const int isol)
{
    if (isol != m_sol)  return 0;
    
    FESolutesMaterialPoint& ps = *(pt.ExtractData<FESolutesMaterialPoint>());
    double ci = ps.m_ci[m_sol];
    double ce = ps.m_ce[m_sol];
    double R = GetFEModel()->GetGlobalConstant("R");
    double T = GetFEModel()->GetGlobalConstant("T");
    double Fc = GetFEModel()->GetGlobalConstant("Fc");
    FEMultiphasic* pMP = m_pReact->m_pMP;
    assert(pMP);
    double ksi = pMP->SBMArealConcentration(pt, m_sbm);

    double dkdc = 0;
    if ((ci > 0) && (ce > 0))
        dkdc = (ci != ce) ? R*T/pow(Fc*m_z,2)*m_g/ksi*(ci*(1-log(ci/ce))-ce)/pow(ci-ce,2)/ci : -R*T*m_g/pow(ci*Fc*m_z,2)/2;
    
    return dkdc;
}

double FEMembraneReactionRateIonChannel::Tangent_ReactionRate_Ce(FEMaterialPoint& pt, const int isol)
{
    if (isol != m_sol)  return 0;
    
    FESolutesMaterialPoint& ps = *(pt.ExtractData<FESolutesMaterialPoint>());
    double ci = ps.m_ci[m_sol];
    double ce = ps.m_ce[m_sol];
    double R = GetFEModel()->GetGlobalConstant("R");
    double T = GetFEModel()->GetGlobalConstant("T");
    double Fc = GetFEModel()->GetGlobalConstant("Fc");
    FEMultiphasic* pMP = m_pReact->m_pMP;
    assert(pMP);
    double ksi = pMP->SBMArealConcentration(pt, m_sbm);

    double dkdc = 0;
    if ((ci > 0) && (ce > 0))
        dkdc = (ci != ce) ? R*T*m_g/ksi/pow(Fc*m_z,2)*(ce*(1+log(ci/ce))-ci)/pow(ci-ce,2)/ce : -R*T*m_g/pow(ci*Fc*m_z,2)/2;
    
    return dkdc;
}
