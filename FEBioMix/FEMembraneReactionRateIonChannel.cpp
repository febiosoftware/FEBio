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
#include "FEMembraneReactionRateIonChannel.h"

// Material parameters for the FEMembraneReactionRateConst material
BEGIN_FECORE_CLASS(FEMembraneReactionRateIonChannel, FEMembraneReactionRate)
	ADD_PARAMETER(m_g, FE_RANGE_GREATER_OR_EQUAL(0.0), "g");
	ADD_PARAMETER(m_sol, "sol");
END_FECORE_CLASS();

bool FEMembraneReactionRateIonChannel::Init()
{
    // reset m_sol to be zero-based
    m_sol -= 1;
    
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
    
    double k = 0;
    if (ci != ce) k = R*T/pow(Fc*m_z,2)*m_g*log(ci/ce)/(ci-ce);
    else if (ce != 0) k = R*T/pow(Fc*m_z,2)*m_g/ce;
    
    return k;
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
    
    double dkdc = 0;
    if ((ci > 0) && (ce > 0)) dkdc = R*T/pow(Fc*m_z,2)*m_g*(ci*(1-log(ci/ce))-ce)/pow(ci-ce,2)/ci;
    
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
    
    double dkdc = 0;
    if ((ci > 0) && (ce > 0)) dkdc = R*T/pow(Fc*m_z,2)*m_g*(ce*(1+log(ci/ce))-ci)/pow(ci-ce,2)/ce;
    
    return dkdc;
}
