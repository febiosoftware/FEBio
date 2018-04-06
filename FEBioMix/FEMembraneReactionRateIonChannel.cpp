//
//  FEMembraneReactionRateIonChannel.cpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 4/6/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#include "FEMembraneReactionRateIonChannel.h"

// Material parameters for the FEMembraneReactionRateConst material
BEGIN_PARAMETER_LIST(FEMembraneReactionRateIonChannel, FEMembraneReactionRate)
ADD_PARAMETER2(m_g, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "g");
ADD_PARAMETER(m_sol, FE_PARAM_INT, "sol");
END_PARAMETER_LIST();

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
