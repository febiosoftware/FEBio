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

double FEMembraneReactionRateIonChannel::ReactionRate(FEMaterialPoint& pt)
{
    FESolutesMaterialPoint& ps = *(pt.ExtractData<FESolutesMaterialPoint>());
    double ci = ps.m_ci[m_sol -1];
    double R = GetFEModel()->GetGlobalConstant("R");
    double T = GetFEModel()->GetGlobalConstant("T");
    double Fc = GetFEModel()->GetGlobalConstant("Fc");
    
    double k = 0;
    if (ci > 0) k = R*T/(Fc*Fc)*m_g/ci;
    
    return k;
}

double FEMembraneReactionRateIonChannel::Tangent_ReactionRate_Ci(FEMaterialPoint& pt, const int isol)
{
    FESolutesMaterialPoint& ps = *(pt.ExtractData<FESolutesMaterialPoint>());
    double ci = ps.m_ci[m_sol -1];
    double R = GetFEModel()->GetGlobalConstant("R");
    double T = GetFEModel()->GetGlobalConstant("T");
    double Fc = GetFEModel()->GetGlobalConstant("Fc");
    
    double dkdc = 0;
    if (ci > 0) dkdc = -R*T/(Fc*Fc)*m_g/(ci*ci);
    
    return dkdc;
}
