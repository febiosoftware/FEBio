//
//  FEMembraneReactionRateVoltageGated.cpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 4/30/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#include "FEMembraneReactionRateVoltageGated.h"

// Material parameters for the FEMembraneReactionRateVoltageGated material
BEGIN_FECORE_CLASS(FEMembraneReactionRateVoltageGated, FEMembraneReactionRate)
	ADD_PARAMETER(m_a, FE_RANGE_GREATER_OR_EQUAL(0.0), "a");
	ADD_PARAMETER(m_b, FE_RANGE_GREATER_OR_EQUAL(0.0), "b");
	ADD_PARAMETER(m_c, FE_RANGE_GREATER_OR_EQUAL(0.0), "c");
	ADD_PARAMETER(m_d, FE_RANGE_GREATER_OR_EQUAL(0.0), "d");
	ADD_PARAMETER(m_sol, "sol");
END_FECORE_CLASS();

bool FEMembraneReactionRateVoltageGated::Init()
{
    // reset m_sol to be zero-based
    m_sol -= 1;
    
    // membrane reaction rate is child of membrane reaction
    FEMembraneReaction* m_MRp = dynamic_cast<FEMembraneReaction*>(GetParent());
    if (m_MRp == nullptr) return false;
    m_z = (m_MRp->FindSoluteData(m_sol))->m_z;
    return true;
}

double FEMembraneReactionRateVoltageGated::ReactionRate(FEMaterialPoint& pt)
{
    FESolutesMaterialPoint& ps = *(pt.ExtractData<FESolutesMaterialPoint>());
    double ci = ps.m_ci[m_sol];
    double ce = ps.m_ce[m_sol];
    if (ce == 0) return 0;
    double x = ci/ce;
    
    double k = (m_c*log(x) + m_d)/(m_a + pow(x,m_b));
    
    return k;
}

double FEMembraneReactionRateVoltageGated::Tangent_ReactionRate_Ci(FEMaterialPoint& pt, const int isol)
{
    if (isol != m_sol)  return 0;
    
    FESolutesMaterialPoint& ps = *(pt.ExtractData<FESolutesMaterialPoint>());
    double ci = ps.m_ci[m_sol];
    double ce = ps.m_ce[m_sol];
    if ((ce == 0) || (ci == 0)) return 0;
    double x = ci/ce;
    double xb = pow(x,m_b);

    double dkdc = (m_a*m_c + xb*(m_c-m_b*m_d) - m_b*m_c*log(x))/pow(m_a+xb,2)/ci;
    
    return dkdc;
}

double FEMembraneReactionRateVoltageGated::Tangent_ReactionRate_Ce(FEMaterialPoint& pt, const int isol)
{
    if (isol != m_sol)  return 0;
    
    FESolutesMaterialPoint& ps = *(pt.ExtractData<FESolutesMaterialPoint>());
    double ci = ps.m_ci[m_sol];
    double ce = ps.m_ce[m_sol];
    if (ce == 0) return 0;
    double x = ci/ce;
    double xb = pow(x,m_b);

    double dkdc = -(m_a*m_c + xb*(m_c-m_b*m_d) - m_b*m_c*log(x))/pow(m_a+xb,2)/ce;

    
    return dkdc;
}
