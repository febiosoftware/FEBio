/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
