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
#include "FEMembraneReactionRateVoltageGated.h"
#include "FESoluteInterface.h"
#include "FESolutesMaterialPoint.h"
#include "FESolute.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>

// Material parameters for the FEMembraneReactionRateVoltageGated material
BEGIN_FECORE_CLASS(FEMembraneReactionRateVoltageGated, FEMembraneReactionRate)
	ADD_PARAMETER(m_a, FE_RANGE_GREATER_OR_EQUAL(0.0), "a");
	ADD_PARAMETER(m_b, FE_RANGE_GREATER_OR_EQUAL(0.0), "b");
	ADD_PARAMETER(m_c, FE_RANGE_GREATER_OR_EQUAL(0.0), "c");
	ADD_PARAMETER(m_d, FE_RANGE_GREATER_OR_EQUAL(0.0), "d");
	ADD_PARAMETER(m_sol, "sol");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMembraneReactionRateVoltageGated::FEMembraneReactionRateVoltageGated(FEModel* pfem) : FEMembraneReactionRate(pfem)
{
    m_a = m_b = m_c = m_d = 0;
    m_sol = m_lid = -1;
    m_z = 0;
}

//-----------------------------------------------------------------------------
bool FEMembraneReactionRateVoltageGated::Init()
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
double FEMembraneReactionRateVoltageGated::ReactionRate(FEMaterialPoint& pt)
{
    FESolutesMaterialPoint& ps = *(pt.ExtractData<FESolutesMaterialPoint>());
    double ci = ps.m_ci[m_lid];
    double ce = ps.m_ce[m_lid];
    if (ce == 0) return 0;
    double x = ci/ce;
    
    double k = (m_c*log(x) + m_d)/(m_a + pow(x,m_b));
    
    return k;
}

//-----------------------------------------------------------------------------
double FEMembraneReactionRateVoltageGated::Tangent_ReactionRate_Ci(FEMaterialPoint& pt, const int isol)
{
    if (isol != m_lid)  return 0;
    
    FESolutesMaterialPoint& ps = *(pt.ExtractData<FESolutesMaterialPoint>());
    double ci = ps.m_ci[m_lid];
    double ce = ps.m_ce[m_lid];
    if ((ce == 0) || (ci == 0)) return 0;
    double x = ci/ce;
    double xb = pow(x,m_b);

    double dkdc = (m_a*m_c + xb*(m_c-m_b*m_d) - m_b*m_c*log(x))/pow(m_a+xb,2)/ci;
    
    return dkdc;
}

//-----------------------------------------------------------------------------
double FEMembraneReactionRateVoltageGated::Tangent_ReactionRate_Ce(FEMaterialPoint& pt, const int isol)
{
    if (isol != m_lid)  return 0;
    
    FESolutesMaterialPoint& ps = *(pt.ExtractData<FESolutesMaterialPoint>());
    double ci = ps.m_ci[m_lid];
    double ce = ps.m_ce[m_lid];
    if (ce == 0) return 0;
    double x = ci/ce;
    double xb = pow(x,m_b);

    double dkdc = -(m_a*m_c + xb*(m_c-m_b*m_d) - m_b*m_c*log(x))/pow(m_a+xb,2)/ce;

    
    return dkdc;
}
