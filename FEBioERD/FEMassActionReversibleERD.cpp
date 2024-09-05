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
#include "FEMassActionReversibleERD.h"
#include "FEElasticReactionDiffusionInterface.h"
#include <FEBioMix/FEBiphasic.h>

BEGIN_FECORE_CLASS(FEMassActionReversibleERD, FEChemicalReactionERD)
// set material properties
ADD_PROPERTY(m_pFwd, "forward_rate", FEProperty::Optional);
ADD_PROPERTY(m_pRev, "reverse_rate", FEProperty::Optional);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMassActionReversibleERD::FEMassActionReversibleERD(FEModel* pfem) : FEChemicalReactionERD(pfem)
{

}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionReversibleERD::FwdReactionSupply(FEMaterialPoint& pt)
{
    // get forward reaction rate
    double k = m_pFwd->ReactionRate(pt);

    // evaluate the reaction molar supply
    double zhat = k;

    // start with contribution from solutes
    int nsol = 0;
    nsol = m_psm->Solutes();
    for (int i = 0; i < nsol; ++i) {
        int vR = m_vR[i];
        if (vR > 0) {
            double c = m_psm->GetActualSoluteConcentration(pt, i);
            zhat *= pow(c, vR);
        }
    }

    return zhat;
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionReversibleERD::RevReactionSupply(FEMaterialPoint& pt)
{
    // get forward reaction rate
    double k = m_pRev->ReactionRate(pt);

    // evaluate the reaction molar supply
    double zhat = k;

    // start with contribution from solutes
    int nsol = m_psm->Solutes();
    for (int i = 0; i < nsol; ++i) {
        int vP = m_vP[i];
        if (vP > 0) {
            double c = m_psm->GetActualSoluteConcentration(pt, i);
            zhat *= pow(c, vP);
        }
    }

    return zhat;
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionReversibleERD::ReactionSupply(FEMaterialPoint& pt)
{
    double zhatF = FwdReactionSupply(pt);
    double zhatR = RevReactionSupply(pt);
    return zhatF - zhatR;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with strain at material point
mat3ds FEMassActionReversibleERD::Tangent_ReactionSupply_Strain(FEMaterialPoint& pt)
{
    // if the reaction supply is insensitive to strain
    return mat3ds(0);

}

//-----------------------------------------------------------------------------
//! tangent of molar supply with  concentration at material point
double FEMassActionReversibleERD::Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol)
{
    // if the reaction supply is insensitive to concentration

    const int nsol = m_nsol;

    // if the derivative is taken with respect to a solid-bound molecule, return 0
    if (sol >= nsol) {
        return 0;
    }

    // forward reaction
    double zhatF = FwdReactionSupply(pt);
    double dzhatFdc = 0.0;
    double c = m_psm->GetActualSoluteConcentration(pt, sol);
    if ((zhatF > 0.0) && (c > 0.0))
        dzhatFdc = (m_vR[sol] / c) * zhatF;

    // reverse reaction
    double zhatR = RevReactionSupply(pt);
    double dzhatRdc = 0.0;
    if ((zhatR > 0.0) && (c > 0.0))
        dzhatRdc = (m_vP[sol] / c) * zhatR;

    return dzhatFdc - dzhatRdc;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with Cauchy stress (sigma) at material point
mat3ds FEMassActionReversibleERD::Tangent_ReactionSupply_Stress(FEMaterialPoint& pt)
{
    // if the reaction supply is insensitive to strain
    return mat3ds(0);

}