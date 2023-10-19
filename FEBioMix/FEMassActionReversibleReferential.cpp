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
#include "FEMassActionReversibleReferential.h"
#include "FESoluteInterface.h"
#include "FEBiphasic.h"

BEGIN_FECORE_CLASS(FEMassActionReversibleReferential, FEChemicalReaction)
// set material properties
ADD_PROPERTY(m_pFwd, "forward_rate", FEProperty::Optional);
ADD_PROPERTY(m_pRev, "reverse_rate", FEProperty::Optional);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMassActionReversibleReferential::FEMassActionReversibleReferential(FEModel* pfem) : FEChemicalReaction(pfem)
{
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionReversibleReferential::FwdReactionSupply(FEMaterialPoint& pt)
{
    FEElasticMaterialPoint* ep = pt.ExtractData<FEElasticMaterialPoint>();
    // get forward reaction rate
    double k = m_pFwd->ReactionRate(pt);

    // evaluate the reaction molar supply
    double zhat = k;

    // start with contribution from solutes
    int nsol = (int)m_psm->Solutes();
    for (int i = 0; i < nsol; ++i) {
        int vR = m_vR[i];
        if (vR > 0) {
            double c = m_psm->GetReferentialSoluteConcentration(pt, i);
            zhat *= pow(c, vR);
        }
    }

    // add contribution of solid-bound molecules
    const int nsbm = (int)m_psm->SBMs();
    for (int i = 0; i < nsbm; ++i) {
        int vR = m_vR[nsol + i];
        if (vR > 0) {
            double c = m_psm->SBMReferentialConcentration(pt, i);
            zhat *= pow(c, vR);
        }
    }

    return zhat / ep->m_J;
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionReversibleReferential::RevReactionSupply(FEMaterialPoint& pt)
{
    FEElasticMaterialPoint* ep = pt.ExtractData<FEElasticMaterialPoint>();
    // get forward reaction rate
    double k = m_pRev->ReactionRate(pt);

    // evaluate the reaction molar supply
    double zhat = k;

    // start with contribution from solutes
    int nsol = m_psm->Solutes();
    for (int i = 0; i < nsol; ++i) {
        int vP = m_vP[i];
        if (vP > 0) {
            double c = m_psm->GetReferentialSoluteConcentration(pt, i);
            zhat *= pow(c, vP);
        }
    }

    // add contribution of solid-bound molecules
    const int nsbm = m_psm->SBMs();
    for (int i = 0; i < nsbm; ++i) {
        int vP = m_vP[nsol + i];
        if (vP > 0) {
            double c = m_psm->SBMReferentialConcentration(pt, i);
            zhat *= pow(c, vP);
        }
    }

    return zhat / ep->m_J;
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionReversibleReferential::ReactionSupply(FEMaterialPoint& pt)
{
    double zhatF = FwdReactionSupply(pt);
    double zhatR = RevReactionSupply(pt);
    return zhatF - zhatR;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with strain at material point
mat3ds FEMassActionReversibleReferential::Tangent_ReactionSupply_Strain(FEMaterialPoint& pt)
{
    //double zhat = ReactionSupply(pt);
    //mat3dd I(1);
    //// evaluate the reaction molar supply
    //FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());
    //double phi0 = pbm->GetReferentialSolidVolumeFraction(pt);
    //FEElasticMaterialPoint& mp = *pt.ExtractData<FEElasticMaterialPoint>();
    //double J = mp.m_J;
    //mat3ds dzde = I * (-J / (J - phi0)) * zhat;
    //return dzde;
    return mat3ds(0.0);
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective pressure at material point
double FEMassActionReversibleReferential::Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt)
{
    //double dzhatFdp = 0;
    //double dzhatRdp = 0;
    //return dzhatFdp - dzhatRdp;
    return 0.0;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective concentration at material point
double FEMassActionReversibleReferential::Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol)
{
    //double dzhatFdc = 0;
    //double dzhatRdc = 0;
    //return dzhatFdc - dzhatRdc;
    return 0.0;
}