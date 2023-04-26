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
#include "FEMassActionForwardReferential.h"
#include "FESoluteInterface.h"
#include "FEBiphasic.h"
#include <FEBioMech/FEKinematicGrowth.h>

BEGIN_FECORE_CLASS(FEMassActionForwardReferential, FEChemicalReaction)
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMassActionForwardReferential::FEMassActionForwardReferential(FEModel* pfem) : FEChemicalReaction(pfem)
{
    // set material properties
    ADD_PROPERTY(m_pFwd, "forward_rate", FEProperty::Optional);
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionForwardReferential::ReactionSupply(FEMaterialPoint& pt)
{
    // get reaction rate
    double kF = m_pFwd->ReactionRate(pt);

    double zhat = kF;

    // evaluate the reaction molar supply
    FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());
    double phi0 = pbm->GetReferentialSolidVolumeFraction(pt);
    FEKinematicMaterialPoint& kp = *pt.ExtractData<FEKinematicMaterialPoint>();

    // start with contribution from solutes
    int nsol = m_psm->Solutes();
    for (int i = 0; i < nsol; ++i) {
        int vR = m_vR[i];
        if (vR > 0) {
            double c = m_psm->GetReferentialSoluteConcentration(pt, i);
            zhat *= pow(c, vR);
        }
    }

    // add contribution of solid-bound molecules
    const int nsbm = m_psm->SBMs();
    for (int i = 0; i < nsbm; ++i) {
        int vR = m_vR[nsol + i];
        if (vR > 0) {
            double c = m_psm->SBMConcentration(pt, i);
            zhat *= pow(c, vR);
        }
    }

    return zhat;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with strain at material point
mat3ds FEMassActionForwardReferential::Tangent_ReactionSupply_Strain(FEMaterialPoint& pt)
{
    const int nsol = m_nsol;
    const int nsbm = (int)m_v.size() - nsol;

    FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());
    FEElasticMaterialPoint& ept = *pt.ExtractData<FEElasticMaterialPoint>();

    double kF = m_pFwd->ReactionRate(pt);
    //mat3ds dkFde = m_pFwd->Tangent_ReactionRate_Strain(pt);
    double zhat = ReactionSupply(pt);
    mat3ds dzhatde = mat3dd(0);
    mat3ds I = mat3dd(1);

    for (int isbm = 0; isbm < nsbm; ++isbm)
        dzhatde -= I * (m_vR[nsol + isbm]);

    dzhatde *= zhat;
    mat3ds dzde = -1.0 * zhat * I;
    return dzhatde;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective pressure at material point
double FEMassActionForwardReferential::Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt)
{
    double dzhatdp = 0;
    return dzhatdp;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective concentration at material point
double FEMassActionForwardReferential::Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol)
{
    double dzhatdc = 0;
    return dzhatdc;
}
