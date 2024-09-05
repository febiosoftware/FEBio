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
#include "FEMassActionForwardERD.h"
#include "FEElasticReactionDiffusionInterface.h"
#include <FEBioMix/FEBiphasic.h>

BEGIN_FECORE_CLASS(FEMassActionForwardERD, FEChemicalReactionERD)
END_FECORE_CLASS();

//! constructor
//-----------------------------------------------------------------------------
FEMassActionForwardERD::FEMassActionForwardERD(FEModel* pfem) : FEChemicalReactionERD(pfem)
{
    // set material properties
    ADD_PROPERTY(m_pFwd, "forward_rate", FEProperty::Optional);
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionForwardERD::ReactionSupply(FEMaterialPoint& pt)
{
    // get reaction rate
    double kF = m_pFwd->ReactionRate(pt);

    // evaluate the reaction molar supply
    double zhat = kF;

    // start with contribution from solutes
    int nsol = m_psm->Solutes();
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
//! tangent of molar supply with strain at material point
mat3ds FEMassActionForwardERD::Tangent_ReactionSupply_Strain(FEMaterialPoint& pt)
{
    // if the reaction supply is insensitive to strain
    return mat3ds(0);


    //for now just return 0
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective concentration at material point
double FEMassActionForwardERD::Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol)
{
    const int nsol = m_nsol;

    double zhat = ReactionSupply(pt);
    double dzhatdc = 0.0;
    double c = m_psm->GetActualSoluteConcentration(pt, sol);
    if ((zhat > 0.0) && (c > 0.0))
        dzhatdc = (m_vR[sol] / c) * zhat;

    return dzhatdc;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with Cauchy stress (sigma) at material point
mat3ds FEMassActionForwardERD::Tangent_ReactionSupply_Stress(FEMaterialPoint& pt)
{
    return m_pFwd->Tangent_ReactionRate_Stress(pt);
    //return mat3ds(0);
}