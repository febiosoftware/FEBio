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
#include "FEChemicalReactionERD.h"
#include <FECore/FEElementTraits.h>
#include <FECore/DOFS.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include "FEElasticReactionDiffusionInterface.h"
#include <FEBioMix/FESolute.h>
#include <stdlib.h>


//=============================================================================
BEGIN_FECORE_CLASS(FEChemicalReactionERD, FEReactionERD)
ADD_PROPERTY(m_vRtmp, "vR", FEProperty::Optional)->SetLongName("Reactants");
ADD_PROPERTY(m_vPtmp, "vP", FEProperty::Optional)->SetLongName("Products");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEChemicalReactionERD::FEChemicalReactionERD(FEModel* pfem) : FEReactionERD(pfem)
{
    // additional initializations
    m_nsol = -1;

    m_pFwd = m_pRev = 0;
}

//-----------------------------------------------------------------------------
bool FEChemicalReactionERD::Init()
{
    // set the parents for the reaction rates
    if (m_pFwd) m_pFwd->m_pReact = this;
    if (m_pRev) m_pRev->m_pReact = this;

    // initialize base class
    if (FEReactionERD::Init() == false) return false;

    // initialize the reaction coefficients
    int isol;

    int nsol = m_psm->Solutes();

    // initialize the stoichiometric coefficients to zero
    m_nsol = nsol;
    m_vR.assign(nsol, 0);
    m_vP.assign(nsol, 0);
    m_v.assign(nsol, 0);

    // create the intmaps
    for (int i = 0; i < m_vRtmp.size(); ++i)
    {
        FEReactionSpeciesRefERD* pvr = m_vRtmp[i];
        SetStoichiometricCoefficient(m_solR, pvr->m_speciesID - 1, pvr->m_v);
    }
    for (int i = 0; i < m_vPtmp.size(); ++i)
    {
        FEReactionSpeciesRefERD* pvp = m_vPtmp[i];
        SetStoichiometricCoefficient(m_solP, pvp->m_speciesID - 1, pvp->m_v);
    }

    // cycle through all the solutes in the mixture and determine
    // if they participate in this reaction
    itrmap it;
    intmap solR = m_solR;
    intmap solP = m_solP;
    for (isol = 0; isol < nsol; ++isol) {
        int sid = m_psm->GetSolute(isol)->GetSoluteID() - 1;
        it = solR.find(sid);
        if (it != solR.end()) m_vR[isol] = it->second;
        it = solP.find(sid);
        if (it != solP.end()) m_vP[isol] = it->second;
    }

    // evaluate the net stoichiometric coefficient
    for (isol = 0; isol < nsol; ++isol) {
        m_v[isol] = m_vP[isol] - m_vR[isol];
    }

    return true;
}

//-----------------------------------------------------------------------------
//! Data serialization
void FEChemicalReactionERD::Serialize(DumpStream& ar)
{
    FEReactionERD::Serialize(ar);

    if (ar.IsShallow() == false)
    {
        if (ar.IsSaving())
        {
            itrmap p;
            ar << m_nsol << m_vR << m_vP << m_v;
            ar << (int)m_solR.size();
            for (p = m_solR.begin(); p != m_solR.end(); ++p) { ar << p->first; ar << p->second; }
            ar << (int)m_solP.size();
            for (p = m_solP.begin(); p != m_solP.end(); ++p) { ar << p->first; ar << p->second; }
        }
        else
        {
            // restore pointers
            if (m_pFwd) m_pFwd->m_pReact = this;
            if (m_pRev) m_pRev->m_pReact = this;

            ar >> m_nsol >> m_vR >> m_vP >> m_v;
            int size, id, vR;
            ar >> size;
            for (int i = 0; i < size; ++i);
            {
                ar >> id; ar >> vR;
                SetStoichiometricCoefficient(m_solR, id, vR);
            }
            ar >> size;
            for (int i = 0; i < size; ++i)
            {
                ar >> id; ar >> vR;
                SetStoichiometricCoefficient(m_solP, id, vR);
            }
        }
    }
}
