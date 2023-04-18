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
#include "FEMembraneReaction.h"
#include <FECore/FEElementTraits.h>
#include <FECore/DOFS.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include "FESoluteInterface.h"
#include "FESolute.h"
#include <stdlib.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEInternalReactantSpeciesRef, FEReactionSpeciesRef)
    ADD_PARAMETER(m_v, "vRi");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEInternalProductSpeciesRef, FEReactionSpeciesRef)
    ADD_PARAMETER(m_v, "vPi");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEExternalReactantSpeciesRef, FEReactionSpeciesRef)
    ADD_PARAMETER(m_v, "vRe");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEExternalProductSpeciesRef, FEReactionSpeciesRef)
    ADD_PARAMETER(m_v, "vPe");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEMembraneReaction, FEReaction)
    ADD_PARAMETER(m_Vovr, "override_vbar")->SetFlags(FE_PARAM_WATCH);
    ADD_PARAMETER(m_Vbar , "Vbar")->SetWatchVariable(&m_Vovr);

    ADD_PROPERTY(m_vRtmp, "vR", FEProperty::Optional)->SetLongName("Membrane reactants");
    ADD_PROPERTY(m_vPtmp, "vP", FEProperty::Optional)->SetLongName("Membrane products");
    ADD_PROPERTY(m_vRitmp, "vRi", FEProperty::Optional)->SetLongName("Inner membrane reactants");
    ADD_PROPERTY(m_vPitmp, "vPi", FEProperty::Optional)->SetLongName("Inner membrane products");;
    ADD_PROPERTY(m_vRetmp, "vRe", FEProperty::Optional)->SetLongName("Outer membrane reactants");
    ADD_PROPERTY(m_vPetmp, "vPe", FEProperty::Optional)->SetLongName("Outer membrane products");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMembraneReaction::FEMembraneReaction(FEModel* pfem) : FEReaction(pfem)
{
    // additional initializations
    m_Vbar = 0.0;
    m_Vovr = false;
    m_nsol = -1;

	m_pFwd = m_pRev = 0;
}

//-----------------------------------------------------------------------------
// Finds the solute id of a solute with given ID nsol.
// This currently returns either nsol if a solute was found or -1 if not
FESoluteData* FEMembraneReaction::GetSolute(int nsol)
{
    FEModel& fem = *GetFEModel();
    int N = fem.GlobalDataItems();
    for (int i=0; i<N; ++i)
    {
        FESoluteData* psd = dynamic_cast<FESoluteData*>(fem.GetGlobalData(i));
        if (psd)
        {
            if (psd->GetID()-1 == nsol) return psd;
        }
    }
    return nullptr;
}

//-----------------------------------------------------------------------------
bool FEMembraneReaction::Init()
{
    // set the parents for the reaction rates
    if (m_pFwd) m_pFwd->m_pReact = this;
    if (m_pRev) m_pRev->m_pReact = this;
    
    // initialize base class
    if (FEReaction::Init() == false) return false;
    
    //************* reactants and products in multiphasic domain **************

    // create the intmaps
    for (int i = 0; i < m_vRtmp.size(); ++i)
    {
        FEReactionSpeciesRef* pvr = m_vRtmp[i];
        if (pvr->IsSolute()) SetStoichiometricCoefficient(m_solR, pvr->m_speciesID - 1, pvr->m_v);
        if (pvr->IsSBM()   ) SetStoichiometricCoefficient(m_sbmR, pvr->m_speciesID - 1, pvr->m_v);
    }
    for (int i = 0; i < m_vPtmp.size(); ++i)
    {
        FEReactionSpeciesRef* pvp = m_vPtmp[i];
        if (pvp->IsSolute()) SetStoichiometricCoefficient(m_solP, pvp->m_speciesID - 1, pvp->m_v);
        if (pvp->IsSBM()   ) SetStoichiometricCoefficient(m_sbmP, pvp->m_speciesID - 1, pvp->m_v);
    }
    for (int i = 0; i < m_vRitmp.size(); ++i)
    {
        FEReactionSpeciesRef* pvr = m_vRitmp[i];
        assert(pvr->IsSolute());
        if (pvr->IsSolute()) SetStoichiometricCoefficient(m_solRi, pvr->m_speciesID - 1, pvr->m_v);
    }
    for (int i = 0; i < m_vPitmp.size(); ++i)
    {
        FEReactionSpeciesRef* pvp = m_vPitmp[i];
        assert(pvp->IsSolute());
        if (pvp->IsSolute()) SetStoichiometricCoefficient(m_solPi, pvp->m_speciesID - 1, pvp->m_v);
    }
    for (int i = 0; i < m_vRetmp.size(); ++i)
    {
        FEReactionSpeciesRef* pvr = m_vRetmp[i];
        assert(pvr->IsSolute());
        if (pvr->IsSolute()) SetStoichiometricCoefficient(m_solRe, pvr->m_speciesID - 1, pvr->m_v);
    }
    for (int i = 0; i < m_vPetmp.size(); ++i)
    {
        FEReactionSpeciesRef* pvp = m_vPetmp[i];
        assert(pvp->IsSolute());
        if (pvp->IsSolute()) SetStoichiometricCoefficient(m_solPe, pvp->m_speciesID - 1, pvp->m_v);
    }

    // initialize the reaction coefficients
    const int nsol = m_psm->Solutes();
    const int nsbm = m_psm->SBMs();
    const int ntot = nsol + nsbm;
    
    // initialize the stoichiometric coefficients to zero
    m_nsol = nsol;
    m_vR.assign(ntot, 0);
    m_vP.assign(ntot, 0);
    m_v.assign(ntot, 0);
    
    // cycle through all the solutes in the mixture and determine
    // if they participate in this reaction
    itrmap it;
    intmap solR = m_solR;
    intmap solP = m_solP;
    for (int isol = 0; isol<nsol; ++isol) {
        int sid = m_psm->GetSolute(isol)->GetSoluteID() - 1;
        it = solR.find(sid);
        if (it != solR.end()) m_vR[isol] = it->second;
        it = solP.find(sid);
        if (it != solP.end()) m_vP[isol] = it->second;
    }
    
    // cycle through all the solid-bound molecules in the mixture
    // and determine if they participate in this reaction
    intmap sbmR = m_sbmR;
    intmap sbmP = m_sbmP;
    for (int isbm = 0; isbm<nsbm; ++isbm) {
        int sid = m_psm->GetSBM(isbm)->GetSBMID() - 1;
        it = sbmR.find(sid);
        if (it != sbmR.end()) m_vR[nsol + isbm] = it->second;
        it = sbmP.find(sid);
        if (it != sbmP.end()) m_vP[nsol + isbm] = it->second;
    }
    
    // evaluate the net stoichiometric coefficient
    for (int itot = 0; itot<ntot; ++itot) {
        m_v[itot] = m_vP[itot] - m_vR[itot];
    }
    
    //********* reactants and products on either side of membrane **********
    // count total number of solutes in model
    DOFS& fedofs = GetFEModel()->GetDOFS();
    int MAX_DDOFS = fedofs.GetVariableSize("shell concentration");
    m_NSOL = MAX_DDOFS;
    m_z.assign(MAX_DDOFS, 0);

    // initialize the stoichiometric coefficients to zero
    m_vRi.assign(MAX_DDOFS, 0); m_vRe.assign(MAX_DDOFS, 0);
    m_vPi.assign(MAX_DDOFS, 0); m_vPe.assign(MAX_DDOFS, 0);
    m_vi.assign(MAX_DDOFS, 0); m_ve.assign(MAX_DDOFS, 0);

    // cycle through all the solutes in the mixture and determine
    // if they participate in this reaction
    for (int ISOL = 0; ISOL<MAX_DDOFS; ++ISOL) {
        it = m_solRi.find(ISOL);
        if (it != m_solRi.end()) m_vRi[ISOL] = it->second;
        it = m_solPi.find(ISOL);
        if (it != m_solPi.end()) m_vPi[ISOL] = it->second;
        it = m_solRe.find(ISOL);
        if (it != m_solRe.end()) m_vRe[ISOL] = it->second;
        it = m_solPe.find(ISOL);
        if (it != m_solPe.end()) m_vPe[ISOL] = it->second;
    }

    // evaluate the net stoichiometric coefficient
    for (int ISOL = 0; ISOL<MAX_DDOFS; ++ISOL) {
        m_vi[ISOL] = m_vPi[ISOL] - m_vRi[ISOL];
        m_ve[ISOL] = m_vPe[ISOL] - m_vRe[ISOL];
        FESoluteData* sd = GetSolute(ISOL);
        m_z[ISOL] = sd->m_z;
    }

    //************** continue with all reactants and products ***************
    
    // evaluate the weighted molar volume of reactants and products
    if (!m_Vovr) {
        m_Vbar = 0;
        for (int isol = 0; isol<nsol; ++isol)
            m_Vbar += m_v[isol] * m_psm->GetSolute(isol)->MolarMass() / m_psm->GetSolute(isol)->Density();
        for (int isbm = 0; isbm<nsbm; ++isbm)
            m_Vbar += m_v[nsol + isbm] * m_psm->GetSBM(isbm)->MolarMass() / m_psm->GetSBM(isbm)->Density();
        for (int ISOL = 0; ISOL<MAX_DDOFS; ++ISOL) {
            FESoluteData* sd = GetSolute(ISOL);
            m_Vbar += m_vi[ISOL] * sd->m_M / sd->m_rhoT;
            m_Vbar += m_ve[ISOL] * sd->m_M / sd->m_rhoT;
        }
    }
    
    // check that the reaction satisfies electroneutrality
    int znet = 0;
    for (int isol = 0; isol<nsol; ++isol)
        znet += m_v[isol] * m_psm->GetSolute(isol)->ChargeNumber();
    for (int isbm = 0; isbm<nsbm; ++isbm)
        znet += m_v[nsol + isbm] * m_psm->GetSBM(isbm)->ChargeNumber();
    for (int ISOL = 0; ISOL<MAX_DDOFS; ++ISOL) {
        znet += m_vi[ISOL] * m_z[ISOL];
        znet += m_ve[ISOL] * m_z[ISOL];
    }
	if (znet != 0) {
		feLogError("membrane reaction must satisfy electroneutrality");
		return false;
	}
    
    return true;
}

//-----------------------------------------------------------------------------
//! Data serialization
void FEMembraneReaction::Serialize(DumpStream& ar)
{
    FEReaction::Serialize(ar);
    
    if (ar.IsShallow() == false)
    {
        if (ar.IsSaving())
        {
            itrmap p;
            ar << m_nsol << m_vR << m_vP << m_v << m_Vovr << m_NSOL;
            ar << (int) m_solR.size();
            for (p = m_solR.begin(); p!=m_solR.end(); ++p) {ar << p->first; ar << p->second;}
            ar << (int) m_solP.size();
            for (p = m_solP.begin(); p!=m_solP.end(); ++p) {ar << p->first; ar << p->second;}
            ar << (int) m_sbmR.size();
            for (p = m_sbmR.begin(); p!=m_sbmR.end(); ++p) {ar << p->first; ar << p->second;}
            ar << (int) m_sbmP.size();
            for (p = m_sbmP.begin(); p!=m_sbmP.end(); ++p) {ar << p->first; ar << p->second;}
            ar << (int) m_solRi.size();
            for (p = m_solRi.begin(); p!=m_solRi.end(); ++p) {ar << p->first; ar << p->second;}
            ar << (int) m_solPi.size();
            for (p = m_solPi.begin(); p!=m_solPi.end(); ++p) {ar << p->first; ar << p->second;}
            ar << (int) m_solRe.size();
            for (p = m_solRe.begin(); p!=m_solRe.end(); ++p) {ar << p->first; ar << p->second;}
            ar << (int) m_solPe.size();
            for (p = m_solPe.begin(); p!=m_solPe.end(); ++p) {ar << p->first; ar << p->second;}
        }
        else
        {
            // restore pointers
            if (m_pFwd) m_pFwd->m_pReact = this;
            if (m_pRev) m_pRev->m_pReact = this;
            
            ar >> m_nsol >> m_vR >> m_vP >> m_v >> m_Vovr >> m_NSOL;
            int size, id, vR;
            ar >> size;
            for (int i=0; i<size; ++i)
            {
                ar >> id; ar >> vR;
                SetStoichiometricCoefficient(m_solR, id, vR);
            }
            ar >> size;
            for (int i=0; i<size; ++i)
            {
                ar >> id; ar >> vR;
                SetStoichiometricCoefficient(m_solP, id, vR);
            }
            ar >> size;
            for (int i=0; i<size; ++i)
            {
                ar >> id; ar >> vR;
                SetStoichiometricCoefficient(m_sbmR, id, vR);
            }
            ar >> size;
            for (int i=0; i<size; ++i)
            {
                ar >> id; ar >> vR;
                SetStoichiometricCoefficient(m_sbmP, id, vR);
            }
            ar >> size;
            for (int i=0; i<size; ++i)
            {
                ar >> id; ar >> vR;
                SetStoichiometricCoefficient(m_solRi, id, vR);
            }
            ar >> size;
            for (int i=0; i<size; ++i)
            {
                ar >> id; ar >> vR;
                SetStoichiometricCoefficient(m_solPi, id, vR);
            }
            ar >> size;
            for (int i=0; i<size; ++i)
            {
                ar >> id; ar >> vR;
                SetStoichiometricCoefficient(m_solRe, id, vR);
            }
            ar >> size;
            for (int i=0; i<size; ++i)
            {
                ar >> id; ar >> vR;
                SetStoichiometricCoefficient(m_solPe, id, vR);
            }
        }
    }
}
