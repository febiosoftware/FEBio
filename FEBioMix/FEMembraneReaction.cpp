//
//  FEMembraneReaction.cpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 3/4/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#include "stdafx.h"
#include "FEMembraneReaction.h"
#include "FECore/FEElementTraits.h"
#include "FECore/DOFS.h"
#include "FECore/FEModel.h"
#include "FECore/FECoreKernel.h"
#include "FEMultiphasic.h"
#include <stdlib.h>


//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEMembraneReaction, FEReaction)
	ADD_PARAMETER(m_Vbar  , "Vbar");
	ADD_PARAMETER(m_vRtmp , "vR"  );
	ADD_PARAMETER(m_vPtmp , "vP"  );
	ADD_PARAMETER(m_vRitmp, "vRi");
	ADD_PARAMETER(m_vPitmp, "vPi");
	ADD_PARAMETER(m_vRetmp, "vRe");
	ADD_PARAMETER(m_vPetmp, "vPe");

	// set material properties
	ADD_PROPERTY(m_pFwd, "forward_rate", FEProperty::Optional);
	ADD_PROPERTY(m_pRev, "reverse_rate", FEProperty::Optional);

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMembraneReaction::FEMembraneReaction(FEModel* pfem) : FEReaction(pfem)
{
    // additional initializations
    m_Vovr = false;

	m_pFwd = m_pRev = 0;
}

//-----------------------------------------------------------------------------
FESoluteData* FEMembraneReaction::FindSoluteData(int nid)
{
    FEModel& fem = *GetFEModel();
    int N = GetFEModel()->GlobalDataItems();
    for (int i=0; i<N; ++i)
    {
        FESoluteData* psd = dynamic_cast<FESoluteData*>(fem.GetGlobalData(i));
        if (psd && (psd->GetID()-1 == nid)) return psd;
    }
    return 0;
}

//-----------------------------------------------------------------------------
bool FEMembraneReaction::Init()
{
    // initialize base class
    FEReaction::Init();
    
    // set the parents for the reaction rates
    if (m_pFwd) m_pFwd->m_pReact = this;
    if (m_pRev) m_pRev->m_pReact = this;
    
    //************* reactants and products in multiphasic domain **************
    
    // initialize the reaction coefficients
    const int nsol = m_pMP->Solutes();
    const int nsbm = m_pMP->SBMs();
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
        int sid = m_pMP->GetSolute(isol)->GetSoluteID();
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
        int sid = m_pMP->GetSBM(isbm)->GetSBMID();
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
    int MAX_CDOFS = fedofs.GetVariableSize("concentration");
    m_NSOL = MAX_CDOFS;

    // initialize the stoichiometric coefficients to zero
    m_vRi.assign(MAX_CDOFS, 0); m_vRe.assign(MAX_CDOFS, 0);
    m_vPi.assign(MAX_CDOFS, 0); m_vPe.assign(MAX_CDOFS, 0);
    m_vi.assign(MAX_CDOFS, 0); m_ve.assign(MAX_CDOFS, 0);
    
    // cycle through all the solutes in the mixture and determine
    // if they participate in this reaction
    for (int ISOL = 0; ISOL<MAX_CDOFS; ++ISOL) {
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
    for (int ISOL = 0; ISOL<MAX_CDOFS; ++ISOL) {
        m_vi[ISOL] = m_vPi[ISOL] - m_vRi[ISOL];
        m_ve[ISOL] = m_vPe[ISOL] - m_vRe[ISOL];
    }

    //************** continue with all reactants and products ***************
    
    // evaluate the weighted molar volume of reactants and products
    if (!m_Vovr) {
        m_Vbar = 0;
        for (int isol = 0; isol<nsol; ++isol)
            m_Vbar += m_v[isol] * m_pMP->GetSolute(isol)->MolarMass() / m_pMP->GetSolute(isol)->Density();
        for (int isbm = 0; isbm<nsbm; ++isbm)
            m_Vbar += m_v[nsol + isbm] * m_pMP->GetSBM(isbm)->MolarMass() / m_pMP->GetSBM(isbm)->Density();
        for (int ISOL = 0; ISOL<MAX_CDOFS; ++ISOL) {
            FESoluteData* sd = FindSoluteData(ISOL);
            m_Vbar += m_vi[ISOL] * sd->m_M / sd->m_rhoT;
            m_Vbar += m_ve[ISOL] * sd->m_M / sd->m_rhoT;
        }
    }
    
    // check that the reaction satisfies electroneutrality
    int znet = 0;
    for (int isol = 0; isol<nsol; ++isol)
        znet += m_v[isol] * m_pMP->GetSolute(isol)->ChargeNumber();
    for (int isbm = 0; isbm<nsbm; ++isbm)
        znet += m_v[nsol + isbm] * m_pMP->GetSBM(isbm)->ChargeNumber();
    for (int ISOL = 0; ISOL<MAX_CDOFS; ++ISOL) {
        FESoluteData* sd = FindSoluteData(ISOL);
        znet += m_vi[ISOL] * sd->m_z;
        znet += m_ve[ISOL] * sd->m_z;
    }
    if (znet != 0) return fecore_error("membrane reaction must satisfy electroneutrality");
    
    return true;
}

//-----------------------------------------------------------------------------
void FEMembraneReaction::SetParameter(FEParam& p)
{
    if (strcmp(p.name(), "Vbar") == 0)
    {
        m_Vovr = true;
    }
}

//-----------------------------------------------------------------------------
bool FEMembraneReaction::SetParameterAttribute(FEParam& p, const char* szatt, const char* szval)
{
    // get number of DOFS
    DOFS& fedofs = GetFEModel()->GetDOFS();
    int MAX_CDOFS = fedofs.GetVariableSize("concentration");
    
    if (strcmp(p.name(), "vR") == 0)
    {
        if (strcmp(szatt, "sbm") == 0)
        {
            int id = atoi(szval) - 1;
            if (id < 0) return false;
            SetStoichiometricCoefficient(m_sbmR, id, m_vRtmp);
            return true;
        }
        if (strcmp(szatt, "sol") == 0)
        {
            int id = atoi(szval) - 1;
            if ((id < 0) || (id >= MAX_CDOFS)) return false;
            SetStoichiometricCoefficient(m_solR, id, m_vRtmp);
            return true;
        }
    }
    else if (strcmp(p.name(), "vP") == 0)
    {
        if (strcmp(szatt, "sbm") == 0)
        {
            int id = atoi(szval) - 1;
            if (id < 0) return false;
            SetStoichiometricCoefficient(m_sbmP, id, m_vPtmp);
            return true;
        }
        if (strcmp(szatt, "sol") == 0)
        {
            int id = atoi(szval) - 1;
            if ((id < 0) || (id >= MAX_CDOFS)) return false;
            SetStoichiometricCoefficient(m_solP, id, m_vPtmp);
            return true;
        }
    }
    else if (strcmp(p.name(), "vRi") == 0)
    {
        if (strcmp(szatt, "sol") == 0)
        {
            int id = atoi(szval) - 1;
            if ((id < 0) || (id >= MAX_CDOFS)) return false;
            SetStoichiometricCoefficient(m_solRi, id, m_vRitmp);
            return true;
        }
    }
    else if (strcmp(p.name(), "vPi") == 0)
    {
        if (strcmp(szatt, "sol") == 0)
        {
            int id = atoi(szval) - 1;
            if ((id < 0) || (id >= MAX_CDOFS)) return false;
            SetStoichiometricCoefficient(m_solPi, id, m_vPitmp);
            return true;
        }
    }
    else if (strcmp(p.name(), "vRe") == 0)
    {
        if (strcmp(szatt, "sol") == 0)
        {
            int id = atoi(szval) - 1;
            if ((id < 0) || (id >= MAX_CDOFS)) return false;
            SetStoichiometricCoefficient(m_solRe, id, m_vRetmp);
            return true;
        }
    }
    else if (strcmp(p.name(), "vPe") == 0)
    {
        if (strcmp(szatt, "sol") == 0)
        {
            int id = atoi(szval) - 1;
            if ((id < 0) || (id >= MAX_CDOFS)) return false;
            SetStoichiometricCoefficient(m_solPe, id, m_vPetmp);
            return true;
        }
    }
    return false;
}

//-----------------------------------------------------------------------------
//! Data serialization
void FEMembraneReaction::Serialize(DumpStream& ar)
{
    FEMaterial::Serialize(ar);
    
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
