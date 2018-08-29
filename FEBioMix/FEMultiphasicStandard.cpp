//
//  FEMultiphasicStandard.cpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 2/23/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEMultiphasicStandard.h"

//-----------------------------------------------------------------------------
//! FEMultiphasicStandard constructor
FEMultiphasicStandard::FEMultiphasicStandard(FEModel* pfem) : FEMultiphasic(pfem)
{
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEMultiphasicStandard::CreateMaterialPointData()
{
	return new FESolutesMaterialPoint(new FEBiphasicMaterialPoint(m_pSolid->CreateMaterialPointData()));
}

//-----------------------------------------------------------------------------
// call this function from shell domains only
void FEMultiphasicStandard::UpdateSolidBoundMolecules(FEMaterialPoint& mp)
{
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    // check if this mixture includes chemical reactions
    int nreact = (int)Reactions();
    int mreact = (int)MembraneReactions();
    if (nreact || mreact) {
        // for chemical reactions involving solid-bound molecules,
        // update their concentration
		// multiphasic material point data
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
		FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
        
        double phi0 = ppt.m_phi0;
        int nsbm = SBMs();
        int nsol = Solutes();
        // create a temporary container for spt.m_sbmr so that this variable remains
        // unchanged as the reaction rates get calculated for each sbm.
        vector<double> sbmr = spt.m_sbmr;
        for (int isbm=0; isbm<nsbm; ++isbm) {
            spt.m_sbmrhat[isbm] = 0;
            // combine the molar supplies from all the reactions
            for (int k=0; k<nreact; ++k) {
                double zetahat = GetReaction(k)->ReactionSupply(mp);
                double v = GetReaction(k)->m_v[nsol+isbm];
                // remember to convert from molar supply to referential mass supply
                spt.m_sbmrhat[isbm] += (pt.m_J-phi0)*SBMMolarMass(isbm)*v*zetahat;
            }
            for (int k=0; k<mreact; ++k) {
                double zetahat = GetMembraneReaction(k)->ReactionSupply(mp);
                double v = GetMembraneReaction(k)->m_v[nsol+isbm];
                // remember to convert from molar supply to referential mass supply
                spt.m_sbmrhat[isbm] += (pt.m_J-phi0)*SBMMolarMass(isbm)*v*zetahat;
            }
            // perform the time integration (midpoint rule)
            sbmr[isbm] = spt.m_sbmrp[isbm] + dt*(spt.m_sbmrhat[isbm]+spt.m_sbmrhatp[isbm])/2;
            // check bounds
            if (sbmr[isbm] < GetSBM(isbm)->m_rhomin)
                sbmr[isbm] = GetSBM(isbm)->m_rhomin;
            if ((GetSBM(isbm)->m_rhomax > 0) && (sbmr[isbm] > GetSBM(isbm)->m_rhomax))
                sbmr[isbm] = GetSBM(isbm)->m_rhomax;
        }
        // now update spt.m_sbmr
        spt.m_sbmr = sbmr;
    }
}

