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
void FEMultiphasicStandard::Init()
{
	FEMultiphasic::Init();
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEMultiphasicStandard::CreateMaterialPointData()
{
	return new FESolutesMaterialPoint(new FEBiphasicMaterialPoint(m_pSolid->CreateMaterialPointData()));
}

//-----------------------------------------------------------------------------
void FEMultiphasicStandard::UpdateSolidBoundMolecules(FEMaterialPoint& mp, const double dt)
{
    // check if this mixture includes chemical reactions
    int nreact = (int)Reactions();
    if (nreact) {
        // for chemical reactions involving solid-bound molecules,
        // update their concentration
		// multiphasic material point data
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
		FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
        
        double phi0 = ppt.m_phi0;
        int nsbm = SBMs();
        int nsol = Solutes();
        for (int isbm=0; isbm<nsbm; ++isbm) {
            spt.m_sbmrhat[isbm] = 0;
            // combine the molar supplies from all the reactions
            for (int k=0; k<nreact; ++k) {
                double zetahat = GetReaction(k)->ReactionSupply(mp);
                double v = GetReaction(k)->m_v[nsol+isbm];
                // remember to convert from molar supply to referential mass supply
                spt.m_sbmrhat[isbm] += (pt.m_J-phi0)*SBMMolarMass(isbm)*v*zetahat;
            }
            // perform the time integration (Euler's method)
            spt.m_sbmr[isbm] = spt.m_sbmrp[isbm] + dt*spt.m_sbmrhat[isbm];
            // check bounds
            if (spt.m_sbmr[isbm] < GetSBM(isbm)->m_rhomin)
                spt.m_sbmr[isbm] = GetSBM(isbm)->m_rhomin;
            if ((GetSBM(isbm)->m_rhomax > 0) && (spt.m_sbmr[isbm] > GetSBM(isbm)->m_rhomax))
                spt.m_sbmr[isbm] = GetSBM(isbm)->m_rhomax;
        }
    }
}
