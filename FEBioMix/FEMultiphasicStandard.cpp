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
        
        // Allocation memory for coefficient matrix, solution vector and right-hand-side
        matrix A(nsbm,nsbm);
        vector<double> x(nsbm), rhs(nsbm);

        // Solve for incremental SBM apparent referential densities for each reaction
        for (int k=0; k<nreact; ++k) {
            double zetahat = GetReaction(k)->ReactionSupply(mp);
            for (int isbm=0; isbm<nsbm; ++isbm) {
                double vi = GetReaction(k)->m_v[nsol+isbm];
                rhs[isbm] = (pt.m_J-phi0)*SBMMolarMass(isbm)*vi*zetahat*dt;
                for (int jsbm=0; jsbm<nsbm; ++jsbm) {
                    A(isbm,jsbm) = SBMMolarMass(isbm)*vi/SBMDensity(jsbm)*zetahat*dt;
                }
                A(isbm,isbm) += 1;
            }
            A.solve(rhs,x);
            // x contains incremental densities for this reaction, now update m_sbmr
            for (int isbm=0; isbm<nsbm; ++isbm) {
                spt.m_sbmr[isbm] += x[isbm];
                // check bounds
                if (spt.m_sbmr[isbm] < GetSBM(isbm)->m_rhomin)
                    spt.m_sbmr[isbm] = GetSBM(isbm)->m_rhomin;
                if ((GetSBM(isbm)->m_rhomax > 0) && (spt.m_sbmr[isbm] > GetSBM(isbm)->m_rhomax))
                    spt.m_sbmr[isbm] = GetSBM(isbm)->m_rhomax;
            }
        }
    }
}
