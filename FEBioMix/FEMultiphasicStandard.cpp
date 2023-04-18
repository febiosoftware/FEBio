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
#include "FEMultiphasicStandard.h"

//-----------------------------------------------------------------------------
//! FEMultiphasicStandard constructor
FEMultiphasicStandard::FEMultiphasicStandard(FEModel* pfem) : FEMultiphasic(pfem)
{
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEMultiphasicStandard::CreateMaterialPointData()
{
	return new FESolutesMaterialPoint(new FEBiphasicMaterialPoint(m_pSolid->CreateMaterialPointData()));
}

//-----------------------------------------------------------------------------
// call this function from shell domains only
void FEMultiphasicStandard::UpdateSolidBoundMolecules(FEMaterialPoint& mp)
{
    double dt = CurrentTimeIncrement();
    
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
        FEShellElement* sel = dynamic_cast<FEShellElement*>(mp.m_elem);
        double h = (sel) ? sel->Evaluate(sel->m_ht, mp.m_index) : 0;   // shell thickness

        double phi0 = ppt.m_phi0t;
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
                spt.m_sbmrhat[isbm] += pt.m_J/h*SBMMolarMass(isbm)*v*zetahat;
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

