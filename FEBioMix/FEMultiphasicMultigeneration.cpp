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
#include "FEMultiphasicMultigeneration.h"

// Material parameters for the FEMultiphasicMultigeneration material
BEGIN_FECORE_CLASS(FEMultiphasicMultigeneration, FEMultiphasic)
    ADD_PARAMETER(m_gtime   , FE_RANGE_GREATER(0.0), "gen_time");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! FEMultiphasicMultigeneration.h constructor
FEMultiphasicMultigeneration::FEMultiphasicMultigeneration(FEModel* pfem) : FEMultiphasic(pfem)
{
    m_gtime = 0;
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEMultiphasicMultigeneration::CreateMaterialPointData()
{
	return new FEMultigenSBMMaterialPoint
    (this, new FESolutesMaterialPoint
     (new FEBiphasicMaterialPoint
      (m_pSolid->CreateMaterialPointData())
      )
     );
}

//--------------------------------------------------------------------------------
// Check if time t constitutes a new generation and return that generation
int FEMultiphasicMultigeneration::CheckGeneration(const double t)
{
    // each generation has the same length of time
    // there is always at least one generation
	int ngen = (int)(t/m_gtime)+1;
	return ngen;
}

//--------------------------------------------------------------------------------
// Return the starting time of a generation
double FEMultiphasicMultigeneration::GetGenerationTime(const int igen)
{
	return (igen-1)*m_gtime;
}

//-----------------------------------------------------------------------------
void FEMultiphasicMultigeneration::UpdateSolidBoundMolecules(FEMaterialPoint& mp)
{
    double dt = CurrentTimeIncrement();
    
    // check if this mixture includes chemical reactions
    int nreact = (int)Reactions();
    int mreact = (int)MembraneReactions();
    if (nreact) {
        // for chemical reactions involving solid-bound molecules,
        // update their concentration
        // multiphasic material point data
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
        FEMultigenSBMMaterialPoint& mpt = *(mp.ExtractData<FEMultigenSBMMaterialPoint>());
        
        double phi0 = ppt.m_phi0t;
        int nsbm = SBMs();
        int nsol = Solutes();
        int ngen = mpt.m_ngen;
        
        for (int isbm=0; isbm<nsbm; ++isbm) {
            
            // initialize referential mass density supply
            spt.m_sbmrhat[isbm] = 0;
            // initialize referential mass density
            spt.m_sbmr[isbm] = spt.m_sbmrp[isbm];
            
            // evaluate mass fraction of each generation before SBM update
            vector<double> mf(ngen,0.0);
            for (int igen=0; igen<ngen; ++igen) {
                mf[igen] = (spt.m_sbmr[isbm] > 0) ? mpt.m_gsbmr[igen][isbm]/spt.m_sbmr[isbm] : 0;
                // initialize generational referential mass density
                mpt.m_gsbmr[igen][isbm] = mpt.m_gsbmrp[igen][isbm];
            }
            
            // for each reaction
            for (int k=0; k<nreact; ++k) {
                
                // evaluate the molar supply for this SBM
                double zetahat = GetReaction(k)->ReactionSupply(mp);
                double v = GetReaction(k)->m_v[nsol+isbm];
                
                // remember to convert from molar supply to referential mass supply
                double sbmrhat = (pt.m_J-phi0)*SBMMolarMass(isbm)*v*zetahat;
                
                // combine the molar supplies from all the reactions
                spt.m_sbmrhat[isbm] += sbmrhat;
                
                // check if mass is added or removed
                if (sbmrhat >= 0) {
                    
                    // mass is added only to the current generation
                    // perform the time integration (Euler's method)
                    double dsbmr = dt*sbmrhat;
                    
                    // add this mass increment to the current generation
                    mpt.m_gsbmr[ngen-1][isbm] += dsbmr;
                    spt.m_sbmr[isbm] += dsbmr;
                    
                    // check bounds
                    if ((GetSBM(isbm)->m_rhomax > 0) && (spt.m_sbmr[isbm] > GetSBM(isbm)->m_rhomax)) {
                        dsbmr = GetSBM(isbm)->m_rhomax - spt.m_sbmr[isbm];
                        mpt.m_gsbmr[ngen-1][isbm] += dsbmr;
                        spt.m_sbmr[isbm] += dsbmr;
                    }
                }
                else
                {
                    
                    // mass is removed from all the generations in proportion to mass fraction
                    // perform the time integration (Euler's method)
                    double dsbmr = dt*sbmrhat;
                    
                    // add this (negative) weighted mass increment to all generations
                    for (int igen=0; igen<ngen; ++igen)
                        mpt.m_gsbmr[igen][isbm] += mf[igen]*dsbmr;
                    spt.m_sbmr[isbm] += dsbmr;
                    
                    // check bounds
                    if (spt.m_sbmr[isbm] < GetSBM(isbm)->m_rhomin) {
                        dsbmr = GetSBM(isbm)->m_rhomin - spt.m_sbmr[isbm];
                        for (int igen=0; igen<ngen; ++igen)
                            mpt.m_gsbmr[igen][isbm] += mf[igen]*dsbmr;
                        spt.m_sbmr[isbm] += dsbmr;
                    }
                }
            }
            // for each membrane reaction
            for (int k=0; k<mreact; ++k) {
                
                // evaluate the molar supply for this SBM
                double zetahat = GetMembraneReaction(k)->ReactionSupply(mp);
                double v = GetMembraneReaction(k)->m_v[nsol+isbm];
                
                // remember to convert from molar supply to referential mass supply
                double sbmrhat = (pt.m_J-phi0)*SBMMolarMass(isbm)*v*zetahat;
                
                // combine the molar supplies from all the reactions
                spt.m_sbmrhat[isbm] += sbmrhat;
                
                // check if mass is added or removed
                if (sbmrhat >= 0) {
                    
                    // mass is added only to the current generation
                    // perform the time integration (Euler's method)
                    double dsbmr = dt*sbmrhat;
                    
                    // add this mass increment to the current generation
                    mpt.m_gsbmr[ngen-1][isbm] += dsbmr;
                    spt.m_sbmr[isbm] += dsbmr;
                    
                    // check bounds
                    if ((GetSBM(isbm)->m_rhomax > 0) && (spt.m_sbmr[isbm] > GetSBM(isbm)->m_rhomax)) {
                        dsbmr = GetSBM(isbm)->m_rhomax - spt.m_sbmr[isbm];
                        mpt.m_gsbmr[ngen-1][isbm] += dsbmr;
                        spt.m_sbmr[isbm] += dsbmr;
                    }
                }
                else
                {
                    
                    // mass is removed from all the generations in proportion to mass fraction
                    // perform the time integration (Euler's method)
                    double dsbmr = dt*sbmrhat;
                    
                    // add this (negative) weighted mass increment to all generations
                    for (int igen=0; igen<ngen; ++igen)
                        mpt.m_gsbmr[igen][isbm] += mf[igen]*dsbmr;
                    spt.m_sbmr[isbm] += dsbmr;
                    
                    // check bounds
                    if (spt.m_sbmr[isbm] < GetSBM(isbm)->m_rhomin) {
                        dsbmr = GetSBM(isbm)->m_rhomin - spt.m_sbmr[isbm];
                        for (int igen=0; igen<ngen; ++igen)
                            mpt.m_gsbmr[igen][isbm] += mf[igen]*dsbmr;
                        spt.m_sbmr[isbm] += dsbmr;
                    }
                }
            }
        }
    }
}

//=============================================================================
//   FEMultigenSBMMaterialPoint
//=============================================================================

//-----------------------------------------------------------------------------
FEMaterialPointData* FEMultigenSBMMaterialPoint::Copy()
{
	FEMultigenSBMMaterialPoint* pt = new FEMultigenSBMMaterialPoint(*this);
    pt->m_ngen = m_ngen;
    pt->m_nsbm = m_nsbm;
	pt->m_pmat = m_pmat;
	pt->m_Fi = m_Fi;
	pt->m_Ji = m_Ji;
    pt->m_gsbmr = m_gsbmr;
    pt->m_gsbmrp = m_gsbmrp;
    pt->m_lsbmr = m_lsbmr;
    pt->m_tgen = m_tgen;
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
void FEMultigenSBMMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);
    ar & m_ngen & m_nsbm & m_tgen;
	ar & m_Fi;
	ar & m_Ji;
	ar & m_gsbmr & m_gsbmrp;
	ar & m_lsbmr;
}

//-----------------------------------------------------------------------------
void FEMultigenSBMMaterialPoint::Init()
{
	FEMaterialPointData::Init();
    
	m_Fi.clear();
	m_Ji.clear();
	m_tgen = 0.0;
        
	FESolutesMaterialPoint& spt = *((*this).ExtractData<FESolutesMaterialPoint>());
	m_nsbm = spt.m_nsbm;
	m_ngen = 1; // first generation starts at t=0
	mat3d Fi(mat3dd(1.0));
	double Ji = 1;
	m_Fi.push_back(Fi);
	m_Ji.push_back(Ji);
	m_lsbmr.resize(m_nsbm,0.0);
	m_gsbmr.push_back(m_lsbmr);
	m_gsbmrp.push_back(m_lsbmr);
}

//-----------------------------------------------------------------------------
void FEMultigenSBMMaterialPoint::Update(const FETimeInfo& timeInfo)
{
	FEMaterialPointData::Update(timeInfo);

	// get the time
	double t = timeInfo.currentTime;
        
	// Check if this constitutes a new generation
	int ngen = m_pmat->CheckGeneration(t);
	t = m_pmat->GetGenerationTime(ngen);
        
	if (t>m_tgen) {
		// new generation
		FEElasticMaterialPoint& pt = *((*this).ExtractData<FEElasticMaterialPoint>());
		FESolutesMaterialPoint& spt = *((*this).ExtractData<FESolutesMaterialPoint>());
            
		// push back F and J to define relative deformation gradient of this generation
		mat3d F = pt.m_F;
		double J = pt.m_J;
		m_Fi.push_back(F.inverse());
		m_Ji.push_back(1.0/J);
		m_lsbmr = spt.m_sbmr;           // sbmr content at start of generation
		vector<double> gsbmr(m_nsbm,0); // incremental sbmr content for this generation
		m_gsbmr.push_back(gsbmr);
		m_gsbmrp.push_back(gsbmr);
		m_tgen = t;
		++m_ngen;
	}
}
