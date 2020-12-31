/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include "FEFluidSolutes.h"
#include <FECore/FECoreKernel.h>
#include <FECore/DumpStream.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEFluidSolutes, FEMaterial)
// material properties
ADD_PARAMETER(m_diffMtmSupp , "dms");
ADD_PROPERTY(m_pFluid, "fluid");
ADD_PROPERTY(m_pOsmC  , "osmotic_coefficient");
ADD_PROPERTY(m_pSolute, "solute"             , FEProperty::Optional);
ADD_PROPERTY(m_pReact , "reaction"           , FEProperty::Optional);
END_FECORE_CLASS();

//============================================================================
// FEFluidSolutesMaterialPoint
//============================================================================
FEFluidSolutesMaterialPoint::FEFluidSolutesMaterialPoint(FEMaterialPoint* pt) : FEMaterialPoint(pt) {}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEFluidSolutesMaterialPoint::Copy()
{
    FEFluidSolutesMaterialPoint* pt = new FEFluidSolutesMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FEFluidSolutesMaterialPoint::Serialize(DumpStream& ar)
{
    FEMaterialPoint::Serialize(ar);
    ar & m_nsol;
    ar & m_c & m_ca & m_gradc & m_j & m_cdot & m_k & m_dkdJ;
    ar & m_dkdc;
}

//-----------------------------------------------------------------------------
void FEFluidSolutesMaterialPoint::Init()
{
    m_nsol = 0;
    m_c.clear();
    m_ca.clear();
    m_gradc.clear();
    m_j.clear();
    m_cdot.clear();
    m_k.clear();
    m_dkdJ.clear();
    m_dkdc.clear();

    FEMaterialPoint::Init();
}

//============================================================================
// FEFluidSolutes
//============================================================================

//-----------------------------------------------------------------------------
//! FEFluidSolutes constructor

FEFluidSolutes::FEFluidSolutes(FEModel* pfem) : FEMaterial(pfem)
{
    m_pFluid = 0;
    m_Rgas = 0; m_Tabs = 0; m_Fc = 0;
    m_diffMtmSupp = 1.0;
    m_pOsmC = 0;
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPoint* FEFluidSolutes::CreateMaterialPointData()
{
    FEFluidMaterialPoint* fpt = new FEFluidMaterialPoint();
    return new FEFluidSolutesMaterialPoint(fpt);
}

//-----------------------------------------------------------------------------
void FEFluidSolutes::AddChemicalReaction(FEChemicalReaction* pcr)
{
    m_pReact.push_back(pcr);
}

//-----------------------------------------------------------------------------
// initialize
bool FEFluidSolutes::Init()
{
    // we first have to set the parent material
    // TODO: This seems redundant since each material already has a pointer to its parent
    for (int i=0; i<Reactions(); ++i)
    {
        m_pReact[i]->m_pFS = this;
    }
    
    // set the solute IDs first, since they are referenced in FESolute::Init()
    for (int i = 0; i<Solutes(); ++i) {
        m_pSolute[i]->SetSoluteLocalID(i);
    }
    
    // call the base class.
    // This also initializes all properties
    if (FEMaterial::Init() == false) return false;
    
    int zmin = 0, zmax = 0;
    
    m_Rgas = GetFEModel()->GetGlobalConstant("R");
    m_Tabs = GetFEModel()->GetGlobalConstant("T");
    m_Fc   = GetFEModel()->GetGlobalConstant("Fc");
    
    if (m_Rgas <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section"); return false; }
    if (m_Tabs <= 0) { feLogError("A positive absolute temperature T must be defined in Globals section");     return false; }
    if ((zmin || zmax) && (m_Fc <= 0)) {
        feLogError("A positive Faraday constant Fc must be defined in Globals section");
        return false;
    }
    
    return true;
}

//-----------------------------------------------------------------------------
void FEFluidSolutes::Serialize(DumpStream& ar)
{
    FEMaterial::Serialize(ar);
    if (ar.IsShallow()) return;
    
    ar & m_Rgas & m_Tabs & m_Fc;
    
    if (ar.IsLoading())
    {
        // restore the m_pMP pointers for reactions
        int NR = (int) m_pReact.size();
        for (int i=0; i<NR; ++i) m_pReact[i]->m_pFS = this;
    }
}

//-----------------------------------------------------------------------------
//! partition coefficient
double FEFluidSolutes::PartitionCoefficient(FEMaterialPoint& pt, const int sol)
{
    
    // solubility
    double khat = m_pSolute[sol]->m_pSolub->Solubility(pt);
    double kappa = khat;
    
    return kappa;
}

//-----------------------------------------------------------------------------
//! partition coefficients and their derivatives
void FEFluidSolutes::PartitionCoefficientFunctions(FEMaterialPoint& mp, vector<double>& kappa,
                                                  vector<double>& dkdJ,
                                                  vector< vector<double> >& dkdc)
{
    int isol, jsol;
    
    FEFluidMaterialPoint& fpt = *(mp.ExtractData<FEFluidMaterialPoint>());
    FEFluidSolutesMaterialPoint& spt = *(mp.ExtractData<FEFluidSolutesMaterialPoint>());
    
    const int nsol = (int)m_pSolute.size();
    
    vector<double> c(nsol);
    vector<int> z(nsol);
    vector<double> khat(nsol);
    vector<double> dkhdJ(nsol);
    vector< vector<double> > dkhdc(nsol, vector<double>(nsol));
    kappa.resize(nsol);
    
    for (isol=0; isol<nsol; ++isol) {
        // get the effective concentration, its gradient and its time derivative
        c[isol] = spt.m_c[isol];
        // evaluate the solubility and its derivatives w.r.t. J and c
        khat[isol] = m_pSolute[isol]->m_pSolub->Solubility(mp);
        dkhdJ[isol] = m_pSolute[isol]->m_pSolub->Tangent_Solubility_Strain(mp);
        for (jsol=0; jsol<nsol; ++jsol) {
            dkhdc[isol][jsol] = m_pSolute[isol]->m_pSolub->Tangent_Solubility_Concentration(mp,jsol);
            dkdc[isol][jsol]=dkhdc[isol][jsol];
        }
        kappa[isol] = khat[isol];
        dkdJ[isol] = dkhdJ[isol];
    }

}

//-----------------------------------------------------------------------------
//! effective concentration
double FEFluidSolutes::Concentration(FEMaterialPoint& pt, const int sol)
{
    FEFluidSolutesMaterialPoint& spt = *pt.ExtractData<FEFluidSolutesMaterialPoint>();
    
    // effective concentration
    double c = spt.m_c[sol];
    
    return c;
}

//-----------------------------------------------------------------------------
//! actual concentration
double FEFluidSolutes::ConcentrationActual(FEMaterialPoint& pt, const int sol)
{
    FEFluidSolutesMaterialPoint& spt = *pt.ExtractData<FEFluidSolutesMaterialPoint>();
    
    // effective concentration
    double ca = spt.m_c[sol];
    
    // partition coefficient
    double kappa = PartitionCoefficient(pt, sol);
    
    ca = kappa*ca;
    
    return ca;
}

//-----------------------------------------------------------------------------
//! actual fluid pressure
double FEFluidSolutes::PressureActual(FEMaterialPoint& pt)
{
    int i;
    
    FEFluidMaterialPoint& fpt = *pt.ExtractData<FEFluidMaterialPoint>();
    const int nsol = (int)m_pSolute.size();
    
    // effective pressure
    double p = fpt.m_pf;
    
    // effective concentration
    vector<double> c(nsol);
    for (i=0; i<nsol; ++i)
        c[i] = Concentration(pt, i);
    
    // osmotic coefficient
    double osmc = m_pOsmC->OsmoticCoefficient(pt);
    
    // actual pressure
    double pa = 0;
    for (i=0; i<nsol; ++i) pa += c[i];
    pa = p + m_Rgas*m_Tabs*osmc*pa;
    
    return pa;
}

//-----------------------------------------------------------------------------
//! Calculate solute molar flux

vec3d FEFluidSolutes::SoluteFlux(FEMaterialPoint& pt, const int sol)
{
    FEFluidSolutesMaterialPoint& spt = *pt.ExtractData<FEFluidSolutesMaterialPoint>();
    FEFluidMaterialPoint& fpt = *pt.ExtractData<FEFluidMaterialPoint>();
    
    // concentration gradient
    vec3d gradc = spt.m_gradc[sol];
    
    // solute free diffusivity
    double D0 = m_pSolute[sol]->m_pDiff->Free_Diffusivity(pt);
    
    double c = spt.m_c[sol];
    vec3d v = fpt.m_vft;
    
    // solute flux j
    vec3d j = -gradc*D0 + v*c;
    
    return j;
}
