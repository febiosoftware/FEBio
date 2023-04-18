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
#include "FESolutesMaterial.h"
#include <FECore/FECoreKernel.h>
#include <FECore/DumpStream.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FESolutesMaterial, FEMaterial)
	ADD_PROPERTY(m_pSolute, "solute", FEProperty::Optional);
    ADD_PROPERTY(m_pOsmC  , "osmotic_coefficient");
    ADD_PROPERTY(m_pReact , "reaction"           , FEProperty::Optional);
END_FECORE_CLASS();

//============================================================================
// FEFluidSolutesMaterialPoint
//============================================================================
FESolutesMaterial::Point::Point(FEMaterialPointData* pt) : FEMaterialPointData(pt)
{
	m_vft = vec3d(0.0, 0.0, 0.0);
	m_JfdotoJf = 0.0;
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FESolutesMaterial::Point::Copy()
{
	FESolutesMaterial::Point* pt = new FESolutesMaterial::Point(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FESolutesMaterial::Point::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);
	ar & m_vft & m_JfdotoJf;
    ar & m_nsol;
    ar & m_c & m_ca & m_gradc & m_j & m_cdot & m_k & m_dkdJ;
    ar & m_dkdc;
}

//-----------------------------------------------------------------------------
void FESolutesMaterial::Point::Init()
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

	FEMaterialPointData::Init();
}

//-----------------------------------------------------------------------------
double FESolutesMaterial::Point::Osmolarity() const
{
    double ew = 0.0;
    for (int isol = 0; isol < (int)m_ca.size(); ++isol)
    {
        ew += m_ca[isol];
    }
    return ew;
}


//============================================================================
// FESolutesMaterial
//============================================================================

//-----------------------------------------------------------------------------
//! FEFluidSolutes constructor

FESolutesMaterial::FESolutesMaterial(FEModel* pfem) : FEMaterial(pfem)
{
    m_Rgas = 0; m_Tabs = 0; m_Fc = 0;
    m_pOsmC = 0;
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPointData* FESolutesMaterial::CreateMaterialPointData()
{
    FEFluidMaterialPoint* fpt = new FEFluidMaterialPoint();
    return new FESolutesMaterial::Point(fpt);
}

//-----------------------------------------------------------------------------
// initialize
bool FESolutesMaterial::Init()
{
    // set the solute IDs first, since they are referenced in FESolute::Init()
    for (int i = 0; i<Solutes(); ++i) {
        m_pSolute[i]->SetSoluteLocalID(i);
    }
    
    // call the base class.
    // This also initializes all properties
    if (FEMaterial::Init() == false) return false;
    
    int zmin = 0, zmax = 0;
    
    m_Rgas = GetGlobalConstant("R");
    m_Tabs = GetGlobalConstant("T");
    m_Fc   = GetGlobalConstant("Fc");
    
    if (m_Rgas <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section"); return false; }
    if (m_Tabs <= 0) { feLogError("A positive absolute temperature T must be defined in Globals section");     return false; }
    if ((zmin || zmax) && (m_Fc <= 0)) {
        feLogError("A positive Faraday constant Fc must be defined in Globals section");
        return false;
    }
    
    return true;
}

//-----------------------------------------------------------------------------
void FESolutesMaterial::Serialize(DumpStream& ar)
{
    FEMaterial::Serialize(ar);
    
    if (ar.IsShallow()) return;
    
    ar & m_Rgas & m_Tabs & m_Fc;
}

//-----------------------------------------------------------------------------
//! partition coefficient
double FESolutesMaterial::PartitionCoefficient(FEMaterialPoint& pt, const int sol)
{
    
    // solubility
    double khat = m_pSolute[sol]->m_pSolub->Solubility(pt);
    double kappa = khat;
    
    return kappa;
}

//-----------------------------------------------------------------------------
//! partition coefficients and their derivatives
void FESolutesMaterial::PartitionCoefficientFunctions(FEMaterialPoint& mp, vector<double>& kappa,
                                                   vector<double>& dkdJ,
                                                   vector< vector<double> >& dkdc)
{
    int isol, jsol;
    
    FEFluidMaterialPoint& fpt = *(mp.ExtractData<FEFluidMaterialPoint>());
    FESolutesMaterial::Point& spt = *mp.ExtractData<FESolutesMaterial::Point>();
    
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
//! actual concentration
double FESolutesMaterial::Concentration(FEMaterialPoint& pt, const int sol)
{
	FESolutesMaterial::Point& spt = *pt.ExtractData<FESolutesMaterial::Point>();
    
    // effective concentration
    double c = spt.m_c[sol];
    
    return c;
}

//-----------------------------------------------------------------------------
//! actual concentration
double FESolutesMaterial::ConcentrationActual(FEMaterialPoint& pt, const int sol)
{
    FESolutesMaterial::Point& spt = *pt.ExtractData<FESolutesMaterial::Point>();
    
    // effective concentration
    double ca = spt.m_c[sol];
    
    // partition coefficient
    double kappa = PartitionCoefficient(pt, sol);
    
    ca = kappa*ca;
    
    return ca;
}

//-----------------------------------------------------------------------------
//! actual fluid pressure
double FESolutesMaterial::PressureActual(FEMaterialPoint& pt)
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

vec3d FESolutesMaterial::SoluteFlux(FEMaterialPoint& pt, const int sol)
{
	FESolutesMaterial::Point& spt = *pt.ExtractData<FESolutesMaterial::Point>();
    
    // concentration gradient
    vec3d gradc = spt.m_gradc[sol];
    
    // solute free diffusivity
    double D0 = m_pSolute[sol]->m_pDiff->Free_Diffusivity(pt);
    
    double c = spt.m_c[sol];
    vec3d v = spt.m_vft;
    
    // solute flux j
    vec3d j = -gradc*D0 + v*c;
    
    return j;
}
