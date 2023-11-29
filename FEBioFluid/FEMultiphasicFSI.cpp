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
#include "FEBiphasicFSI.h"
#include "FEMultiphasicFSI.h"
#include "FEFluidFSI.h"
#include <FECore/FECoreKernel.h>
#include <FECore/DumpStream.h>
#include "FECore/FEModel.h"
#include <FECore/log.h>
#include <FECore/tools.h>
#include <complex>
using namespace std;

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEMultiphasicFSI, FEBiphasicFSI)
ADD_PARAMETER(m_penalty, FE_RANGE_GREATER_OR_EQUAL(0.0), "penalty"      );
ADD_PARAMETER(m_diffMtmSupp , "dms");
ADD_PARAMETER(m_cFr    , "fixed_charge_density");
// material properties
ADD_PROPERTY(m_pOsmC  , "osmotic_coefficient");
ADD_PROPERTY(m_pSolute, "solute"             , FEProperty::Optional);
ADD_PROPERTY(m_pReact , "reaction"           , FEProperty::Optional);
END_FECORE_CLASS();

//============================================================================
// FEFSIMaterialPoint
//============================================================================
FEMultiphasicFSIMaterialPoint::FEMultiphasicFSIMaterialPoint(FEMaterialPointData* pt) : FEMaterialPointData(pt) {}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEMultiphasicFSIMaterialPoint::Copy()
{
    FEMultiphasicFSIMaterialPoint* pt = new FEMultiphasicFSIMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FEMultiphasicFSIMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);
    ar & m_nsol & m_psi & m_Ie & m_cF & m_pe;
    ar & m_c & m_ca & m_gradc & m_j & m_cdot & m_k & m_dkdJ;
    ar & m_dkdc;
}

//-----------------------------------------------------------------------------
void FEMultiphasicFSIMaterialPoint::Init()
{
    m_nsol = 0;
    m_psi = m_cF = 0;
    m_pe = 0;
    m_Ie = vec3d(0,0,0);
    m_c.clear();
    m_ca.clear();
    m_gradc.clear();
    m_j.clear();
    m_cdot.clear();
    m_k.clear();
    m_dkdJ.clear();
    m_dkdJJ.clear();
    m_dkdc.clear();
    m_dkdJc.clear();
    m_dkdcc.clear();
    
	FEMaterialPointData::Init();
}

//-----------------------------------------------------------------------------
double FEMultiphasicFSIMaterialPoint::Osmolarity() const
{
    double ew = 0.0;
    for (int isol = 0; isol < (int)m_ca.size(); ++isol)
    {
        ew += m_ca[isol];
    }
    return ew;
}

//============================================================================
// FEFluidFSI
//============================================================================

//-----------------------------------------------------------------------------
//! FEFluidFSI constructor

FEMultiphasicFSI::FEMultiphasicFSI(FEModel* pfem) : FEBiphasicFSI(pfem)
{
    m_Rgas = 0; m_Tabs = 0; m_Fc = 0;
    m_cFr = 0.0;
    m_diffMtmSupp = 1.0;
    m_penalty = 1;
    m_pOsmC = 0;
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPointData* FEMultiphasicFSI::CreateMaterialPointData()
{
    FEFluidMaterialPoint* fpt = new FEFluidMaterialPoint(m_pSolid->CreateMaterialPointData());
    FEFSIMaterialPoint* fst = new FEFSIMaterialPoint(fpt);
    FEBiphasicFSIMaterialPoint* bfpt = new FEBiphasicFSIMaterialPoint(fst);
    FEMultiphasicFSIMaterialPoint* mfpt = new FEMultiphasicFSIMaterialPoint(bfpt);
    
    return mfpt;
}

//-----------------------------------------------------------------------------
// initialize
bool FEMultiphasicFSI::Init()
{
    // set the solute IDs first, since they are referenced in FESolute::Init()
    for (int i = 0; i<Solutes(); ++i) {
        m_pSolute[i]->SetSoluteLocalID(i);
    }
    
    // call the base class.
    // This also initializes all properties
    if (FEBiphasicFSI::Init() == false) return false;
    
    // Determine how to solve for the electric potential psi
    int isol;
    int zmin = 0, zmax = 0, z;
    for (isol=0; isol<(int)m_pSolute.size(); ++isol) {
        z = m_pSolute[isol]->ChargeNumber();
        if (z < zmin) zmin = z;
        if (z > zmax) zmax = z;
    }
    m_zmin = zmin;
    m_ndeg = zmax - zmin;    // polynomial degree
    
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
void FEMultiphasicFSI::Serialize(DumpStream& ar)
{
    FEMaterial::Serialize(ar);
    if (ar.IsShallow()) return;
    
    ar & m_Rgas & m_Tabs & m_Fc;
    ar & m_zmin & m_ndeg;
}

//-----------------------------------------------------------------------------
//! The stress of a poro-elastic material is the sum of the fluid stress
//! and the elastic stress. Note that this function is declared in the base class
//! so you do not have to reimplement it in a derived class, unless additional
//! pressure terms are required.

mat3ds FEMultiphasicFSI::Stress(FEMaterialPoint& mp)
{
    // calculate solid material stress
    mat3ds s = m_pSolid->Stress(mp);
    
    // add fluid stress
    s += m_pFluid->GetViscous()->Stress(mp);
    
    //add actual pressure
    s += -mat3dd(1.0)*PressureActual(mp);
    
    return s;
}

//-----------------------------------------------------------------------------
//! Return the permeability tensor as a double array

void FEMultiphasicFSI::Diffusivity(double d[3][3], FEMaterialPoint& pt, int sol)

{
    mat3ds dt = m_pSolute[sol]->m_pDiff->Diffusivity(pt);
    
    d[0][0] = dt.xx();
    d[1][1] = dt.yy();
    d[2][2] = dt.zz();
    d[0][1] = d[1][0] = dt.xy();
    d[1][2] = d[2][1] = dt.yz();
    d[2][0] = d[0][2] = dt.xz();
    
}

//-----------------------------------------------------------------------------
mat3ds FEMultiphasicFSI::Diffusivity(FEMaterialPoint& mp, int sol)
{
    return m_pSolute[sol]->m_pDiff->Diffusivity(mp);
}

//-----------------------------------------------------------------------------
tens4dmm FEMultiphasicFSI::Diffusivity_Tangent_Strain(FEMaterialPoint& mp, int isol)
{
    //Return 0 if diffusivity is 0 to avoid NAN
    if (m_pSolute[isol]->m_pDiff->Diffusivity(mp).xx() == 0.0 && m_pSolute[isol]->m_pDiff->Diffusivity(mp).xy() == 0.0 && m_pSolute[isol]->m_pDiff->Diffusivity(mp).xz() == 0.0 && m_pSolute[isol]->m_pDiff->Diffusivity(mp).yy() == 0.0 && m_pSolute[isol]->m_pDiff->Diffusivity(mp).yz() == 0.0 && m_pSolute[isol]->m_pDiff->Diffusivity(mp).zz() == 0.0)
        return tens4dmm(0.0);
    else
        return m_pSolute[isol]->m_pDiff->Tangent_Diffusivity_Strain(mp);
}

//-----------------------------------------------------------------------------
mat3ds FEMultiphasicFSI::Diffusivity_Tangent_Concentration(FEMaterialPoint& mp, int isol, int jsol)
{
    //Return 0 if diffusivity is 0 to avoid NAN
    if (m_pSolute[isol]->m_pDiff->Diffusivity(mp).xx() == 0.0 && m_pSolute[isol]->m_pDiff->Diffusivity(mp).xy() == 0.0 && m_pSolute[isol]->m_pDiff->Diffusivity(mp).xz() == 0.0 && m_pSolute[isol]->m_pDiff->Diffusivity(mp).yy() == 0.0 && m_pSolute[isol]->m_pDiff->Diffusivity(mp).yz() == 0.0 && m_pSolute[isol]->m_pDiff->Diffusivity(mp).zz() == 0.0)
        return mat3ds(0.0);
    else
        return m_pSolute[isol]->m_pDiff->Tangent_Diffusivity_Concentration(mp, jsol);
}

//-----------------------------------------------------------------------------
mat3ds FEMultiphasicFSI::InvDiffusivity(FEMaterialPoint& mp, int sol)
{
    //Return 0 for inverse diffusivity when diffusivity is set to 0.
    //Acts as if diffusivity if infinite.
    if (m_pSolute[sol]->m_pDiff->Diffusivity(mp).xx() == 0.0 && m_pSolute[sol]->m_pDiff->Diffusivity(mp).xy() == 0.0 && m_pSolute[sol]->m_pDiff->Diffusivity(mp).xz() == 0.0 && m_pSolute[sol]->m_pDiff->Diffusivity(mp).yy() == 0.0 && m_pSolute[sol]->m_pDiff->Diffusivity(mp).yz() == 0.0 && m_pSolute[sol]->m_pDiff->Diffusivity(mp).zz() == 0.0)
        return mat3ds(0.0);
    else
        return m_pSolute[sol]->m_pDiff->Diffusivity(mp).inverse();
}

//-----------------------------------------------------------------------------
//! Electric potential
double FEMultiphasicFSI::ElectricPotential(FEMaterialPoint& pt, const bool eform)
{
    // check if solution is neutral
    if (m_ndeg == 0) {
        if (eform) return 1.0;
        else return 0.0;
    }
    
    int i, j;
    
    // if not neutral, solve electroneutrality polynomial for zeta
    FEMultiphasicFSIMaterialPoint& set = *pt.ExtractData<FEMultiphasicFSIMaterialPoint>();
    const int nsol = (int)m_pSolute.size();
    double cF = FixedChargeDensity(pt);
    
    vector<double> c(nsol);        // effective concentration
    vector<double> khat(nsol);    // solubility
    vector<int> z(nsol);        // charge number
    for (i=0; i<nsol; ++i) {
        c[i] = set.m_c[i];
        khat[i] = m_pSolute[i]->m_pSolub->Solubility(pt);
        z[i] = m_pSolute[i]->ChargeNumber();
    }
    
    // evaluate polynomial coefficients
    const int n = m_ndeg;
    vector<double> a(n+1,0);
    if (m_zmin < 0) {
        for (i=0; i<nsol; ++i) {
            j = z[i] - m_zmin;
            a[j] += z[i]*khat[i]*c[i];
        }
        a[-m_zmin] = cF;
    } else {
        for (i=0; i<nsol; ++i) {
            j = z[i];
            a[j] += z[i]*khat[i]*c[i];
        }
        a[0] = cF;
    }
    
    // solve polynomial
    double psi = set.m_psi;        // use previous solution as initial guess
    double zeta = exp(-m_Fc*psi/m_Rgas/m_Tabs);
    if (!solvepoly(n, a, zeta)) {
        zeta = 1.0;
    }
    
    // Return exponential (non-dimensional) form if desired
    if (eform) return zeta;
    
    // Otherwise return dimensional value of electric potential
    psi = -m_Rgas*m_Tabs/m_Fc*log(zeta);
    
    return psi;
}

//-----------------------------------------------------------------------------
//! partition coefficient
double FEMultiphasicFSI::PartitionCoefficient(FEMaterialPoint& pt, const int sol)
{
    
    // solubility
    double khat = m_pSolute[sol]->m_pSolub->Solubility(pt);
    // charge number
    int z = m_pSolute[sol]->ChargeNumber();
    // electric potential
    double zeta = ElectricPotential(pt, true);
    double zz = pow(zeta, z);
    // partition coefficient
    double kappa = zz*khat;
    
    return kappa;
}

//-----------------------------------------------------------------------------
//! partition coefficients and their derivatives
void FEMultiphasicFSI::PartitionCoefficientFunctions(FEMaterialPoint& mp, vector<double>& kappa,
                                                   vector<double>& dkdJ,
                                                   vector< vector<double> >& dkdc)
{
    //TODO: Include dkdcc and dkdJc
    
    int isol, jsol, ksol;
    
    FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
    FEFluidMaterialPoint& fpt = *(mp.ExtractData<FEFluidMaterialPoint>());
    FEBiphasicFSIMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicFSIMaterialPoint>());
    FEMultiphasicFSIMaterialPoint& spt = *(mp.ExtractData<FEMultiphasicFSIMaterialPoint>());
    
    const int nsol = (int)m_pSolute.size();
    
    vector<double> c(nsol);
    vector<int> z(nsol);
    vector<double> khat(nsol);
    vector<double> dkhdJ(nsol);
    vector<double> dkhdJJ(nsol);
    vector< vector<double> > dkhdc(nsol, vector<double>(nsol));
    vector< vector<double> > dkhdJc(nsol, vector<double>(nsol));
    vector< vector< vector<double> > > dkhdcc(nsol, dkhdc);    // use dkhdc to initialize only
    vector<double> zz(nsol);
    kappa.resize(nsol);
    
    double den = 0;
    double num = 0;
    double zeta = ElectricPotential(mp, true);
    
    for (isol=0; isol<nsol; ++isol) {
        // get the effective concentration, its gradient and its time derivative
        c[isol] = spt.m_c[isol];
        // get the charge number
        z[isol] = m_pSolute[isol]->ChargeNumber();
        // evaluate the solubility and its derivatives w.r.t. J and c
        khat[isol] = m_pSolute[isol]->m_pSolub->Solubility(mp);
        dkhdJ[isol] = m_pSolute[isol]->m_pSolub->Tangent_Solubility_Strain(mp);
        dkhdJJ[isol] = m_pSolute[isol]->m_pSolub->Tangent_Solubility_Strain_Strain(mp);
        for (jsol=0; jsol<nsol; ++jsol) {
            dkhdc[isol][jsol] = m_pSolute[isol]->m_pSolub->Tangent_Solubility_Concentration(mp,jsol);
            dkhdJc[isol][jsol] = m_pSolute[isol]->m_pSolub->Tangent_Solubility_Strain_Concentration(mp,jsol);
            for (ksol=0; ksol<nsol; ++ksol) {
                dkhdcc[isol][jsol][ksol] =
                m_pSolute[isol]->m_pSolub->Tangent_Solubility_Concentration_Concentration(mp,jsol,ksol);
            }
        }
        zz[isol] = pow(zeta, z[isol]);
        kappa[isol] = zz[isol]*khat[isol];
        den += SQR(z[isol])*kappa[isol]*c[isol];
        num += pow((double)z[isol],3)*kappa[isol]*c[isol];
    }
    
    // get the charge density and its derivatives
    double J = ept.m_J;
    double phi0 = ppt.m_phi0;
    double cF = FixedChargeDensity(mp);
    double dcFdJ = -cF/(J - phi0);
    double dcFdJJ = 2*cF/SQR(J-phi0);
    
    // evaluate electric potential (nondimensional exponential form) and its derivatives
    // also evaluate partition coefficients and their derivatives
    double zidzdJ = 0;
    double zidzdJJ = 0, zidzdJJ1 = 0, zidzdJJ2 = 0;
    vector<double> zidzdc(nsol,0);
    vector<double> zidzdJc(nsol,0), zidzdJc1(nsol,0), zidzdJc2(nsol,0);
    vector< vector<double> > zidzdcc(nsol, vector<double>(nsol,0));
    vector< vector<double> > zidzdcc1(nsol, vector<double>(nsol,0));
    vector<double> zidzdcc2(nsol,0);
    double zidzdcc3 = 0;
    
    if (den > 0) {
        
        for (isol=0; isol<nsol; ++isol)
            zidzdJ += z[isol]*zz[isol]*dkhdJ[isol]*c[isol];
        zidzdJ = -(dcFdJ+zidzdJ)/den;
        
        for (isol=0; isol<nsol; ++isol) {
            for (jsol=0; jsol<nsol; ++jsol) {
                zidzdJJ1 += SQR(z[jsol])*c[jsol]*(z[jsol]*zidzdJ*kappa[jsol]+zz[jsol]*dkhdJ[jsol]);
                zidzdJJ2 += z[jsol]*zz[jsol]*c[jsol]*(zidzdJ*z[jsol]*dkhdJ[jsol]+dkhdJJ[jsol]);
                zidzdc[isol] += z[jsol]*zz[jsol]*dkhdc[jsol][isol]*c[jsol];
            }
            zidzdc[isol] = -(z[isol]*kappa[isol]+zidzdc[isol])/den;
            zidzdcc3 += pow(double(z[isol]),3)*kappa[isol]*c[isol];
        }
        zidzdJJ = zidzdJ*(zidzdJ-zidzdJJ1/den)-(dcFdJJ+zidzdJJ2)/den;
        
        for (isol=0; isol<nsol; ++isol) {
            for (jsol=0; jsol<nsol; ++jsol) {
                zidzdJc1[isol] += SQR(z[jsol])*c[jsol]*(zidzdc[isol]*z[jsol]*kappa[jsol]+zz[jsol]*dkhdc[jsol][isol]);
                zidzdJc2[isol] += z[jsol]*zz[jsol]*c[jsol]*(zidzdc[isol]*z[jsol]*dkhdJ[jsol]+dkhdJc[jsol][isol]);
                zidzdcc2[isol] += SQR(z[jsol])*zz[jsol]*c[jsol]*dkhdc[jsol][isol];
                for (ksol=0; ksol<nsol; ++ksol)
                    zidzdcc1[isol][jsol] += z[ksol]*zz[ksol]*c[ksol]*dkhdcc[ksol][isol][jsol];
            }
            zidzdJc[isol] = zidzdJ*(zidzdc[isol]-(SQR(z[isol])*kappa[isol] + zidzdJc1[isol])/den)
            -(z[isol]*zz[isol]*dkhdJ[isol] + zidzdJc2[isol])/den;
        }
        
        for (isol=0; isol<nsol; ++isol) {
            for (jsol=0; jsol<nsol; ++jsol) {
                zidzdcc[isol][jsol] = zidzdc[isol]*zidzdc[jsol]*(1 - zidzdcc3/den)
                - zidzdcc1[isol][jsol]/den
                - z[isol]*(z[isol]*kappa[isol]*zidzdc[jsol]+zz[isol]*dkhdc[isol][jsol])/den
                - z[jsol]*(z[jsol]*kappa[jsol]*zidzdc[isol]+zz[jsol]*dkhdc[jsol][isol])/den
                - zidzdc[jsol]*zidzdcc2[isol]/den
                - zidzdc[isol]*zidzdcc2[jsol]/den;
            }
        }
    }
    
    dkdJ.resize(nsol);
    dkdc.resize(nsol, vector<double>(nsol,0));
    
    for (isol=0; isol<nsol; ++isol) {
        dkdJ[isol] = zz[isol]*dkhdJ[isol]+z[isol]*kappa[isol]*zidzdJ;
        for (jsol=0; jsol<nsol; ++jsol) {
            dkdc[isol][jsol] = zz[isol]*dkhdc[isol][jsol]+z[isol]*kappa[isol]*zidzdc[jsol];
        }
    }
}

//-----------------------------------------------------------------------------
//! Current density
vec3d FEMultiphasicFSI::CurrentDensity(FEMaterialPoint& pt)
{
    int i;
    const int nsol = (int)m_pSolute.size();
    
    vector<vec3d> j(nsol);
    vector<int> z(nsol);
    vec3d Ie(0,0,0);
    for (i=0; i<nsol; ++i) {
        j[i] = SoluteFlux(pt, i);
        z[i] = m_pSolute[i]->ChargeNumber();
        Ie += j[i]*z[i];
    }
    Ie *= m_Fc;
    
    return Ie;
}

//-----------------------------------------------------------------------------
//! actual concentration
double FEMultiphasicFSI::ConcentrationActual(FEMaterialPoint& pt, const int sol)
{
    FEMultiphasicFSIMaterialPoint& spt = *pt.ExtractData<FEMultiphasicFSIMaterialPoint>();
    
    // effective concentration
    double c = spt.m_c[sol];
    
    // partition coefficient
    double kappa = PartitionCoefficient(pt, sol);
    
    double ca = kappa*c;
    
    return ca;
}

//-----------------------------------------------------------------------------
//! actual fluid pressure
double FEMultiphasicFSI::PressureActual(FEMaterialPoint& pt)
{
    int i;
    
    FEFluidMaterialPoint& fpt = *pt.ExtractData<FEFluidMaterialPoint>();
    const int nsol = (int)m_pSolute.size();
    
    // effective pressure
    double p = Fluid()->Pressure(pt);
    
    // actual concentration
    vector<double> c(nsol);
    for (i=0; i<nsol; ++i)
        c[i] = ConcentrationActual(pt, i);
    
    // osmotic coefficient
    double osmc = m_pOsmC->OsmoticCoefficient(pt);
    
    // actual pressure
    double ca = 0;
    for (i=0; i<nsol; ++i) ca += c[i];
    double pa = p + m_Rgas*m_Tabs*osmc*ca;
    
    return pa;
}

//-----------------------------------------------------------------------------
//! Fixed charge density in current configuration
double FEMultiphasicFSI::FixedChargeDensity(FEMaterialPoint& pt)
{
    FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    FEBiphasicFSIMaterialPoint& bt = *pt.ExtractData<FEBiphasicFSIMaterialPoint>();
    FEMultiphasicFSIMaterialPoint& spt = *pt.ExtractData<FEMultiphasicFSIMaterialPoint>();
    
    // relative volume
    double J = et.m_J;
    double phi0 = bt.m_phi0;
    double phif = Porosity(pt);
    
    double cFr = m_cFr(pt);
    double cF = cFr*(1-phi0)/(J - phi0);
    
    return cF;
}

//-----------------------------------------------------------------------------
//! Calculate solute molar flux

vec3d FEMultiphasicFSI::SoluteFlux(FEMaterialPoint& pt, const int sol)
{
    FEMultiphasicFSIMaterialPoint& spt = *pt.ExtractData<FEMultiphasicFSIMaterialPoint>();
    FEFSIMaterialPoint& fpt = *pt.ExtractData<FEFSIMaterialPoint>();
    
    // concentration gradient
    vec3d gradc = spt.m_gradc[sol];
    
    // solute properties
    mat3ds D = Diffusivity(pt, sol);
    double khat = m_pSolute[sol]->m_pSolub->Solubility(pt);
    int z = m_pSolute[sol]->ChargeNumber();
    double zeta = ElectricPotential(pt, true);
    double zz = pow(zeta, z);
    double kappa = zz*khat;
    double d0 = GetSolute(sol)->m_pDiff->Free_Diffusivity(pt);
    
    double c = spt.m_c[sol];
    vec3d w = fpt.m_w;
    double phif = Porosity(pt);
    
    // solute flux j
    vec3d j = D*(-gradc*phif + w*c/d0)*kappa;
    
    return j;
}

//-----------------------------------------------------------------------------
void FEMultiphasicFSI::AddChemicalReaction(FEChemicalReaction* pcr)
{
    m_pReact.push_back(pcr);
}

//-----------------------------------------------------------------------------
double FEMultiphasicFSI::GetReferentialFixedChargeDensity(const FEMaterialPoint& mp)
{
    const FEElasticMaterialPoint* ept = (mp.ExtractData<FEElasticMaterialPoint >());
    const FEMultiphasicFSIMaterialPoint* mfspt = mp.ExtractData<FEMultiphasicFSIMaterialPoint>();
    const FEBiphasicFSIMaterialPoint* bfspt = mp.ExtractData<FEBiphasicFSIMaterialPoint>();
    double cf = (ept->m_J - bfspt->m_phi0) * mfspt->m_cF / (1 - bfspt->m_phi0);
    return cf;
}
