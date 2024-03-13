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
#include "FEMultiphasic.h"
#include <FECore/FEModel.h>
#include <FECore/FECoreKernel.h>
#include <FECore/log.h>
#include <FECore/tens4d.h>
#include <FECore/tools.h>
#include <complex>
using namespace std;

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

// Material parameters for the FEMultiphasic material
BEGIN_FECORE_CLASS(FEMultiphasic, FEMaterial)
	ADD_PARAMETER(m_phi0   , FE_RANGE_CLOSED     (0.0, 1.0), "phi0"         );
	ADD_PARAMETER(m_rhoTw  , FE_RANGE_GREATER_OR_EQUAL(0.0), "fluid_density");
	ADD_PARAMETER(m_penalty, FE_RANGE_GREATER_OR_EQUAL(0.0), "penalty"      );
	ADD_PARAMETER(m_cFr    , "fixed_charge_density");

	// define the material properties
	ADD_PROPERTY(m_pSolid , "solid"              , FEProperty::Required | FEProperty::TopLevel);
	ADD_PROPERTY(m_pPerm  , "permeability"       );
	ADD_PROPERTY(m_pOsmC  , "osmotic_coefficient");
	ADD_PROPERTY(m_pSupp  , "solvent_supply"     , FEProperty::Optional);
	ADD_PROPERTY(m_pSolute, "solute"             , FEProperty::Optional);
	ADD_PROPERTY(m_pSBM   , "solid_bound"        , FEProperty::Optional);
	ADD_PROPERTY(m_pReact , "reaction"           , FEProperty::Optional);
    ADD_PROPERTY(m_pMReact, "membrane_reaction"  , FEProperty::Optional);

	ADD_PROPERTY(m_Q, "mat_axis", FEProperty::Optional);

END_FECORE_CLASS();

//=============================================================================
//   FEMultiphasic
//=============================================================================

//-----------------------------------------------------------------------------
//! FEMultiphasic constructor
FEMultiphasic::FEMultiphasic(FEModel* pfem) : FEMaterial(pfem)
{	
	m_rhoTw = 0;
	m_Rgas = 0; m_Tabs = 0; m_Fc = 0;
	m_penalty = 1;

	m_pSolid = 0;
	m_pPerm = 0;
	m_pOsmC = 0;
	m_pSupp = 0;

	m_bool_refC = false;
}

//-----------------------------------------------------------------------------
void FEMultiphasic::AddSolidBoundMolecule(FESolidBoundMolecule* psbm)
{
	m_pSBM.push_back(psbm);
}

//-----------------------------------------------------------------------------
void FEMultiphasic::AddChemicalReaction(FEChemicalReaction* pcr)
{
	m_pReact.push_back(pcr);
}

//-----------------------------------------------------------------------------
void FEMultiphasic::AddMembraneReaction(FEMembraneReaction* pcr)
{
    m_pMReact.push_back(pcr);
}

//-----------------------------------------------------------------------------
//! Returns the local ID of the SBM, given the global ID.
//! \param nid global ID (one - based)
//! \return the local ID (zero-based index) or -1 if not found.
int FEMultiphasic::FindLocalSBMID(int nid)
{
	int lsbm = -1;
	int nsbm = (int) SBMs();
	for (int isbm=0; isbm<nsbm; ++isbm) {
		if (m_pSBM[isbm]->GetSBMID() == nid) {
			lsbm = isbm;
			break;
		}
	}
	return lsbm;
}

//-----------------------------------------------------------------------------
bool FEMultiphasic::Init()
{
	// set the solute IDs first, since they are referenced in FESolute::Init()
	for (int i = 0; i<Solutes(); ++i) {
		m_pSolute[i]->SetSoluteLocalID(i);
	}

    if (m_pSolid->Init() == false) return false;
    if (m_pPerm->Init() == false) return false;
    if (m_pOsmC->Init() == false) return false;
    if (m_pSupp && (m_pSupp->Init() == false)) return false;
    for (int i=0; i<Solutes(); ++i)
        if (m_pSolute[i]->Init() == false) return false;
    for (int i=0; i<SBMs(); ++i)
        if (m_pSBM[i]->Init() == false) return false;
    for (int i=0; i<Reactions(); ++i)
        if (m_pReact[i]->Init() == false) return false;
    for (int i=0; i<MembraneReactions(); ++i)
        if (m_pMReact[i]->Init() == false) return false;

	// call the base class.
	// This also initializes all properties
	if (FEMaterial::Init() == false) return false;

	// Determine how to solve for the electric potential psi
	int isol;
	int zmin = 0, zmax = 0, z;
	for (isol=0; isol<(int)m_pSolute.size(); ++isol) {
		z = m_pSolute[isol]->ChargeNumber();
		if (z < zmin) zmin = z;
		if (z > zmax) zmax = z;
	}
	m_zmin = zmin;
	m_ndeg = zmax - zmin;	// polynomial degree

	m_Rgas = GetFEModel()->GetGlobalConstant("R");
	m_Tabs = GetFEModel()->GetGlobalConstant("T");
	m_Fc   = GetFEModel()->GetGlobalConstant("Fc");
	m_bool_refC = GetFEModel()->GetGlobalConstant("referential_concentration");


	if (m_Rgas <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section"); return false; }
	if (m_Tabs <= 0) { feLogError("A positive absolute temperature T must be defined in Globals section");	 return false; }
	if ((zmin || zmax) && (m_Fc <= 0)) {
		feLogError("A positive Faraday constant Fc must be defined in Globals section");
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
// update specialized material points
void FEMultiphasic::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp)
{
    m_pSolid->UpdateSpecializedMaterialPoints(mp, tp);
    m_pPerm->UpdateSpecializedMaterialPoints(mp, tp);
    m_pOsmC->UpdateSpecializedMaterialPoints(mp, tp);
    if (m_pSupp) m_pSupp->UpdateSpecializedMaterialPoints(mp, tp);
    for (int i=0; i<Solutes(); ++i)
        m_pSolute[i]->UpdateSpecializedMaterialPoints(mp, tp);
    for (int i=0; i<SBMs(); ++i)
        m_pSBM[i]->UpdateSpecializedMaterialPoints(mp, tp);
    for (int i=0; i<Reactions(); ++i)
        m_pReact[i]->UpdateSpecializedMaterialPoints(mp, tp);
    for (int i=0; i<MembraneReactions(); ++i)
        m_pMReact[i]->UpdateSpecializedMaterialPoints(mp, tp);
}

//-----------------------------------------------------------------------------
void FEMultiphasic::Serialize(DumpStream& ar)
{
	FEMaterial::Serialize(ar);
	if (ar.IsShallow()) return;

	ar & m_Rgas & m_Tabs & m_Fc;
	ar & m_zmin & m_ndeg;
}

//-----------------------------------------------------------------------------
//! Solid referential apparent density
double FEMultiphasic::SolidReferentialApparentDensity(FEMaterialPoint& pt)
{
	FEBiphasicMaterialPoint& pet = *pt.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
		
	// evaluate referential apparent density of base solid
	double density = m_pSolid->Density(pt);
	double rhosr = pet.m_phi0t*density;

	// add contribution from solid-bound molecules
	for (int isbm=0; isbm<(int)spt.m_sbmr.size(); ++isbm)
		rhosr += spt.m_sbmr[isbm];
	
	return rhosr;
}

//! Return solid referential apparent density
double FEMultiphasic::GetReferentialSolidVolumeFraction(const FEMaterialPoint& pt)
{
    const FEBiphasicMaterialPoint& bt = *pt.ExtractData<FEBiphasicMaterialPoint>();
    return bt.m_phi0t;
}

//-----------------------------------------------------------------------------
//! Evaluate and return solid referential volume fraction
double FEMultiphasic::SolidReferentialVolumeFraction(FEMaterialPoint& pt)
{
	// get referential apparent density of base solid (assumed constant)
	double phisr = m_phi0(pt);
    
    // add contribution from solid-bound 'solutes'
    FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    const int nsol = (int)m_pSolute.size();
    double f = 0;
    for (int isol=0; isol<nsol; ++isol)
        if (spt.m_bsb[isol]) f += spt.m_ca[isol]*m_pSolute[isol]->MolarMass()/m_pSolute[isol]->Density();
    phisr = (phisr + et.m_J*f)/(1+f);
    
	// add contribution from solid-bound molecules
	for (int isbm=0; isbm<(int)m_pSBM.size(); ++isbm)
		phisr += SBMReferentialVolumeFraction(pt, isbm);
    
    FEBiphasicMaterialPoint& bt = *pt.ExtractData<FEBiphasicMaterialPoint>();
    bt.m_phi0t = phisr;
    
	return phisr;
}

//-----------------------------------------------------------------------------
//! Evaluate and return tangent of solid referential volume fraction w.r.t. relative volume
double FEMultiphasic::TangentSRVFStrain(FEMaterialPoint& pt)
{
    // get referential apparent density of base solid (assumed constant)
    double phis0 = m_phi0(pt);
    double phisr = SolidReferentialVolumeFraction(pt);
    
    // add contribution from solid-bound 'solutes'
    FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    double J = et.m_J;
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    const int nsol = (int)m_pSolute.size();
    double f = 0;
    double s = 0;
    for (int isol=0; isol<nsol; ++isol) {
        if (spt.m_bsb[isol]) {
            f += spt.m_ca[isol]*m_pSolute[isol]->MolarMass()/m_pSolute[isol]->Density();
            s += spt.m_ca[isol]*m_pSolute[isol]->MolarMass()/m_pSolute[isol]->Density()*spt.m_dkdJ[isol]/spt.m_k[isol];
        }
    }
    double d = 1+f;
    double dphisrdJ = f/d + (J-phis0)*s/(d*d);
    
    // add contribution from solid-bound molecules
    f = s = 0;
    for (int isbm=0; isbm<(int)m_pSBM.size(); ++isbm) {
        f += spt.m_sbmr[isbm]/m_pSBM[isbm]->Density();
    }
    dphisrdJ += f/(J-phisr + f);
    
    return dphisrdJ;
}

//-----------------------------------------------------------------------------
//! evaluate and return tangent of  solid referential volume fraction w.r.t. to concentration
double FEMultiphasic::TangentSRVFConcentration(FEMaterialPoint& pt, const int sol)
{
    // get referential apparent density of base solid (assumed constant)
    double phis0 = m_phi0(pt);
    double phisr = SolidReferentialVolumeFraction(pt);
    
    // add contribution from solid-bound 'solutes'
    FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    double J = et.m_J;
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    const int nsol = (int)m_pSolute.size();
    double f = 0;
    if (spt.m_bsb[sol]) {
        for (int isol=0; isol<nsol; ++isol) {
            if (spt.m_bsb[isol]) {
                f += spt.m_ca[isol]*spt.m_dkdc[isol][sol]*m_pSolute[isol]->MolarMass()/m_pSolute[isol]->Density();
                if (isol == sol) f += spt.m_k[isol]*m_pSolute[isol]->MolarMass()/m_pSolute[isol]->Density();
            }
        }
    }
    double dphisrdc = pow(J-phisr,2)/(J-phis0)*f;
    
    return dphisrdc;
}


//-----------------------------------------------------------------------------
//! Porosity in current configuration
double FEMultiphasic::Porosity(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& bt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	
	// solid referential volume fraction
	double phisr = bt.m_phi0t;
    
	// relative volume
	double J = et.m_J;
	double J_eval = (m_bool_refC) ? 1.0 : J;
	
	double phiw = 1 - phisr/J_eval;
	// check for pore collapse
	// TODO: throw an error if pores collapse
	phiw = (phiw > 0) ? phiw : 0;
	
	return phiw;
}

//-----------------------------------------------------------------------------
//! Fixed charge density in current configuration
double FEMultiphasic::FixedChargeDensity(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& bt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	
	// relative volume
	double J = et.m_J;
	double phi0 = bt.m_phi0t;
	double ce = 0;

	// add contribution from charged solid-bound molecules
	for (int isbm=0; isbm<(int)m_pSBM.size(); ++isbm)
		ce += SBMChargeNumber(isbm)*spt.m_sbmr[isbm]/SBMMolarMass(isbm);
    
    double cFr = m_cFr(pt);
	double cF = (cFr*(1-bt.m_phi0)+ce)/(J-phi0);
    
    // add contribution from solid-bound 'solutes'
    const int nsol = (int)m_pSolute.size();
    for (int isol=0; isol<nsol; ++isol)
        if (spt.m_bsb[isol]) cF += spt.m_ca[isol]*m_pSolute[isol]->ChargeNumber();

	return cF;
}

//-----------------------------------------------------------------------------
//! Electric potential
double FEMultiphasic::ElectricPotential(FEMaterialPoint& pt, const bool eform)
{
	// check if solution is neutral
	if (m_ndeg == 0) {
		if (eform) return 1.0;
		else return 0.0;
	}
	
	int i, j;
	
	// if not neutral, solve electroneutrality polynomial for zeta
	FESolutesMaterialPoint& set = *pt.ExtractData<FESolutesMaterialPoint>();
	const int nsol = (int)m_pSolute.size();
	double cF = FixedChargeDensity(pt);

	vector<double> c(nsol,0);       // effective concentration
	vector<double> khat(nsol,1);	// solubility
	vector<int> z(nsol,0);          // charge number
	for (i=0; i<nsol; ++i) {
        if (!set.m_bsb[i]) {
            c[i] = set.m_c[i];
            khat[i] = m_pSolute[i]->m_pSolub->Solubility(pt);
            z[i] = m_pSolute[i]->ChargeNumber();
        }
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
	double psi = set.m_psi;		// use previous solution as initial guess
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
double FEMultiphasic::PartitionCoefficient(FEMaterialPoint& pt, const int sol)
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
void FEMultiphasic::PartitionCoefficientFunctions(FEMaterialPoint& mp, vector<double>& kappa,
                                                  vector<double>& dkdJ,
                                                  vector< vector<double> >& dkdc,
                                                  vector< vector<double> >& dkdr,
                                                  vector< vector<double> >& dkdJr,
                                                  vector< vector< vector<double> > >& dkdrc)
{
	int isol, jsol, ksol;
    int isbm;
	
	FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
	FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
	FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
	
	const int nsol = (int)m_pSolute.size();
    const int nsbm = (int)m_pSBM.size();
	
	vector<double> c(nsol);
	vector<int> z(nsol);
	vector<double> khat(nsol);
	vector<double> dkhdJ(nsol);
	vector<double> dkhdJJ(nsol);
	vector< vector<double> > dkhdc(nsol, vector<double>(nsol));
	vector< vector<double> > dkhdJc(nsol, vector<double>(nsol));
	vector< vector< vector<double> > > dkhdcc(nsol, dkhdc);	// use dkhdc to initialize only
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
	double phi0 = ppt.m_phi0t;
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
    vector<double> zidzdr(nsbm,0);
    vector<double> zidzdJr(nsbm,0);
	vector< vector<double> > zidzdrc(nsbm, vector<double>(nsol,0));
	dkdr.resize(nsol, vector<double>(nsbm));
	dkdJr.resize(nsol, vector<double>(nsbm));
	dkdrc.resize(nsol, zidzdrc);	// use zidzdrc for initialization only
    
	if (den > 0) {
		
        for (isbm=0; isbm<nsbm; ++isbm) {
            zidzdr[isbm] = -(cF/SBMDensity(isbm) + SBMChargeNumber(isbm)/SBMMolarMass(isbm))/(J-phi0)/den;
            
            for (isol=0; isol<nsol; ++isol) {
                zidzdJr[isbm] += SQR(z[isol])*dkdJ[isol]*c[isol];
            }
            zidzdJr[isbm] = 1/(J-phi0) + zidzdJr[isbm]/den;
            zidzdJr[isbm] = (zidzdJr[isbm] + zidzdJ)*zidzdr[isbm];
            zidzdJr[isbm] += cF/SBMDensity(isbm)/SQR(J-phi0)/den;
            
            for (isol=0; isol<nsol; ++isol) {
                zidzdrc[isbm][isol] = SQR(z[isol])*kappa[isol];
                for (jsol=0; jsol<nsol; ++jsol)
                    zidzdrc[isbm][isol] += SQR(z[jsol])*zz[jsol]*c[jsol]*dkhdc[jsol][isol];
                zidzdrc[isbm][isol] = zidzdr[isbm]*(zidzdc[isol]*(1+num/den) - zidzdrc[isbm][isol]/den);
            }
        }
	}
    
    for (isbm=0; isbm<nsbm; ++isbm) {
        for (isol=0; isol<nsol; ++isol) {
            dkdr[isol][isbm] = z[isol]*kappa[isol]*zidzdr[isbm];
            dkdJr[isol][isbm] = z[isol]*((dkhdJ[isol]*zz[isol]+(z[isol]-1)*kappa[isol]*zidzdJ)*zidzdr[isbm]
                                         +kappa[isol]*zidzdJr[isbm]);
            for (jsol=0; jsol<nsol; ++jsol) {
                dkdrc[isol][isbm][jsol] = z[isol]*(zz[isol]*dkhdc[isol][jsol]*zidzdr[isbm]
                                                   +(z[isol]-1)*kappa[isol]*zidzdr[isbm]*zidzdc[jsol]
                                                   +kappa[isol]*zidzdrc[isbm][jsol]);
            }
        }
    }
}

//-----------------------------------------------------------------------------
//! actual concentration
double FEMultiphasic::Concentration(FEMaterialPoint& pt, const int sol)
{
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	
	// effective concentration
	double c = spt.m_c[sol];
	
	// partition coefficient
	double kappa = PartitionCoefficient(pt, sol);
	
	// actual concentration
	double ca = kappa*c;
	
	return ca;
}

//-----------------------------------------------------------------------------
//! The stress of a triphasic material is the sum of the fluid pressure
//! and the elastic stress. Note that this function is declared in the base class
//! so you do not have to reimplement it in a derived class, unless additional
//! pressure terms are required.

mat3ds FEMultiphasic::Stress(FEMaterialPoint& mp)
{
	// calculate solid material stress
	mat3ds s = m_pSolid->Stress(mp);
	
	// fluid pressure
	double p = Pressure(mp);
	
	// add fluid pressure
	s.xx() -= p;
	s.yy() -= p;
	s.zz() -= p;
	
	return s;
}

//-----------------------------------------------------------------------------
//! The tangent is the elastic tangent. Note
//! that this function is declared in the base class, so you don't have to 
//! reimplement it unless additional tangent components are required.

tens4ds FEMultiphasic::Tangent(FEMaterialPoint& mp)
{
	int i;
	
	FEElasticMaterialPoint& ept = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& bpt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	const int nsol = (int)m_pSolute.size();
	
	// call solid tangent routine
	tens4ds C = m_pSolid->Tangent(mp);
	double D[6][6] = {0};
	C.extract(D);
	
	// relative volume and solid volume fraction
	double J = ept.m_J;
	double phi0 = bpt.m_phi0t;
	
	// get the charge density and its derivatives
	double cF = FixedChargeDensity(mp);
	double dcFdJ = -cF/(J - phi0);
	
	// fluid pressure and solute concentration
	double p = Pressure(mp);
	
	// get remaining variables
	double zeta = ElectricPotential(mp, true);
	vector<double> c(nsol);
	vector<int> z(nsol);
	vector<double> khat(nsol);
	vector<double> dkhdJ(nsol);
	vector<double> zz(nsol);
	vector<double> kappa(nsol);
	double den = 0;
	for (i=0; i<nsol; ++i) {
		c[i] = spt.m_c[i];
		z[i] = m_pSolute[i]->ChargeNumber();
		khat[i] = m_pSolute[i]->m_pSolub->Solubility(mp);
		dkhdJ[i] = m_pSolute[i]->m_pSolub->Tangent_Solubility_Strain(mp);
		zz[i] = pow(zeta, z[i]);
		kappa[i] = zz[i]*khat[i];
		den += SQR(z[i])*kappa[i]*c[i];
	}
	
	// evaluate electric potential (nondimensional exponential form) and its derivatives
	// also evaluate partition coefficients and their derivatives
	double zidzdJ = 0;
	if (den > 0) {
		zidzdJ = dcFdJ;
		for (i=0; i<nsol; ++i)
			zidzdJ += z[i]*zz[i]*dkhdJ[i]*c[i];
		zidzdJ = -zidzdJ/den;
	}
	vector<double> dkdJ(nsol);
	for (i=0; i<nsol; ++i) dkdJ[i] = zz[i]*dkhdJ[i]+z[i]*kappa[i]*zidzdJ;
	
	// osmotic coefficient and its derivative w.r.t. strain
	double osmc = m_pOsmC->OsmoticCoefficient(mp);
	double dodJ = m_pOsmC->Tangent_OsmoticCoefficient_Strain(mp);
	
	double dp = 0;
	for (i=0; i<nsol; ++i) {
		dp += c[i]*(osmc*dkdJ[i]+dodJ*kappa[i]);
	}
	dp *= m_Rgas*m_Tabs*J;
	
	// adjust tangent for pressures
	D[0][0] -= -p + dp;
	D[1][1] -= -p + dp;
	D[2][2] -= -p + dp;
	
	D[0][1] -= p + dp; D[1][0] -= p + dp;
	D[1][2] -= p + dp; D[2][1] -= p + dp;
	D[0][2] -= p + dp; D[2][0] -= p + dp;
	
	D[3][3] -= -p;
	D[4][4] -= -p;
	D[5][5] -= -p;
	
	return tens4ds(D);
}

//-----------------------------------------------------------------------------
//! Calculate fluid flux

vec3d FEMultiphasic::FluidFlux(FEMaterialPoint& pt)
{
	int i;
	
	FEBiphasicMaterialPoint& ppt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	const int nsol = (int)m_pSolute.size();
	vector<double> c(nsol);
	vector<vec3d> gradc(nsol);
	vector<mat3ds> D(nsol);
	vector<double> D0(nsol);
	vector<double> khat(nsol);
	vector<int> z(nsol);
	vector<double> zz(nsol);
	vector<double> kappa(nsol);
	
	// fluid volume fraction (porosity) in current configuration
	double phiw = Porosity(pt);
	
	// pressure gradient
	vec3d gradp = ppt.m_gradp;
	
	// electric potential
	double zeta = ElectricPotential(pt, true);
	
	
	for (i=0; i<nsol; ++i) {
		// concentration
		c[i] = spt.m_c[i];

		// concentration gradient
		gradc[i] = spt.m_gradc[i];
		
		// solute diffusivity in mixture
		D[i] = m_pSolute[i]->m_pDiff->Diffusivity(pt);
		
		// solute free diffusivity
		D0[i] = m_pSolute[i]->m_pDiff->Free_Diffusivity(pt);
		
		// solubility
		khat[i] = m_pSolute[i]->m_pSolub->Solubility(pt);
		z[i] = m_pSolute[i]->ChargeNumber();
		zz[i] = pow(zeta, z[i]);
		kappa[i] = zz[i]*khat[i];
	}
	
	// identity matrix
	mat3dd I(1);
	
	// hydraulic permeability
	mat3ds kt = m_pPerm->Permeability(pt);
	
	// effective hydraulic permeability
	mat3ds ke;
	ke.zero();
	for (i=0; i<nsol; ++i)
		ke += (I-D[i]/D0[i])*(kappa[i]*c[i]/D0[i]);
	ke = (kt.inverse() + ke*(m_Rgas*m_Tabs/phiw)).inverse();
	
	// fluid flux w
	vec3d w(0,0,0);
	for (i=0; i<nsol; ++i)
		w += (D[i]*gradc[i])*(kappa[i]/D0[i]);
	w = -(ke*(gradp + w*m_Rgas*m_Tabs));
	
	return w;
}

//-----------------------------------------------------------------------------
//! Calculate solute molar flux

vec3d FEMultiphasic::SoluteFlux(FEMaterialPoint& pt, const int sol)
{
	FEBiphasicMaterialPoint& bpt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	
	// fluid volume fraction (porosity) in current configuration
	double phiw = Porosity(pt);
	
	// concentration
	double c = spt.m_c[sol];
	
	// concentration gradient
	vec3d gradc = spt.m_gradc[sol];
	
	// solute diffusivity in mixture
	mat3ds D = m_pSolute[sol]->m_pDiff->Diffusivity(pt);
	
	// solute free diffusivity
	double D0 = m_pSolute[sol]->m_pDiff->Free_Diffusivity(pt);
	
	// solubility
	double khat = m_pSolute[sol]->m_pSolub->Solubility(pt);
	int z = m_pSolute[sol]->ChargeNumber();
	double zeta = ElectricPotential(pt, true);
	double zz = pow(zeta, z);
	double kappa = zz*khat;
	
	// fluid flux w
	vec3d w = bpt.m_w;
	
	// solute flux j
	vec3d j = (D*(w*(c/D0) - gradc*phiw))*kappa;
	
	return j;
}

//-----------------------------------------------------------------------------
//! actual fluid pressure
double FEMultiphasic::Pressure(FEMaterialPoint& pt)
{
	FEBiphasicMaterialPoint& ppt = *pt.ExtractData<FEBiphasicMaterialPoint>();
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	const int nsol = (int)m_pSolute.size();
	
	// effective pressure
	double p = ppt.m_p;
	
	// osmolarity
    double c = spt.Osmolarity();
	
	// osmotic coefficient
	double osmc = m_pOsmC->OsmoticCoefficient(pt);
	
	// actual pressure
	double pa = p + m_Rgas*m_Tabs*osmc*c;
	
	return pa;
}

//-----------------------------------------------------------------------------
//! Current density
vec3d FEMultiphasic::CurrentDensity(FEMaterialPoint& pt)
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
//! Evaluate effective permeability
mat3ds FEMultiphasic::EffectivePermeability(FEMaterialPoint& pt)
{
    // evaluate the hydraulic permeability
    mat3ds K = GetPermeability()->Permeability(pt);

    const int nsol = Solutes();
    
    // if there are no solutes in this mixture, we're done
    if (nsol == 0) return K;

    // initialize effective permeability
    mat3ds Ke = K.inverse();
    
    // fluid volume fraction (porosity) in current configuration
    double phiw = Porosity(pt);
    double tmp = m_Rgas*m_Tabs/phiw;
    mat3dd I(1.0);

    FESolutesMaterialPoint&  spt = *(pt.ExtractData<FESolutesMaterialPoint >());

    // add solute contributions (but not 'solid-bound' solutes)
    for (int isol=0; isol<nsol; ++isol) {
        if (!spt.m_bsb[isol]) {
            // concentration
            double ca = spt.m_ca[isol];
            // solute diffusivity in mixture
            mat3ds D = m_pSolute[isol]->m_pDiff->Diffusivity(pt);
            // solute free diffusivity
            double D0 = m_pSolute[isol]->m_pDiff->Free_Diffusivity(pt);
            
            Ke += (I - D/D0)*(tmp*ca/D0);
        }
    }
        
    return Ke.inverse();
}

//-----------------------------------------------------------------------------
//! Evaluate tangent of effective permeability w.r.t. strain
tens4dmm FEMultiphasic::TangentPermeabilityStrain(FEMaterialPoint& pt, const mat3ds& Ke)
{
	if (m_bool_refC)
		return tens4dmm(0.0);

	// get the hydraulic permeability strain tangent
    tens4dmm dKdE = GetPermeability()->Tangent_Permeability_Strain(pt);
    
    const int nsol = Solutes();
    
    // if there are no solutes in this mixture, we're done
    if (nsol == 0) return dKdE;
    
    // evaluate the inverse of the hydraulic permeability
    mat3ds Ki = GetPermeability()->Permeability(pt).inverse();
    
    // fluid volume fraction (porosity) in current configuration
    double phiw = Porosity(pt);
    mat3dd I(1.0);
    
    FEElasticMaterialPoint&  ept = *(pt.ExtractData<FEElasticMaterialPoint >());
    FESolutesMaterialPoint&  spt = *(pt.ExtractData<FESolutesMaterialPoint >());
    
    tens4dmm dKedE;
    dKedE.zero();
    
    // add solute contributions
    for (int isol=0; isol<nsol; ++isol) {
        // concentration
        double ca = spt.m_ca[isol];
        // solute free diffusivity
        double D0 = m_pSolute[isol]->m_pDiff->Free_Diffusivity(pt);
        // solute diffusivity in mixture, normalized by D0
        mat3ds D = m_pSolute[isol]->m_pDiff->Diffusivity(pt)/D0;
        // solute diffusiviety strain tangent, normalized by D0
        tens4dmm dDdE = m_pSolute[isol]->m_pDiff->Tangent_Diffusivity_Strain(pt)/D0;
        
        dKedE += (dDdE + (dyad4mm(I, D) + dyad4mm(D,I) - dyad4mm(I,I))*2
        - dyad1mm(D,I) + dyad1mm(I-D,I)*(1./phiw - ept.m_J*spt.m_dkdJ[isol]/spt.m_k[isol]))*ca/D0;
    }
    
    dKedE = dKedE*(m_Rgas*m_Tabs/phiw) + ddot(dyad2mm(Ki,Ki), dKdE);
    
    return ddot(dyad2mm(Ke, Ke), dKedE);
}

//-----------------------------------------------------------------------------
//! Evaluate tangent of effective permeability w.r.t. concentration
mat3ds FEMultiphasic::TangentPermeabilityConcentration(FEMaterialPoint& pt, const int sol, const mat3ds& Ke)
{
    mat3ds dKedc(0,0,0,0,0,0);

	if (m_bool_refC)
		return dKedc;
    
    const int nsol = Solutes();
    
    if (nsol == 0) return dKedc;
    
    // fluid volume fraction (porosity) in current configuration
    double phiw = Porosity(pt);
    mat3dd I(1.0);
    
    FESolutesMaterialPoint&  spt = *(pt.ExtractData<FESolutesMaterialPoint >());
    
    // add solute contributions
    for (int isol=0; isol<nsol; ++isol) {
        // concentration
        double ca = spt.m_ca[isol];
        // solute free diffusivity
        double D0 = m_pSolute[isol]->m_pDiff->Free_Diffusivity(pt);
        // solute diffusivity in mixture, normalized by D0
        mat3ds D = m_pSolute[isol]->m_pDiff->Diffusivity(pt)/D0;
        // solute free diffusivity concentration tangent
        double dD0dc = m_pSolute[isol]->m_pDiff->Tangent_Free_Diffusivity_Concentration(pt, sol);
        // solute diffusivity concentration tangent
        mat3ds dDdc = m_pSolute[isol]->m_pDiff->Tangent_Diffusivity_Concentration(pt,sol);
        
        double kd = (isol == sol) ? 1 : 0;
        dKedc += (I-D)*((spt.m_dkdc[isol][sol]*spt.m_c[isol] + spt.m_k[isol]*kd - dD0dc*ca/D0)/D0)
        - (dDdc - D*dD0dc)*(ca/D0/D0);
    }
    
    dKedc *= m_Rgas*m_Tabs/phiw;
    
    return -(Ke*dKedc*Ke).sym();
}

double FEMultiphasic::GetReferentialFixedChargeDensity(const FEMaterialPoint& mp)
{
	const FEElasticMaterialPoint* ept = (mp.ExtractData<FEElasticMaterialPoint >());
	const FEBiphasicMaterialPoint* bpt = (mp.ExtractData<FEBiphasicMaterialPoint>());
	const FESolutesMaterialPoint* spt = (mp.ExtractData<FESolutesMaterialPoint >());
	double cf = (ept->m_J - bpt->m_phi0t) * spt->m_cF / (1 - bpt->m_phi0);
	return cf;
}
