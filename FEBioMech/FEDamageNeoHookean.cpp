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
#include "FEDamageNeoHookean.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEDamageNeoHookean, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E");
	ADD_PARAMETER(m_v, FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v");
	ADD_PARAMETER(m_alpha, FE_RANGE_GREATER_OR_EQUAL(0.0), "a");
	ADD_PARAMETER(m_beta , FE_RANGE_CLOSED(0.0, 1.0), "b");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// Constructor
FEDamageNeoHookean::FEDamageNeoHookean(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_E = 0;
	m_v = 0;

	m_alpha = 0.014;
	m_beta = 0.34;
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPointData* FEDamageNeoHookean::CreateMaterialPointData()
{
	return new FEDamageMaterialPoint(new FEElasticMaterialPoint);
}

//-----------------------------------------------------------------------------
// Initialization routine and parameter checking
bool FEDamageNeoHookean::Init()
{
	if (FEElasticMaterial::Init() == false) return false;

	// calculate Lame parameters
	m_lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	m_mu  = 0.5*m_E/(1+m_v);

	return true;
}

//-----------------------------------------------------------------------------
//! Calculate the stress. This happens in two phases. First, we calculate 
//! the stress for the undamaged material. Second, we update the damage
//! parameter and correct the stress accordingly.
mat3ds FEDamageNeoHookean::Stress(FEMaterialPoint& mp)
{
	// --- A. Calculate neo-Hookean stress ----
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double detF = pt.m_J;
	double detFi = 1.0/detF;
	double lndetF = log(detF);

	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();

	// Identity
	mat3dd I(1);

	// calculate stress
	mat3ds s = (b - I)*(m_mu*detFi) + I*(m_lam*lndetF*detFi);

	// --- B. Calculate the damage reduction factor ---
	double g = Damage(mp);

	return s*g;
}

//-----------------------------------------------------------------------------
// Calculate damage reduction factor 
double FEDamageNeoHookean::Damage(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate right Cauchy-Green tensor
	mat3ds C = pt.RightCauchyGreen();

	// Invariants
	double I1 = C.tr();
	double J = pt.m_J;

	// strain-energy value
	double lnJ = log(J);
	double SEF = 0.5*m_mu*(I1 - 3) - m_mu*lnJ + 0.5*m_lam*(lnJ*lnJ);

	// get the damage material point data
	FEDamageMaterialPoint& dp = *mp.ExtractData<FEDamageMaterialPoint>();

	// calculate trial-damage parameter
	dp.m_Etrial = sqrt(2.0*fabs(SEF));

	// calculate damage parameter
	double Es = max(dp.m_Etrial, dp.m_Emax);

	// calculate reduction parameter
	double g = 1.0;
	if (fabs(Es) > 1e-12) g = m_beta + (1.0 - m_beta)*(1.0 - exp(-Es/m_alpha))/(Es/m_alpha);
	else g = 1.0 - 0.5*(1.0 - m_beta)/m_alpha*Es;

	dp.m_D = 1-g;
	return g;
}

//-----------------------------------------------------------------------------
// Calculate tangent. I'm not sure if the tangent needs to be modified for the 
// damage model For now, I don't modify it.
tens4ds FEDamageNeoHookean::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	double detF = pt.m_J;

	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);

	double lam1 = lam / detF;
	double mu1  = (mu - lam*log(detF)) / detF;
	
	double D[6][6] = {0};
	D[0][0] = lam1+2.*mu1; D[0][1] = lam1       ; D[0][2] = lam1       ;
	D[1][0] = lam1       ; D[1][1] = lam1+2.*mu1; D[1][2] = lam1       ;
	D[2][0] = lam1       ; D[2][1] = lam1       ; D[2][2] = lam1+2.*mu1;
	D[3][3] = mu1;
	D[4][4] = mu1;
	D[5][5] = mu1;

	double g = Damage(mp);

	return tens4ds(D)*g;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FEDamageNeoHookean::StrainEnergyDensity(FEMaterialPoint& mp)
{
	// --- A. Calculate neo-Hookean strain energy density ----
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	double J = pt.m_J;
    double lnJ = log(J);
    
	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();
    
    double I1 = b.tr();
    
    double sed = m_mu*(0.5*(I1-3) - lnJ) + 0.5*m_lam*lnJ*lnJ;
    
	// --- B. Calculate the damage reduction factor ---
	double g = Damage(mp);
    
    return sed*g;
}
