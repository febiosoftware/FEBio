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
#include "FECarterHayes.h"
#include "FEMultiphasic.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FECarterHayes, FEElasticMaterial)
	ADD_PARAMETER(m_E0  , FE_RANGE_GREATER(0.0)         , "E0"   )->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_rho0, FE_RANGE_GREATER(0.0)         , "rho0" )->setUnits(UNIT_DENSITY);
	ADD_PARAMETER(m_g   , FE_RANGE_GREATER_OR_EQUAL(0.0), "gamma");
	ADD_PARAMETER(m_v   , FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v"    );
	ADD_PARAMETER(m_sbm , "sbm");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
bool FECarterHayes::Init()
{
	if (FEElasticMaterial::Init() == false) return false;
	
	// get the parent material which must be a multiphasic material
    FEMultiphasic* pMP = GetAncestor()->ExtractProperty<FEMultiphasic>();
	if (pMP == 0) {
		feLogError("Parent material must be multiphasic");
		return false;
	}

	// extract the local id of the SBM whose density controls Young's modulus from the global id
	m_lsbm = pMP->FindLocalSBMID(m_sbm);
	if (m_lsbm == -1) {
		feLogError("Invalid value for sbm");
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
void FECarterHayes::Serialize(DumpStream& ar)
{
	FEElasticMaterial::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_lsbm;
}

//-----------------------------------------------------------------------------
//! Create material point data
FEMaterialPointData* FECarterHayes::CreateMaterialPointData()
{
	return new FERemodelingMaterialPoint(new FEElasticMaterialPoint);
}

//-----------------------------------------------------------------------------
double FECarterHayes::StrainEnergy(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
	double detF = pt.m_J;
	double lndetF = log(detF);
	
	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();
	double I1 = b.tr();
	
	// lame parameters
	double rhor = spt.m_sbmr[m_lsbm];
	double m_E = YoungModulus(rhor);
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);
	
	double sed = mu*((I1-3)/2 - lndetF)+lam*lndetF*lndetF/2;
	
	return sed;
}

//-----------------------------------------------------------------------------
mat3ds FECarterHayes::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
	double detF = pt.m_J;
	double detFi = 1.0/detF;
	double lndetF = log(detF);

	// evaluate the strain energy
	FERemodelingMaterialPoint& rpt = *mp.ExtractData<FERemodelingMaterialPoint>();
	rpt.m_sed = StrainEnergy(mp);
	
	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();
	
	// lame parameters
	double rhor = spt.m_sbmr[m_lsbm];
	double m_E = YoungModulus(rhor);
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);
	
	// Identity
	mat3dd I(1);
	
	// calculate stress
	mat3ds s = (b - I)*(mu*detFi) + I*(lam*lndetF*detFi);
	
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FECarterHayes::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
	// deformation gradient
	double detF = pt.m_J;
	
	// lame parameters
	double rhor = spt.m_sbmr[m_lsbm];
	double m_E = YoungModulus(rhor);
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
	
	return tens4ds(D);
}

//-----------------------------------------------------------------------------
//! calculate tangent of strain energy density with mass density
double FECarterHayes::Tangent_SE_Density(FEMaterialPoint& mp)
{
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	double rhor = spt.m_sbmr[m_lsbm];
    return StrainEnergy(mp)*m_g/rhor;
}

//-----------------------------------------------------------------------------
//! calculate tangent of stress with mass density
mat3ds FECarterHayes::Tangent_Stress_Density(FEMaterialPoint& mp)
{
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	double rhor = spt.m_sbmr[m_lsbm];
    return Stress(mp)*m_g/rhor;
}

