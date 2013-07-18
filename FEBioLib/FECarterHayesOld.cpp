#include "stdafx.h"
#include "FECarterHayesOld.h"

// register the material with the framework
REGISTER_MATERIAL(FECarterHayesOld, "Carter-Hayes (old)");

// define the material parameters
BEGIN_PARAMETER_LIST(FECarterHayesOld, FEElasticMaterial)
ADD_PARAMETER(m_c, FE_PARAM_DOUBLE, "c");
ADD_PARAMETER(m_g, FE_PARAM_DOUBLE, "gamma");
ADD_PARAMETER(m_v, FE_PARAM_DOUBLE, "v");
END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// FECarterHayes
//////////////////////////////////////////////////////////////////////

void FECarterHayesOld::Init()
{
	FEElasticMaterial::Init();
	
	if (m_c <= 0) throw MaterialError("Invalid value for c");
	if (m_g < 0) throw MaterialError("Invalid value for gamma");
	if (!IN_RIGHT_OPEN_RANGE(m_v, -1.0, 0.5)) throw MaterialRangeError("v", -1.0, 0.5, true, false);
}

double FECarterHayesOld::StrainEnergy(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	double detF = pt.J;
	double lndetF = log(detF);
	
	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();
	double I1 = b.tr();
	
	// lame parameters
	double rhor = pt.rhor;
	double m_E = YoungModulus(rhor);
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);
	
	double sed = mu*((I1-3)/2 - lndetF)+lam*lndetF*lndetF/2;
	
	return sed;
}

mat3ds FECarterHayesOld::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	double detF = pt.J;
	double detFi = 1.0/detF;
	double lndetF = log(detF);
	
	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();
	
	// lame parameters
	double rhor = pt.rhor;
	double m_E = YoungModulus(rhor);
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);
	
	// Identity
	mat3dd I(1);
	
	// calculate stress
	mat3ds s = (b - I)*(mu*detFi) + I*(lam*lndetF*detFi);
	
	return s;
}

tens4ds FECarterHayesOld::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	double detF = pt.J;
	
	// lame parameters
	double rhor = pt.rhor;
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

//! calculate tangent of strain energy density with mass density
double FECarterHayesOld::Tangent_SE_Density(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    return StrainEnergy(mp)*m_g/pt.rhor;
}

//! calculate tangent of stress with mass density
mat3ds FECarterHayesOld::Tangent_Stress_Density(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    return Stress(mp)*m_g/pt.rhor;
}

