#include "stdafx.h"
#include "FEHolmesMow.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEHolmesMow, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E");
	ADD_PARAMETER(m_v, FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v");
	ADD_PARAMETER(m_b, FE_RANGE_GREATER_OR_EQUAL(0.0), "beta");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
bool FEHolmesMow::Validate()
{
	if (FEElasticMaterial::Validate() == false) return false;
	
	// Lame coefficients
	lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	mu  = 0.5*m_E/(1+m_v);
	Ha = lam + 2*mu;	

	return true;
}

//-----------------------------------------------------------------------------
mat3ds FEHolmesMow::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double detF = pt.m_J;
	double detFi = 1.0/detF;
	
	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen(); //(F*F.transpose()).sym();
	mat3ds b2 = b*b;
	mat3ds identity(1.,1.,1.,0.,0.,0.);

	// calculate invariants of B
	double I1 = b.tr();
	double I2 = (I1*I1 - b2.tr())/2.;
	double I3 = b.det();

	// Exponential term
	double eQ = exp(m_b*((2*mu-lam)*(I1-3) + lam*(I2-3))/Ha)/pow(I3,m_b);
	
	// calculate stress
	mat3ds s = 0.5*detFi*eQ*((2*mu+lam*(I1-1))*b - lam*b2 - Ha*identity);

	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEHolmesMow::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double detF = pt.m_J;
	double detFi = 1.0/detF;
	
	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen(); //(F*F.transpose()).sym();
	mat3ds b2 = b*b;
	mat3ds identity(1.,1.,1.,0.,0.,0.);
	
	// calculate invariants of B
	double I1 = b.tr();
	double I2 = (I1*I1 - b2.tr())*0.5;
	double I3 = b.det();
	
	// Exponential term
	double eQ = exp(m_b*((2*mu-lam)*(I1-3) + lam*(I2-3))/Ha)/pow(I3,m_b);
	
	// calculate stress
	mat3ds s = 0.5*detFi*eQ*((2*mu+lam*(I1-1))*b - lam*b2 - Ha*identity);
	
	// calculate elasticity tensor
	tens4ds c = 4.*m_b/Ha*detF/eQ*dyad1s(s) 
	+ detFi*eQ*(lam*(dyad1s(b) - dyad4s(b)) + Ha*dyad4s(identity));
	
	return c;
}

//-----------------------------------------------------------------------------
double FEHolmesMow::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen(); //(F*F.transpose()).sym();
	mat3ds b2 = b*b;
    
	// calculate invariants of B
	double I1 = b.tr();
	double I2 = (I1*I1 - b2.tr())/2.;
	double I3 = b.det();
    
	// Exponential term
	double eQ = exp(m_b*((2*mu-lam)*(I1-3) + lam*(I2-3))/Ha)/pow(I3,m_b);
	
	// calculate strain energy density
	double sed = Ha/(4*m_b)*(eQ - 1);
	
	return sed;
}
