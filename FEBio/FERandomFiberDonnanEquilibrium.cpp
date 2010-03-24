#include "stdafx.h"
#include "FERandomFiberDonnanEquilibrium.h"

// The following file contains the integration points and weights
// for the integration over a unit sphere in spherical coordinates
#include "geodesic.h"

// we store the cos and sin of the angles here
int FERandomFiberDonnanEquilibrium::m_nres = 0;
double FERandomFiberDonnanEquilibrium::m_cth[NSTH];
double FERandomFiberDonnanEquilibrium::m_sth[NSTH];
double FERandomFiberDonnanEquilibrium::m_cph[NSTH];
double FERandomFiberDonnanEquilibrium::m_sph[NSTH];
double FERandomFiberDonnanEquilibrium::m_w[NSTH];

// register the material with the framework
REGISTER_MATERIAL(FERandomFiberDonnanEquilibrium, "EFD Donnan equilibrium");

// define the material parameters
BEGIN_PARAMETER_LIST(FERandomFiberDonnanEquilibrium, FEElasticMaterial)
ADD_PARAMETER(m_phiwr, FE_PARAM_DOUBLE, "phiw0");
ADD_PARAMETER(m_cFr, FE_PARAM_DOUBLE, "cF0");
ADD_PARAMETER(m_bosm, FE_PARAM_DOUBLE, "bosm");
ADD_PARAMETER(m_Rgas, FE_PARAM_DOUBLE, "R");
ADD_PARAMETER(m_Tabs, FE_PARAM_DOUBLE, "T");
ADD_PARAMETERV(m_beta, FE_PARAM_DOUBLEV, 3, "beta");
ADD_PARAMETERV(m_ksi , FE_PARAM_DOUBLEV, 3, "ksi" );
END_PARAMETER_LIST();

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//////////////////////////////////////////////////////////////////////
// FERandomFiberDonnanEquilibrium
//////////////////////////////////////////////////////////////////////

FERandomFiberDonnanEquilibrium::FERandomFiberDonnanEquilibrium()
{
	static bool bfirst = true;
	
	if (bfirst)
	{
		// select the integration rule
		const int nint = (m_nres == 0? NSTL  : NSTH  );
		const double* phi = (m_nres == 0? PHIL  : PHIH  );
		const double* the = (m_nres == 0? THETAL: THETAH);
		const double* w   = (m_nres == 0? AREAL : AREAH );
		
		for (int n=0; n<nint; ++n)
		{
			m_cth[n] = cos(the[n]);
			m_sth[n] = sin(the[n]);
			m_cph[n] = cos(phi[n]);
			m_sph[n] = sin(phi[n]);
			m_w[n] = w[n];
		}
		
		bfirst = false;
	}
	
}

void FERandomFiberDonnanEquilibrium::Init()
{
	FEElasticMaterial::Init();
	
	FEDonnanEquilibriumInit(m_phiwr, m_cFr, m_Rgas, m_Tabs, m_bosm);
	FEEllipsoidalFiberDistributionInit(m_ksi, m_beta);
}

mat3ds FERandomFiberDonnanEquilibrium::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	const int nint = (m_nres == 0? NSTL  : NSTH  );
	
	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
			
	// --- M A T R I X   C O N T R I B U T I O N ---

	mat3ds s = FEDonnanEquilibriumStress(m_phiwr, m_cFr, m_Rgas, m_Tabs, m_bosm, J);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// calculate right Cauchy-Green tensor
	mat3ds C = pt.RightCauchyGreen();
	
	// get the element's local coordinate system
	mat3d& Q = pt.Q;
	
	// evaluate stress and add it to matrix contribution
	s += FEEllipsoidalFiberDistributionStress(m_ksi, m_beta, nint, m_cth, m_sth, 
											  m_cph, m_sph, m_w, J, F, Q);
		
	return s;
}

tens4ds FERandomFiberDonnanEquilibrium::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	const int nint = (m_nres == 0? NSTL  : NSTH  );
	
	// --- M A T R I X   C O N T R I B U T I O N ---

	tens4ds c = FEDonnanEquilibriumTangent(m_phiwr, m_cFr, m_Rgas, m_Tabs, m_bosm, J);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// get the element's local coordinate system
	mat3d& Q = pt.Q;
	
	// evaluate stress and add it to matrix contribution
	c += FEEllipsoidalFiberDistributionTangent(m_ksi, m_beta, nint, m_cth, m_sth, 
											   m_cph, m_sph, m_w, J, F, Q);

	return c;
}

double FERandomFiberDonnanEquilibrium::BulkModulus()
{
	// Evaluate bulk modulus in reference configuration
	// Only the matrix contributes to the bulk modulus
	return FEDonnanEquilibriumBulkModulus(m_phiwr, m_cFr, m_Rgas, m_Tabs, m_bosm);
}
