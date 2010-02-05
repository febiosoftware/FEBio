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
REGISTER_MATERIAL(FERandomFiberDonnanEquilibrium, "random fiber Donnan equilibrium");

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
	
	if (m_phiwr < 0 || m_phiwr > 1) throw MaterialError("phiw0 must be between 0. and 1.");
	if (m_Rgas < 0) throw MaterialError("R must be positive.");
	if (m_Tabs < 0) throw MaterialError("T must be positive.");
	if (m_ksi[0] < 0) throw MaterialError("ksi1 must be positive.");
	if (m_ksi[1] < 0) throw MaterialError("ksi2 must be positive.");
	if (m_ksi[2] < 0) throw MaterialError("ksi3 must be positive.");
	if (m_beta[0] < 2) throw MaterialError("beta1 must be greater than 2.");
	if (m_beta[1] < 2) throw MaterialError("beta1 must be greater than 2.");
	if (m_beta[2] < 2) throw MaterialError("beta1 must be greater than 2.");
}

mat3ds FERandomFiberDonnanEquilibrium::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	const int nint = (m_nres == 0? NSTL  : NSTH  );
	
	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
			
	// --- M A T R I X   C O N T R I B U T I O N ---

	// calculate fixed charge density in current configuration
	double cF = m_phiwr*m_cFr/(J-1+m_phiwr);
	
	// calculate osmotic pressure
	double p = m_Rgas*m_Tabs*(sqrt(cF*cF+m_bosm*m_bosm) - m_bosm);

	// calculate T = -p*I
	mat3dd I(1.0);	// identity tensor
	mat3ds s = -p*I;
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// calculate right Cauchy-Green tensor
	mat3ds C = pt.RightCauchyGreen();
	
	// get the element's local coordinate system
	mat3d& Q = pt.Q;
	
	// loop over all integration points
	double ksi, beta;
	vec3d n0e, n0a, nt;
	double In, Wl;
	const double eps = 0;
	for (int n=0; n<nint; ++n)
	{
		// set the global fiber direction in reference configuration
		n0e.x = m_cth[n]*m_sph[n];
		n0e.y = m_sth[n]*m_sph[n];
		n0e.z = m_cph[n];
		
		// Calculate In = n0e*C*n0e
		In = n0e*(C*n0e);
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// get the global spatial fiber direction in current configuration
			nt = (F*n0e)/sqrt(In);
			
			// calculate the outer product of nt
			mat3ds N = dyad(nt);
			
			// get the local material fiber direction in reference configuration
			n0a = Q*n0e;
			
			// calculate material coefficients
			ksi  = 1.0 / sqrt(SQR(n0a.x / m_ksi [0]) + SQR(n0a.y / m_ksi [1]) + SQR(n0a.z / m_ksi [2]));
			beta = 1.0 / sqrt(SQR(n0a.x / m_beta[0]) + SQR(n0a.y / m_beta[1]) + SQR(n0a.z / m_beta[2]));
			
			// calculate strain energy derivative
			Wl = beta*ksi*pow(In - 1.0, beta-1.0);
			
			// calculate the stress
			s += N*(2.0/J*In*Wl*m_w[n]);
		}
		
	}
	
	// --- END FIBER CONTRIBUTION --
		
	return s;
}

tens4ds FERandomFiberDonnanEquilibrium::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	const int nint = (m_nres == 0? NSTL  : NSTH  );
	
	// right Cauchy-Green tensor: C = Ft*F
	mat3ds C = pt.RightCauchyGreen();
		
	// --- M A T R I X   C O N T R I B U T I O N ---
	
	// calculate fixed charge density in current configuration
	double cF = m_phiwr*m_cFr/(J-1+m_phiwr);
	
	// calculate osmotic pressure
	double tosm = sqrt(cF*cF+m_bosm*m_bosm);	// tissue osmolarity
	double p = m_Rgas*m_Tabs*(tosm - m_bosm);	// osmotic pressure
	
	// calculate derivative of osmotic pressure w.r.t. J
	double bpi = m_Rgas*m_Tabs*J*cF*cF/(J-1+m_phiwr)/tosm;

	mat3dd I(1.0);	// Identity
	
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);

	// calculate tangent osmotic modulus
	tens4ds c = bpi*IxI + p*(2.0*I4 - IxI);
	
	// --- F I B E R   C O N T R I B U T I O N ---
	
	// get the element's local coordinate system
	mat3d& Q = pt.Q;
	
	// loop over all integration points
	double ksi, beta;
	vec3d n0e, n0a, nt;
	double In, Wll;
	const double eps = 0;
	tens4ds cf, cfw; cf.zero();
	mat3ds N2;
	tens4ds N4;
	for (int n=0; n<nint; ++n)
	{
		// set the global fiber direction in reference configuration
		n0e.x = m_cth[n]*m_sph[n];
		n0e.y = m_sth[n]*m_sph[n];
		n0e.z = m_cph[n];
		
		// Calculate In = n0e*C*n0e
		In = n0e*(C*n0e);
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// get the global spatial fiber direction in current configuration
			nt = (F*n0e)/sqrt(In);
			
			// calculate the outer product of nt
			N2 = dyad(nt);
			
			// get the local material fiber direction in reference configuration
			n0a = Q*n0e;
			
			// calculate material coefficients in local fiber direction
			ksi  = 1.0 / sqrt(SQR(n0a.x / m_ksi [0]) + SQR(n0a.y / m_ksi [1]) + SQR(n0a.z / m_ksi [2]));
			beta = 1.0 / sqrt(SQR(n0a.x / m_beta[0]) + SQR(n0a.y / m_beta[1]) + SQR(n0a.z / m_beta[2]));
			
			// calculate strain energy derivative
			Wll = beta*(beta-1.0)*ksi*pow(In - 1.0, beta-2.0);
			
			N2 = dyad(nt);
			N4 = dyad1s(N2);
						
			c += N4*(4.0/J*In*In*Wll*m_w[n]);
		}
	}
	
	return c;
}
