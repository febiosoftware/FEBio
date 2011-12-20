#include "stdafx.h"
#include "FEEFDUncoupled.h"

// The following file contains the integration points and weights
// for the integration over a unit sphere in spherical coordinates
#include "geodesic.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

// we store the cos and sin of the angles here
int FEEFDUncoupled::m_nres = 0;
double FEEFDUncoupled::m_cth[NSTH];
double FEEFDUncoupled::m_sth[NSTH];
double FEEFDUncoupled::m_cph[NSTH];
double FEEFDUncoupled::m_sph[NSTH];
double FEEFDUncoupled::m_w[NSTH];

// define the material parameters
BEGIN_PARAMETER_LIST(FEEFDUncoupled, FEUncoupledMaterial)
	ADD_PARAMETERV(m_beta, FE_PARAM_DOUBLEV, 3, "beta");
	ADD_PARAMETERV(m_ksi , FE_PARAM_DOUBLEV, 3, "ksi" );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// FEEFDUncoupled
//-----------------------------------------------------------------------------

void FEEFDUncoupled::Init()
{
	if (m_unstable) throw MaterialError("This fibrous material is unstable (collapses on itself) when used alone.  Combine it in a solid mixture with a material that can serve as a ground matrix.");
	if (m_ksi[0] < 0) throw MaterialError("ksi1 must be positive.");
	if (m_ksi[1] < 0) throw MaterialError("ksi2 must be positive.");
	if (m_ksi[2] < 0) throw MaterialError("ksi3 must be positive.");
	if (m_beta[0] < 2) throw MaterialError("beta1 must be greater than 2.");
	if (m_beta[1] < 2) throw MaterialError("beta1 must be greater than 2.");
	if (m_beta[2] < 2) throw MaterialError("beta1 must be greater than 2.");
	
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

//-----------------------------------------------------------------------------
mat3ds FEEFDUncoupled::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	double J = pt.J;
	// deviatoric deformation gradient
	mat3d F = pt.F*pow(J,-1.0/3.0);
	
	// get the element's local coordinate system
	mat3d QT = (pt.Q).transpose();
	
	// loop over all integration points
	double ksi, beta;
	vec3d n0e, n0a, nt;
	double In, Wl;
	const double eps = 0;
	mat3ds C = pt.DevRightCauchyGreen();
	mat3ds s;
	s.zero();
	
	const int nint = (m_nres == 0? NSTL  : NSTH  );
	
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
			n0a = QT*n0e;
			
			// calculate material coefficients
			ksi  = 1.0 / sqrt(SQR(n0a.x / m_ksi [0]) + SQR(n0a.y / m_ksi [1]) + SQR(n0a.z / m_ksi [2]));
			beta = 1.0 / sqrt(SQR(n0a.x / m_beta[0]) + SQR(n0a.y / m_beta[1]) + SQR(n0a.z / m_beta[2]));
			
			// calculate strain energy derivative
			Wl = beta*ksi*pow(In - 1.0, beta-1.0);
			
			// calculate the stress
			s += N*(2.0/J*In*Wl*m_w[n]);
		}
		
	}
	return s.dev();
}

//-----------------------------------------------------------------------------
tens4ds FEEFDUncoupled::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	double J = pt.J;
	// deviatoric deformation gradient
	mat3d F = pt.F*pow(J,-1.0/3.0);
	
	// get the element's local coordinate system
	mat3d QT = (pt.Q).transpose();
	
	// loop over all integration points
	double ksi, beta;
	vec3d n0e, n0a, nt;
	double In, Wl, Wll;
	const double eps = 0;
	mat3ds C = pt.DevRightCauchyGreen();
	mat3ds s;
	tens4ds cf, cfw; cf.zero();
	mat3ds N2;
	mat3dd I(1);
	tens4ds N4;
	tens4ds c;

	s.zero();
	c.zero();
	
	const int nint = (m_nres == 0? NSTL  : NSTH  );
	
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
			n0a = QT*n0e;
			
			// calculate material coefficients in local fiber direction
			ksi  = 1.0 / sqrt(SQR(n0a.x / m_ksi [0]) + SQR(n0a.y / m_ksi [1]) + SQR(n0a.z / m_ksi [2]));
			beta = 1.0 / sqrt(SQR(n0a.x / m_beta[0]) + SQR(n0a.y / m_beta[1]) + SQR(n0a.z / m_beta[2]));
			
			// calculate strain energy derivative
			Wl = beta*ksi*pow(In - 1.0, beta-1.0);
			Wll = beta*(beta-1.0)*ksi*pow(In - 1.0, beta-2.0);
			
			// calculate the stress
			s += N2*(2.0/J*In*Wl*m_w[n]);

			// calculate the elasticity tensor
			N4 = dyad1s(N2);
			
			c += N4*(4.0/J*In*In*Wll*m_w[n]);
		}
	}
	
	// This is the final value of the elasticity tensor
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	c += - 1./3.*(ddots(c,IxI) - IxI*(c.tr()/3.))
	+ 2./3.*((I4-IxI/3.)*s.tr()-dyad1s(s.dev(),I));
	
	return c;
}
