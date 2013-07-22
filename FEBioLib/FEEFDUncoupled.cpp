#include "stdafx.h"
#include "FEEFDUncoupled.h"

// The following file contains the integration points and weights
// for the integration over a unit sphere in spherical coordinates
#include "geodesic.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
// FEEFDUncoupled
//-----------------------------------------------------------------------------

// register the material with the framework
REGISTER_MATERIAL(FEEFDUncoupled, "EFD uncoupled");

// define the material parameters
BEGIN_PARAMETER_LIST(FEEFDUncoupled, FEUncoupledMaterial)
	ADD_PARAMETERV(m_beta, FE_PARAM_DOUBLEV, 3, "beta");
	ADD_PARAMETERV(m_ksi , FE_PARAM_DOUBLEV, 3, "ksi" );
END_PARAMETER_LIST();

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
}

//-----------------------------------------------------------------------------
mat3ds FEEFDUncoupled::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	double J = pt.m_J;
	// deviatoric deformation gradient
	mat3d F = pt.m_F*pow(J,-1.0/3.0);
	
	// get the element's local coordinate system
	mat3d Q = pt.m_Q;
	
	// loop over all integration points
	vec3d n0e, n0a, n0q, nt;
	double In, Wl;
	const double eps = 0;
	mat3ds C = pt.DevRightCauchyGreen();
	mat3ds s;
	s.zero();

	const int li[4][3] = {
		{ 1, 1, 1},
		{-1, 1, 1},
		{-1,-1, 1},
		{ 1,-1, 1}
	};

	const int nint = 45;
	for (int n=0; n<nint; ++n)
	{
		// set the global fiber direction in material coordinate system
		n0a.x = XYZ2[n][0];
		n0a.y = XYZ2[n][1];
		n0a.z = XYZ2[n][2];
		double wn = XYZ2[n][3];

		// calculate material coefficients
		// TODO: There is an obvious optimization opportunity here, since the values of ksi
		//       and beta can be precalculated and reused. I have not done this yet since I
		//       need to figure out how to initialize the material parameters for each time
		//       step (instead of once at the start) in case the values depend on load curves.
		double ksi  = 1.0 / sqrt(SQR(n0a.x / m_ksi [0]) + SQR(n0a.y / m_ksi [1]) + SQR(n0a.z / m_ksi [2]));
		double beta = 1.0 / sqrt(SQR(n0a.x / m_beta[0]) + SQR(n0a.y / m_beta[1]) + SQR(n0a.z / m_beta[2]));

		// loop over the four quadrants
		for (int l=0; l<4; ++l)
		{
			n0q = vec3d(li[l][0]*n0a.x, li[l][1]*n0a.y, li[l][2]*n0a.z);

			// rotate to reference configuration
			n0e = Q*n0q;

			// get the global spatial fiber direction in current configuration
			nt = F*n0e;

			// Calculate In = n0e*C*n0e
			In = nt*nt;
		
			// only take fibers in tension into consideration
			if (In > 1. + eps)
			{
				// calculate the outer product of nt
				mat3ds N = dyad(nt);
			
				// calculate strain energy derivative
				Wl = beta*ksi*pow(In - 1.0, beta-1.0);
			
				// calculate the stress
				s += N*(Wl*wn);
			}
		}
	}
	// don't forget to multiply by two to include the other half-sphere
	return s.dev()*(4.0/J);
}


//-----------------------------------------------------------------------------
tens4ds FEEFDUncoupled::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	double J = pt.m_J;
	// deviatoric deformation gradient
	mat3d F = pt.m_F*pow(J,-1.0/3.0);
	
	// get the element's local coordinate system
	mat3d Q = pt.m_Q;
	
	// loop over all integration points
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

	const int li[4][3] = {
		{ 1, 1, 1},
		{-1, 1, 1},
		{-1,-1, 1},
		{ 1,-1, 1}
	};
	
	const int nint = 45;
	for (int n=0; n<nint; ++n)
	{
		// set the global fiber direction in material coordinate system
		n0a.x = XYZ2[n][0];
		n0a.y = XYZ2[n][1];
		n0a.z = XYZ2[n][2];
		double wn = XYZ2[n][3];

		// calculate material coefficients
		double ksi  = 1.0 / sqrt(SQR(n0a.x / m_ksi [0]) + SQR(n0a.y / m_ksi [1]) + SQR(n0a.z / m_ksi [2]));
		double beta = 1.0 / sqrt(SQR(n0a.x / m_beta[0]) + SQR(n0a.y / m_beta[1]) + SQR(n0a.z / m_beta[2]));

		for (int l=0; l<4; ++l)
		{
			vec3d n0q = vec3d(li[l][0]*n0a.x, li[l][1]*n0a.y, li[l][2]*n0a.z);

			// rotate to reference configuration
			n0e = Q*n0q;

			// get the global spatial fiber direction in current configuration
			nt = F*n0e;

			// Calculate In = n0e*C*n0e
			In = nt*nt;
	
			// only take fibers in tension into consideration
			if (In > 1. + eps)
			{
			
				// calculate strain energy derivative
				Wl = beta*ksi*pow(In - 1.0, beta-1.0);
				Wll = beta*(beta-1.0)*ksi*pow(In - 1.0, beta-2.0);
			
				// calculate the outer product of nt
				N2 = dyad(nt);
				N4 = dyad1s(N2);

				// calculate the stress
				s += N2*(2.0/J*Wl*wn);

				// calculate tangent
				c += N4*(4.0/J*Wll*wn);
			}
		}
	}

	// don't forget to multiply by two to include the other half-sphere
	s *= 2.0;
	c *= 2.0;
	
	// This is the final value of the elasticity tensor
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	c += ((I4+IxI/3.0)*s.tr() - dyad1s(I,s))*(2./3.)
	- (ddots(IxI, c)-IxI*(c.tr()/3.))/3.;
	
	return c;
}
