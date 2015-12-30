#include "stdafx.h"
#include "FEEllipsoidalFiberDistribution.h"
#include "FEContinuousFiberDistribution.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
// FEEllipsoidalFiberDistribution
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_PARAMETER_LIST(FEEllipsoidalFiberDistribution, FEElasticMaterial)
	ADD_PARAMETERV2(m_beta, FE_PARAM_DOUBLE, 3, FE_RANGE_GREATER_OR_EQUAL(2.0), "beta");
	ADD_PARAMETERV2(m_ksi , FE_PARAM_DOUBLE, 3, FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi" );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
mat3ds FEEllipsoidalFiberDistribution::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;

	// get the element's local coordinate system
	mat3d Q = pt.m_Q;

	// loop over all integration points
	vec3d n0e, n0a, n0q, nt;
	double In, Wl;
	const double eps = 0;
	mat3ds s;
	s.zero();

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
		
		// --- quadrant 1,1,1 ---

		// rotate to reference configuration
		n0e = Q*n0a;

		// get the global spatial fiber direction in current configuration
		nt = F*n0e;

		// Calculate In = n0e*C*n0e
		In = nt*nt;
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate strain energy derivative
			Wl = beta*ksi*pow(In - 1.0, beta-1.0);

			// calculate the stress
			s += dyad(nt)*(Wl*wn);
		}

		// --- quadrant -1,1,1 ---
		n0q = vec3d(-n0a.x, n0a.y, n0a.z);

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

			// calculate the stress
			s += dyad(nt)*(Wl*wn);
		}

		// --- quadrant -1,-1,1 ---
		n0q = vec3d(-n0a.x, -n0a.y, n0a.z);

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

			// calculate the stress
			s += dyad(nt)*(Wl*wn);
		}

		// --- quadrant 1,-1,1 ---
		n0q = vec3d(n0a.x, -n0a.y, n0a.z);

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

			// calculate the stress
			s += dyad(nt)*(Wl*wn);
		}
	}

	// we multiply by two to add contribution from other half-sphere
	return s*(4.0/J);
}

//-----------------------------------------------------------------------------
tens4ds FEEllipsoidalFiberDistribution::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;

	// get the element's local coordinate system
	mat3d Q = pt.m_Q;

	// loop over all integration points
	vec3d n0e, n0a, n0q, nt;
	double In, Wll;
	const double eps = 0;
	tens4ds cf, cfw; cf.zero();
	mat3ds N2;
	tens4ds N4;
	tens4ds c;
	c.zero();
	
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

		// --- quadrant 1,1,1 ---

		// rotate to reference configuration
		n0e = Q*n0a;

		// get the global spatial fiber direction in current configuration
		nt = F*n0e;

		// Calculate In = n0e*C*n0e
		In = nt*nt;
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate strain energy derivative
			Wll = beta*(beta-1.0)*ksi*pow(In - 1.0, beta-2.0);

			N2 = dyad(nt);
			N4 = dyad1s(N2);
			
			c += N4*(Wll*wn);
		}

		// --- quadrant -1,1,1 ---

		n0q = vec3d(-n0a.x, n0a.y, n0a.z);

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
			Wll = beta*(beta-1.0)*ksi*pow(In - 1.0, beta-2.0);
			
			N2 = dyad(nt);
			N4 = dyad1s(N2);
			
			c += N4*(Wll*wn);
		}

		// --- quadrant -1,-1,1 ---

		n0q = vec3d(-n0a.x, -n0a.y, n0a.z);

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
			Wll = beta*(beta-1.0)*ksi*pow(In - 1.0, beta-2.0);
			
			N2 = dyad(nt);
			N4 = dyad1s(N2);
			
			c += N4*(Wll*wn);
		}

		// --- quadrant 1,-1,1 ---

		n0q = vec3d(n0a.x, -n0a.y, n0a.z);

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
			Wll = beta*(beta-1.0)*ksi*pow(In - 1.0, beta-2.0);
			
			N2 = dyad(nt);
			N4 = dyad1s(N2);
			
			c += N4*(Wll*wn);
		}
	}

	// multiply by two to integrate over other half of sphere
	return c*(2.0*4.0/J);
}

//-----------------------------------------------------------------------------
double FEEllipsoidalFiberDistribution::StrainEnergyDensity(FEMaterialPoint& mp)
{
    double sed = 0.0;
    
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
    
	// get the element's local coordinate system
	mat3d Q = pt.m_Q;
    
	// loop over all integration points
	vec3d n0e, n0a, n0q, nt;
	double In, W;
	const double eps = 0;
    
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
		
		// --- quadrant 1,1,1 ---
        
		// rotate to reference configuration
		n0e = Q*n0a;
        
		// get the global spatial fiber direction in current configuration
		nt = F*n0e;
        
		// Calculate In = n0e*C*n0e
		In = nt*nt;
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate strain energy
			W = ksi*pow(In - 1.0, beta);
			sed += W*wn;
		}
        
		// --- quadrant -1,1,1 ---
		n0q = vec3d(-n0a.x, n0a.y, n0a.z);
        
		// rotate to reference configuration
		n0e = Q*n0q;
        
		// get the global spatial fiber direction in current configuration
		nt = F*n0e;
        
		// Calculate In = n0e*C*n0e
		In = nt*nt;
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate strain energy
			W = ksi*pow(In - 1.0, beta);
			sed += W*wn;
		}
        
		// --- quadrant -1,-1,1 ---
		n0q = vec3d(-n0a.x, -n0a.y, n0a.z);
        
		// rotate to reference configuration
		n0e = Q*n0q;
        
		// get the global spatial fiber direction in current configuration
		nt = F*n0e;
        
		// Calculate In = n0e*C*n0e
		In = nt*nt;
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate strain energy
			W = ksi*pow(In - 1.0, beta);
			sed += W*wn;
		}
        
		// --- quadrant 1,-1,1 ---
		n0q = vec3d(n0a.x, -n0a.y, n0a.z);
        
		// rotate to reference configuration
		n0e = Q*n0q;
        
		// get the global spatial fiber direction in current configuration
		nt = F*n0e;
        
		// Calculate In = n0e*C*n0e
		In = nt*nt;
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate strain energy
			W = ksi*pow(In - 1.0, beta);
			sed += W*wn;
		}
	}
    
	// multiply by two to integrate over other half of sphere
    return sed*2.0;
}

//-----------------------------------------------------------------------------
// FEEllipsoidalFiberDistributionOld
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_PARAMETER_LIST(FEEllipsoidalFiberDistributionOld, FEElasticMaterial)
	ADD_PARAMETERV2(m_beta, FE_PARAM_DOUBLE, 3, FE_RANGE_GREATER_OR_EQUAL(2.0), "beta");
	ADD_PARAMETERV2(m_ksi , FE_PARAM_DOUBLE, 3, FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi" );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// FEEllipsoidalFiberDistributionOld
//-----------------------------------------------------------------------------

bool FEEllipsoidalFiberDistributionOld::Init()
{
	if (FEElasticMaterial::Init() == false) return false;

	m_bfirst = true;
	
	if (m_bfirst)
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
		
		m_bfirst = false;
	}

	return true;
}

//-----------------------------------------------------------------------------
mat3ds FEEllipsoidalFiberDistributionOld::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the element's local coordinate system
	mat3d QT = (pt.m_Q).transpose();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;

	// loop over all integration points
	double ksi, beta;
	vec3d n0e, n0a, nt;
	double In, Wl;
	const double eps = 0;
	mat3ds C = pt.RightCauchyGreen();
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
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEEllipsoidalFiberDistributionOld::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the element's local coordinate system
	mat3d QT = (pt.m_Q).transpose();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;

	// loop over all integration points
	double ksi, beta;
	vec3d n0e, n0a, nt;
	double In, Wll;
	const double eps = 0;
	tens4ds cf, cfw; cf.zero();
	mat3ds N2;
	tens4ds N4;
	tens4ds c;
	c.zero();
	
	// right Cauchy-Green tensor: C = Ft*F
	mat3ds C = (F.transpose()*F).sym();

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
			Wll = beta*(beta-1.0)*ksi*pow(In - 1.0, beta-2.0);
			
			N2 = dyad(nt);
			N4 = dyad1s(N2);
			
			c += N4*(4.0/J*In*In*Wll*m_w[n]);
		}
	}
	
	return c;
}

//-----------------------------------------------------------------------------
double FEEllipsoidalFiberDistributionOld::StrainEnergyDensity(FEMaterialPoint& mp)
{
    double sed = 0.0;
    
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the element's local coordinate system
	mat3d QT = (pt.m_Q).transpose();
	
	// deformation gradient
	mat3d &F = pt.m_F;
    
	// loop over all integration points
	double ksi, beta;
	vec3d n0e, n0a, nt;
	double In, W;
	const double eps = 0;
	mat3ds C = pt.RightCauchyGreen();
    
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
			
			// get the local material fiber direction in reference configuration
			n0a = QT*n0e;
			
			// calculate material coefficients
			ksi  = 1.0 / sqrt(SQR(n0a.x / m_ksi [0]) + SQR(n0a.y / m_ksi [1]) + SQR(n0a.z / m_ksi [2]));
			beta = 1.0 / sqrt(SQR(n0a.x / m_beta[0]) + SQR(n0a.y / m_beta[1]) + SQR(n0a.z / m_beta[2]));
			
			// calculate strain energy density
			W = ksi*pow(In - 1.0, beta);
			sed += W*m_w[n];
		}
		
	}
    
    return sed;
}
