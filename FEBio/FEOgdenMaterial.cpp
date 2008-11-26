#include "stdafx.h"
#include "FEOgdenMaterial.h"

REGISTER_MATERIAL(FEOgdenMaterial, "Ogden");

BEGIN_PARAMETER_LIST(FEOgdenMaterial, FEIncompressibleMaterial);
	ADD_PARAMETER(m_c[0], FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(m_c[1], FE_PARAM_DOUBLE, "c2");
	ADD_PARAMETER(m_c[2], FE_PARAM_DOUBLE, "c3");
	ADD_PARAMETER(m_c[3], FE_PARAM_DOUBLE, "c4");
	ADD_PARAMETER(m_c[4], FE_PARAM_DOUBLE, "c5");
	ADD_PARAMETER(m_c[5], FE_PARAM_DOUBLE, "c6");
	ADD_PARAMETER(m_m[0], FE_PARAM_DOUBLE, "m1");
	ADD_PARAMETER(m_m[1], FE_PARAM_DOUBLE, "m2");
	ADD_PARAMETER(m_m[2], FE_PARAM_DOUBLE, "m3");
	ADD_PARAMETER(m_m[3], FE_PARAM_DOUBLE, "m4");
	ADD_PARAMETER(m_m[4], FE_PARAM_DOUBLE, "m5");
	ADD_PARAMETER(m_m[5], FE_PARAM_DOUBLE, "m6");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor
FEOgdenMaterial::FEOgdenMaterial() : FEIncompressibleMaterial(FE_OGDEN_MATERIAL)
{
	for (int i=0; i<MAX_TERMS; ++i)
	{
		m_c[i] = 0;
		m_m[i] = 1;
	}

	m_eps = 1e-12;
}

//-----------------------------------------------------------------------------
//! data initialization and checking
void FEOgdenMaterial::Init()
{
	for (int i=0; i<MAX_TERMS; ++i) 
		if (m_m[i] == 0) throw MaterialError("Invalid value for m%d", i+1);
}

//-----------------------------------------------------------------------------

void FEOgdenMaterial::EigenValues(mat3ds& A, double l[3], vec3d r[3], const double eps)
{
	A.eigen(l, r);

	// correct for numerical inaccuracy
	double d01 = fabs(l[0] - l[1]);
	double d12 = fabs(l[1] - l[2]);
	double d02 = fabs(l[0] - l[2]);

	if (d01 < eps) l[1] = l[0]; //= 0.5*(l[0]+l[1]);
	if (d02 < eps) l[2] = l[0]; //= 0.5*(l[0]+l[2]);
	if (d12 < eps) l[2] = l[1]; //= 0.5*(l[1]+l[2]);

/*
	// the following logic was "borrowed" from nike3d.
	// in this case, we slighly modify the eigenvalues
	// so that we only deal with the distinct case.
	if (d01 < 1e-15) l[1] = l[0] + 1e-14;
	if (d02 < 1e-15) l[2] = l[0] + 2e-14;
	if (d12 < 1e-15)
	{
		if (l[1] > l[0])
			l[2] = l[1] + 3e-14;
		else
			l[2] = l[1] - 3e-14;
	}
*/
}

//-----------------------------------------------------------------------------
//! Calculates the Cauchy stress
mat3ds FEOgdenMaterial::Stress(FEMaterialPoint &mp)
{
	double beta, beta1, beta3, D;
	int i, j, k, l;

	mat3dd I(1.0);	// identity tensor

	// extract elastic material data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// average pressure
	double p = pt.avgp;

	// jacobian
	double J = pt.J;
	double Jm13 = pow(J, -1.0/3.0);

	// get the left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();

	// get the eignvalues and eigenvectors of b
	double lam2[3];	// these are the squares of the eigenvalues of F
	vec3d ev[3];
	EigenValues(b, lam2, ev, m_eps);

	// get the deviatoric eigenvalues of b
	double ld1 = Jm13*sqrt(lam2[0]);
	double ld2 = Jm13*sqrt(lam2[1]);
	double ld3 = Jm13*sqrt(lam2[2]);

	// calculate the powers of (deviatoric) eigenvalues
	double lamd[3][MAX_TERMS];
	for (i=0; i<MAX_TERMS; ++i)
	{
		lamd[0][i] = pow(ld1, m_m[i]);
		lamd[1][i] = pow(ld2, m_m[i]);
		lamd[2][i] = pow(ld3, m_m[i]);
	}

	// stress deviator
	mat3ds sd;
	sd.zero();

	mat3ds m;

	// consider the three different cases
	if ((fabs(lam2[0]-lam2[1])>=m_eps) && (fabs(lam2[0] - lam2[2])>=m_eps) && (fabs(lam2[1] - lam2[2])>=m_eps))
	{
		// all three are different
		for (i=0; i<3; ++i)
		{
			j = (i+1)%3;
			k = (j+1)%3;

			beta = 0;
			for (l=0; l<MAX_TERMS; ++l) beta += m_c[l]/m_m[l]*(lamd[i][l] - (lamd[0][l]+lamd[1][l]+lamd[2][l])/3.0);

			D = (lam2[i] - lam2[j])*(lam2[i] - lam2[k]);
			m = (((b - I*lam2[j])*(b - I*lam2[k]))) / D;

			sd += m*(beta/J);
		}
	}
	else if ((fabs(lam2[0] - lam2[1]) >= m_eps) || 
			 (fabs(lam2[0] - lam2[2]) >= m_eps) ||
			 (fabs(lam2[1] - lam2[2]) >= m_eps ))
	{
		// two are different
		if      (fabs(lam2[0] - lam2[1]) < m_eps) i=2;
		else if (fabs(lam2[0] - lam2[2]) < m_eps) i=1;
		else i=0;

		j = (i+1)%3;
		k = (j+1)%3;

		beta1 = beta3 = 0;
		for (l=0; l<MAX_TERMS; ++l)
		{
			beta1 += m_c[l]/m_m[l]*(lamd[j][l] - (lamd[0][l]+lamd[1][l]+lamd[2][l])/3.0);
			beta3 += m_c[l]/m_m[l]*(lamd[i][l] - (lamd[0][l]+lamd[1][l]+lamd[2][l])/3.0);
		}

		D = (lam2[i] - lam2[j])*(lam2[i] - lam2[k]);
		m = (((b - I*lam2[j])*(b - I*lam2[k]))) / D;

		sd += I*(beta1/J) + m*((beta3 - beta1)/J);
	}

	return (I*p+sd.dev());
}

//-----------------------------------------------------------------------------
//! Calculates the spatial tangent
void FEOgdenMaterial::Tangent(double D[6][6], FEMaterialPoint& mp)
{
	int i, j, k, l, n, i0, i1, j0, j1;
	double Di, Dpi, beta, gab;
	double beta1, beta3, g11, g33, g13;

	// get the elastic material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();

	// get the jacobian
	double J = pt.J;
	double Jm13 = pow(J, -1.0/3.0);

	// get the eigenvalues and vectors of b
	double lam2[3];
	vec3d ev[3];
	EigenValues(b, lam2, ev, m_eps);

	// get the eigenvalues of F
	double lam[3];
	lam[0] = sqrt(lam2[0]);
	lam[1] = sqrt(lam2[1]);
	lam[2] = sqrt(lam2[2]);

	// get the deviatoric eigenvalues of F
	double ld[3];
	ld[0] = Jm13*lam[0];
	ld[1] = Jm13*lam[1];
	ld[2] = Jm13*lam[2];

	// invariants
	double I1 = lam2[0]+lam2[1]+lam2[2];
	double I3 = lam2[0]*lam2[1]*lam2[2];

	// calculate the powers of eigenvalues
	double lamd[3][MAX_TERMS];
	for (i=0; i<MAX_TERMS; ++i)
	{
		lamd[0][i] = pow(ld[0], m_m[i]);
		lamd[1][i] = pow(ld[1], m_m[i]);
		lamd[2][i] = pow(ld[2], m_m[i]);
	}

	memset(D, 0, sizeof(double)* 6*6);

	mat3dd I(1.0);	// unit matrix

	// 4-th order unit tensor
	double I4[6][6] = {
		{ 1, 0, 0, 0, 0, 0 },
		{ 0, 1, 0, 0, 0, 0 },
		{ 0, 0, 1, 0, 0, 0 },
		{ 0, 0, 0, .5,  0,  0 },
		{ 0, 0, 0,  0, .5,  0 },
		{ 0, 0, 0,  0,  0, .5 }
	};

	// b tensor
	double B[3][3];
	B[0][0] = b(0,0); B[0][1] = b(0,1); B[0][2] = b(0,2);
	B[1][0] = b(1,0); B[1][1] = b(1,1); B[1][2] = b(1,2);
	B[2][0] = b(2,0); B[2][1] = b(2,1); B[2][2] = b(2,2);

	// Ib tensor
	// Note that we only assign the upper triangular part
	double Ib[6][6] = {0};
	Ib[0][0] = B[0][0]*B[0][0];	// = Ib[0][0][0][0]
	Ib[1][1] = B[1][1]*B[1][1];	// = Ib[1][1][1][1]
	Ib[2][2] = B[2][2]*B[2][2];	// = Ib[2][2][2][2]

	Ib[0][1] = B[0][1]*B[0][1]; // = Ib[0][0][1][1]
	Ib[0][2] = B[0][2]*B[0][2]; // = Ib[0][0][2][2]
	Ib[1][2] = B[1][2]*B[1][2]; // = Ib[1][1][2][2]

	Ib[3][3] = 0.5*(B[0][0]*B[1][1] + B[0][1]*B[1][0]); // = Ib[0][1][0][1]
	Ib[4][4] = 0.5*(B[1][1]*B[2][2] + B[1][2]*B[2][1]); // = Ib[1][2][1][2]
	Ib[5][5] = 0.5*(B[0][0]*B[2][2] + B[0][2]*B[2][0]); // = Ib[0][2][0][2]

	Ib[0][3] = 0.5*(B[0][0]*B[0][1] + B[0][1]*B[0][0]); // = Ib[0][0][0][1]
	Ib[0][4] = 0.5*(B[0][1]*B[0][2] + B[0][2]*B[0][1]); // = Ib[0][0][1][2]
	Ib[0][5] = 0.5*(B[0][0]*B[0][2] + B[0][2]*B[0][0]); // = Ib[0][0][0][2]

	Ib[1][3] = 0.5*(B[1][0]*B[1][1] + B[1][1]*B[1][0]); // = Ib[1][1][0][1]
	Ib[1][4] = 0.5*(B[1][1]*B[1][2] + B[1][2]*B[1][1]); // = Ib[1][1][1][2]
	Ib[1][5] = 0.5*(B[1][0]*B[1][2] + B[1][2]*B[1][0]); // = Ib[1][1][0][2]

	Ib[2][3] = 0.5*(B[2][0]*B[2][1] + B[2][1]*B[2][0]); // = Ib[2][2][0][1]
	Ib[2][4] = 0.5*(B[2][1]*B[2][2] + B[2][2]*B[2][1]); // = Ib[2][2][1][2]
	Ib[2][5] = 0.5*(B[2][0]*B[2][2] + B[2][2]*B[2][0]); // = Ib[2][2][0][2]

	Ib[3][4] = 0.5*(B[0][1]*B[1][2] + B[0][2]*B[1][1]); // = Ib[0][1][1][2]
	Ib[3][5] = 0.5*(B[0][0]*B[1][2] + B[0][2]*B[1][0]); // = Ib[0][1][0][2]
	Ib[4][5] = 0.5*(B[1][0]*B[2][2] + B[1][2]*B[2][0]); // = Ib[1][2][0][2]

	// 2nd order identity
	double I2[3][3] = {{1,0,0},{0,1,0},{0,0,1}};

	// IxI tensor
	double IxI[6][6] = {
		{ 1, 1, 1, 0, 0, 0},
		{ 1, 1, 1, 0, 0, 0},
		{ 1, 1, 1, 0, 0, 0},
		{ 0, 0, 0, 0, 0, 0},
		{ 0, 0, 0, 0, 0, 0},
		{ 0, 0, 0, 0, 0, 0}
	};

	mat3ds ma, mb;
	double dm[6][6] = {0};

	double BxB[6][6], BxM[6][6], MxB[6][6], MxM[6][6], IxM[6][6], MxI[6][6];

	const int LUT[6][2] = {{0,0}, {1,1}, {2,2}, {0,1}, {1,2}, {0,2}};

	// consider the three distinct cases
	if ((fabs(lam2[0] - lam2[1]) < m_eps) && (fabs(lam2[0] - lam2[2]) < m_eps))
	{
		// all three are the same
		double g = 0;
		for (i=0; i<MAX_TERMS; ++i) g += m_c[i]*lamd[0][i];

		D[0][0] += g/J*(2.0/3.0);	// C[0][0][0][0]
		D[1][1] += g/J*(2.0/3.0);	// C[1][1][1][1]
		D[2][2] += g/J*(2.0/3.0);	// C[2][2][2][2]

		D[0][1] += g/J*(-1.0/3.0);	// C[0][0][1][1]
		D[0][2] += g/J*(-1.0/3.0);	// C[0][0][2][2]
		D[1][2] += g/J*(-1.0/3.0);	// C[1][1][2][2]

		D[3][3] += g/J*(0.5);		// C[0][1][0][1]
		D[4][4] += g/J*(0.5);		// C[1][2][1][2]
		D[5][5] += g/J*(0.5);		// C[0][2][0][2]
	}
	else if ((fabs(lam2[0] - lam2[1]) >= m_eps) && (fabs(lam2[0] - lam2[2]) >= m_eps) && (fabs(lam2[1] - lam2[2]) >= m_eps))
	{
		// all three are different
		for (i=0; i<3; ++i)
		{
			j = (i+1)%3;
			k = (j+1)%3;

			// D prime i
			Dpi = 8*lam[i]*lam[i]*lam[i]-2*I1*lam[i]-2*I3/(lam[i]*lam[i]*lam[i]);

			// D i
			Di = (lam2[i] - lam2[j])*(lam2[i] - lam2[k]);

			// the matrix mi
			ma = ((b - I*lam2[j])*(b - I*lam2[k]))/Di;

			for (j=0; j<6; ++j)
			{
				i0 = LUT[j][0];
				i1 = LUT[j][1];
				for (k=j; k<6; ++k)
				{
					j0 = LUT[k][0];
					j1 = LUT[k][1];

					BxB[j][k] = B [i0][i1]*B [j0][j1];
					BxM[j][k] = B [i0][i1]*ma(j0, j1);
					MxB[j][k] = ma(i0 ,i1)*B [j0][j1];
					IxM[j][k] = I2[i0][i1]*ma(j0, j1);
					MxI[j][k] = ma(i0 ,i1)*I2[j0][j1];
					MxM[j][k] = ma(i0, i1)*ma(j0, j1);
				}
			}

			// the tensor dmi
			for (j=0; j<6; ++j)
				for (k=j; k<6; ++k)
				{
					dm[j][k] = 1.0/Di*(Ib[j][k] - BxB[j][k] + I3/lam2[i]*(IxI[j][k] - I4[j][k])) \
						+ 1.0/Di*(lam2[i]*(BxM[j][k]+ MxB[j][k]) - 0.5*Dpi*lam[i]*MxM[j][k]) \
						- 1.0/Di*(I3/lam2[i]*(IxM[j][k] + MxI[j][k]));
				}

			beta = 0;
			for (l=0; l<MAX_TERMS; ++l) beta += m_c[l]/m_m[l]*(lamd[i][l] - (lamd[0][l] + lamd[1][l] + lamd[2][l])/3.0);


			// add everything together
			for (j=0; j<6; ++j)
			{
				for (k=j; k<6; ++k)
				{
					D[j][k] += 2.0/J*beta*dm[j][k];
				}
			}

			for (n=0; n<3; ++n)
			{
				j = (n+1)%3;
				k = (j+1)%3;

				Di = (lam2[n] - lam2[j])*(lam2[n] - lam2[k]);
				mb = ((b - I*lam2[j])*(b - I*lam2[k]))/Di;

				gab = 0;
				if (i==n)
				{
					for (l=0; l<MAX_TERMS; ++l) gab += m_c[l]*(lamd[i][l]/3.0 + (lamd[0][l]+lamd[1][l]+lamd[2][l])/9.0);
				}
				else
				{
					for (l=0; l<MAX_TERMS; ++l) gab += m_c[l]*(-lamd[i][l]/3.0 -lamd[n][l]/3.0 + (lamd[0][l]+lamd[1][l]+lamd[2][l])/9.0);
				}

				for (j=0; j<6; ++j)
				{
					i0 = LUT[j][0];
					i1 = LUT[j][1];

					for (k=j; k<6; ++k)
					{
						j0 = LUT[k][0];
						j1 = LUT[k][1];

						D[j][k] += 1.0/J*gab*ma(i0,i1)*mb(j0,j1);
					}
				}
			}
		}
	}
	else
	{
		// two are the same
		if      (fabs(lam2[0] - lam2[1]) <= m_eps) i = 2;
		else if (fabs(lam2[0] - lam2[2]) <= m_eps) i = 1;
		else i = 0;
		j = (i+1)%3;
		k = (j+1)%3;

		// D prime i
		Dpi = 8*lam[i]*lam[i]*lam[i]-2*I1*lam[i]-2*I3/(lam[i]*lam[i]*lam[i]);

		// D i
		Di = (lam2[i] - lam2[j])*(lam2[i] - lam2[k]);

		// beta i
		beta1 = beta3 = 0;
		for (l=0; l<MAX_TERMS; ++l) 
		{
			beta1 += m_c[l]/m_m[l]*(lamd[i][l] - (lamd[0][l] + lamd[1][l] + lamd[2][l])/3.0);
			beta3 += m_c[l]/m_m[l]*(lamd[j][l] - (lamd[0][l] + lamd[1][l] + lamd[2][l])/3.0);
		}

		g11 = g33 = g13 = 0;
		for (l=0; l<MAX_TERMS; ++l) g11 += m_c[l]*(lamd[i][l]/3.0 + (lamd[0][l]+lamd[1][l]+lamd[2][l])/9.0);
		for (l=0; l<MAX_TERMS; ++l) g33 += m_c[l]*(lamd[j][l]/3.0 + (lamd[0][l]+lamd[1][l]+lamd[2][l])/9.0);
		for (l=0; l<MAX_TERMS; ++l) g13 += m_c[l]*(-lamd[i][l]/3.0 - lamd[j][l]/3.0 + (lamd[0][l]+lamd[1][l]+lamd[2][l])/9.0);

		// the matrix mi
		ma = ((b - I*lam2[j])*(b - I*lam2[k]))/Di;

		for (j=0; j<6; ++j)
		{
			i0 = LUT[j][0];
			i1 = LUT[j][1];
			for (k=j; k<6; ++k)
			{
				j0 = LUT[k][0];
				j1 = LUT[k][1];

				BxB[j][k] = B[i0][i1]*B[j0][j1];
				BxM[j][k] = B[i0][i1]*ma(j0, j1);
				MxB[j][k] = ma(i0,i1)*B[j0][j1];
				MxM[j][k] = ma(i0,i1)*ma(j0, j1);
				IxM[j][k] = I2[i0][i1]*ma(j0, j1);
				MxI[j][k] = ma(i0,i1)*I2[j0][j1];
			}
		}

		// the tensor dmi
		for (j=0; j<6; ++j)
			for (k=j; k<6; ++k)
			{
				dm[j][k] = 1.0/Di*(Ib[j][k] - BxB[j][k] + I3/lam2[i]*(IxI[j][k] - I4[j][k])) \
					+ 1.0/Di*(lam2[i]*(BxM[j][k]+ MxB[j][k]) - 0.5*Dpi*lam[i]*MxM[j][k]) \
					- 1.0/Di*(I3/lam2[i]*(IxM[j][k] + MxI[j][k]));
			}

		// add everything together
		for (j=0; j<6; ++j)
			for (k=j; k<6; ++k)
			{
				D[j][k] += (-2.0/J)*beta1*I4[j][k] + (2.0/J)*(beta3 - beta1)*dm[j][k] \
					+ (g11/J)*(IxI[j][k] - IxM[j][k] - MxI[j][k] + MxM[j][k]) \
					+ (g33/J)*MxM[j][k] \
					+ (g13/J)*(MxI[j][k] - MxM[j][k] + IxM[j][k] - MxM[j][k]);
			}
	}

	// Add the volumetric contribution
	double p = pt.avgp;
	for (j=0; j<6; ++j)
		for (k=j; k<6; ++k)
			D[j][k] += p*(IxI[j][k] - 2.0*I4[j][k]);

	// set symmetric components 
	for (j=0; j<6; ++j)
		for (k=j+1; k<6; ++k)
			D[k][j] = D[j][k];
}
