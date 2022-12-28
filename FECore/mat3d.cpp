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
#include "mat3d.h"
#include "eig3.h"
#include "sys.h"

#define ROTATE(a, i, j, k, l) g=a[i][j]; h=a[k][l];a[i][j]=g-s*(h+g*tau); a[k][l] = h + s*(g - h*tau);

void mat3ds::eigen(double l[3], vec3d r[3]) const
{
	const int NMAX = 50;
	double sm, tresh, g, h, t, c, tau, s, th;
	int i, j, k;

	// copy the matrix components since we will be overwriting them
	double a[3][3] = {
			{m[XX], m[XY], m[XZ]},
			{m[XY], m[YY], m[YZ]}, 
			{m[XZ], m[YZ], m[ZZ]}
	};

	// the v matrix contains the eigen vectors
	// intialize to identity
	double v[3][3] = {
		{ 1, 0, 0 },
		{ 0, 1, 0 },
		{ 0, 0, 1 }
	};

	// initialize b and d to the diagonal of a
	double b[3] = {a[0][0], a[1][1], a[2][2]};
	double d[3] = {a[0][0], a[1][1], a[2][2]};
	double z[3] = {0};

	const double eps = 0;//1.0e-15;

	// loop
	int n, nrot = 0;
	for (n=0; n<NMAX; ++n)
	{
		// sum off-diagonal elements
		sm = fabs(a[0][1]) + fabs(a[0][2]) + fabs(a[1][2]);
		if (sm <= eps) break;

		// set the treshold
		if (n < 3) tresh = 0.2*sm/9.0; else tresh = 0.0;

		// loop over off-diagonal elements
		for (i=0; i<2; ++i)
		{
			for (j=i+1; j<3; ++j)
			{
				g = 100.0*fabs(a[i][j]);

				// after four sweeps, skip the rotation if the off-diagonal element is small
				if ((n > 3) && ((fabs(d[i])+g) == fabs(d[i]))
							&& ((fabs(d[j])+g) == fabs(d[j])))
				{
					a[i][j] = 0.0;
				}
				else if (fabs(a[i][j]) > tresh)
				{
					h = d[j] - d[i];
					if ((fabs(h)+g) == fabs(h))
						t = a[i][j]/h;
					else
					{
						th = 0.5*h/a[i][j];
						t = 1.0/(fabs(th) + sqrt(1+th*th));
						if (th < 0.0) t = -t;
					}

					c = 1.0/sqrt(1.0 + t*t);
					s = t*c;
					tau = s/(1.0+c);
					h = t*a[i][j];
					z[i] -= h;
					z[j] += h;
					d[i] -= h;
					d[j] += h;
					a[i][j] = 0;

					for (k=  0; k<=i-1; ++k) { ROTATE(a, k, i, k, j) }
					for (k=i+1; k<=j-1; ++k) { ROTATE(a, i, k, k, j) }
					for (k=j+1; k<   3; ++k) { ROTATE(a, i, k, j, k) }
					for (k=  0; k<   3; ++k) { ROTATE(v, k, i, k, j) }
					++nrot;
				}
			}
		}

		for (i=0; i<3; ++i) 
		{
			b[i] += z[i];
			d[i] = b[i];
			z[i] = 0.0;
		}
	}

	// we sure we converged
	assert(n < NMAX);

	// copy eigenvalues
	l[0] = d[0];
	l[1] = d[1];
	l[2] = d[2];

	// copy eigenvectors
	if (r)
	{
		r[0].x = v[0][0]; r[0].y = v[1][0]; r[0].z = v[2][0];
		r[1].x = v[0][1]; r[1].y = v[1][1]; r[1].z = v[2][1];
		r[2].x = v[0][2]; r[2].y = v[1][2]; r[2].z = v[2][2];
	}
}

//-----------------------------------------------------------------------------
// Calculate the eigenvalues of A using an analytical expression for the 
// eigen values.
void mat3ds::exact_eigen(double l[3]) const
{
	const double S3 = sqrt(3.0);
	const double S2 = sqrt(2.0);

	mat3ds S = dev();
	double nS = S.norm();
	mat3ds T = (S.sqr()).dev();
	double nT = T.norm();

	double D = nS * nT;
	if (D > 0.0)
	{
		double w = S.dotdot(T) / D;	if (w > 1.0) w = 1.0; if (w < -1.0) w = -1.0;
		double t = asin(w) / 3.0;
		double r = S.norm();
		double z = tr() / S3;

		l[0] = z / S3 + (r / S2)*(sin(t) / S3 + cos(t));
		l[1] = z / S3 - (S2 / S3)*r*sin(t);
		l[2] = z / S3 + (r / S2)*(sin(t) / S3 - cos(t));
	}
	else
	{
		l[0] = l[1] = l[2] = 0.0;
	}
}

//-----------------------------------------------------------------------------
// Calculate the eigenvalues and eigenvectors of A using the method of
// Connelly Barnes ( http://barnesc.blogspot.com/2007/02/eigenvectors-of-3x3-symmetric-matrix.html )
void mat3ds::eigen2(double l[3], vec3d r[3]) const
{
    double A[3][3] = {xx(), xy(), xz(), xy(), yy(), yz(), xz(), yz(), zz()};
    double V[3][3];
    if (ISNAN(tr())) return;
    eigen_decomposition(A, V, l);
    if (r) {
        r[0] = vec3d(V[0][0],V[1][0],V[2][0]);
        r[1] = vec3d(V[0][1],V[1][1],V[2][1]);
        r[2] = vec3d(V[0][2],V[1][2],V[2][2]);
    }
}

//-----------------------------------------------------------------------------
// calculates the unique right polar decomposition F = R*U
void mat3d::right_polar(mat3d& R, mat3ds& U) const
{
	const mat3d& F = *this;
	mat3ds C = (F.transpose()*F).sym();

	double l[3];
	vec3d v[3];
	C.eigen2(l, v);

	U = dyad(v[0])*sqrt(l[0]) + dyad(v[1])*sqrt(l[1]) + dyad(v[2])*sqrt(l[2]);
	R = F*U.inverse();
}

//-----------------------------------------------------------------------------
// calculates the unique left polar decomposition F = V*R
void mat3d::left_polar(mat3ds& V, mat3d& R) const
{
	const mat3d& F = *this;
	mat3ds b = (F*F.transpose()).sym();

	double l[3];
	vec3d v[3];
	b.eigen2(l, v);

	V = dyad(v[0])*sqrt(l[0]) + dyad(v[1])*sqrt(l[1]) + dyad(v[2])*sqrt(l[2]);
	R = V.inverse()*F;
}

//-----------------------------------------------------------------------------
// the "max shear" value
double mat3ds::max_shear() const
{
	double e[3];
	exact_eigen(e);

	double t1 = fabs(0.5f*(e[1] - e[2]));
	double t2 = fabs(0.5f*(e[2] - e[0]));
	double t3 = fabs(0.5f*(e[0] - e[1]));

	// TODO: is this necessary? I think the values returned
	//       by Principals are already ordered.
	double tmax = t1;
	if (t2 > tmax) tmax = t2;
	if (t3 > tmax) tmax = t3;

	return tmax;
}


//-----------------------------------------------------------------------------
void mat3fs::Principals(float e[3]) const
{
	const static float ONETHIRD = 1.f / 3.f;

	// pressure
	float p = -(x + y + z) * ONETHIRD;

	DeviatoricPrincipals(e);

	e[0] -= p;
	e[1] -= p;
	e[2] -= p;
}

//-----------------------------------------------------------------------------
void mat3fs::DeviatoricPrincipals(float e[3]) const
{
	const static float ONETHIRD = 1.f / 3.f;

	// pressure
	float p = -(x + y + z) * ONETHIRD;

	// deviatoric stresses
	float dev[3];
	dev[0] = x + p;
	dev[1] = y + p;
	dev[2] = z + p;

	// invariants
	float I[3];
	I[0] = dev[0] + dev[1] + dev[2]; // = tr[s']

	I[1] = 0.5f * (dev[0] * dev[0]
		+ dev[1] * dev[1]
		+ dev[2] * dev[2])
		+ xy * xy
		+ yz * yz
		+ xz * xz; // = s':s'

	I[2] = -dev[0] * dev[1] * dev[2] - 2.0f * xy * yz * xz
		+ dev[0] * yz * yz
		+ dev[1] * xz * xz
		+ dev[2] * xy * xy; // = -det(s')

						  // check to see if we can have a non-zero divisor, if no
						  // set principal stress to 0
	if (I[1] != 0)
	{
		double a = -0.5 * sqrt(27.0 / I[1]) * I[2] / I[1];
		if (a < 0)
			a = MAX(a, -1.0);
		else
			a = MIN(a, 1.0);
		double w = acos(a) * ONETHIRD;
		double val = 2.0 * sqrt(I[1] * ONETHIRD);
		e[0] = (float)(val * cos(w));
		w = w - 2.0 * PI * ONETHIRD;
		e[1] = (float)(val * cos(w));
		w = w + 4.0 * PI * ONETHIRD;
		e[2] = (float)(val * cos(w));
	}
	else
	{
		e[0] = e[1] = e[2] = 0.f;
	}
}

//-----------------------------------------------------------------------------
float mat3fs::MaxShear() const
{
	float e[3];
	Principals(e);

	float t1 = (float)fabs(0.5f * (e[1] - e[2]));
	float t2 = (float)fabs(0.5f * (e[2] - e[0]));
	float t3 = (float)fabs(0.5f * (e[0] - e[1]));

	// TODO: is this necessary? I think the values returned
	//       by Principals are already ordered.
	float tmax = t1;
	if (t2 > tmax) tmax = t2;
	if (t3 > tmax) tmax = t3;

	return tmax;
}

//-----------------------------------------------------------------------------
#define ROTATE(a, i, j, k, l) g=a[i][j]; h=a[k][l];a[i][j]=g-s*(h+g*tau); a[k][l] = h + s*(g - h*tau);
#define SWAPF(a, b) { float t = a; a = b; b = t; }
#define SWAPV(a, b) { vec3f t = a; a = b; b = t; }

void mat3fs::eigen(vec3f e[3], float l[3]) const
{
	const int NMAX = 50;
	double sm, tresh, g, h, t, c, tau, s, th;
	int i, j, k;

	// copy the Matrix components since we will be overwriting them
	double a[3][3] = {
		{ x , xy, xz },
		{ xy, y , yz },
		{ xz, yz, z }
	};

	// the v Matrix contains the eigen vectors
	// intialize to identity
	double v[3][3] = {
		{ 1, 0, 0 },
		{ 0, 1, 0 },
		{ 0, 0, 1 }
	};

	// initialize b and d to the diagonal of a
	double b[3] = { a[0][0], a[1][1], a[2][2] };
	double d[3] = { a[0][0], a[1][1], a[2][2] };
	double z[3] = { 0 };

	const double eps = 0;//1.0e-15;

						 // loop
	int n, nrot = 0;
	for (n = 0; n < NMAX; ++n)
	{
		// sum off-diagonal elements
		sm = fabs(a[0][1]) + fabs(a[0][2]) + fabs(a[1][2]);
		if (sm <= eps) break;

		// set the treshold
		if (n < 3) tresh = 0.2 * sm / 9.0; else tresh = 0.0;

		// loop over off-diagonal elements
		for (i = 0; i < 2; ++i)
		{
			for (j = i + 1; j < 3; ++j)
			{
				g = 100.0 * fabs(a[i][j]);

				// after four sweeps, skip the rotation if the off-diagonal element is small
				if ((n > 3) && ((fabs(d[i]) + g) == fabs(d[i]))
					&& ((fabs(d[j]) + g) == fabs(d[j])))
				{
					a[i][j] = 0.0;
				}
				else if (fabs(a[i][j]) > tresh)
				{
					h = d[j] - d[i];
					if ((fabs(h) + g) == fabs(h))
						t = a[i][j] / h;
					else
					{
						th = 0.5 * h / a[i][j];
						t = 1.0 / (fabs(th) + sqrt(1 + th * th));
						if (th < 0.0) t = -t;
					}

					c = 1.0 / sqrt(1.0 + t * t);
					s = t * c;
					tau = s / (1.0 + c);
					h = t * a[i][j];
					z[i] -= h;
					z[j] += h;
					d[i] -= h;
					d[j] += h;
					a[i][j] = 0;

					for (k = 0; k <= i - 1; ++k) { ROTATE(a, k, i, k, j) }
					for (k = i + 1; k <= j - 1; ++k) { ROTATE(a, i, k, k, j) }
					for (k = j + 1; k < 3; ++k) { ROTATE(a, i, k, j, k) }
					for (k = 0; k < 3; ++k) { ROTATE(v, k, i, k, j) }
					++nrot;
				}
			}
		}

		for (i = 0; i < 3; ++i)
		{
			b[i] += z[i];
			d[i] = b[i];
			z[i] = 0.0;
		}
	}

	// we sure we converged
	assert(n < NMAX);

	// copy eigenvalues
	l[0] = (float)d[0];
	l[1] = (float)d[1];
	l[2] = (float)d[2];

	// copy eigenvectors
	e[0].x = (float)v[0][0]; e[0].y = (float)v[1][0]; e[0].z = (float)v[2][0];
	e[1].x = (float)v[0][1]; e[1].y = (float)v[1][1]; e[1].z = (float)v[2][1];
	e[2].x = (float)v[0][2]; e[2].y = (float)v[1][2]; e[2].z = (float)v[2][2];

	// we still need to sort the eigenvalues
	if (l[1] > l[0]) { SWAPF(l[0], l[1]); SWAPV(e[0], e[1]); }
	if (l[2] > l[0]) { SWAPF(l[0], l[2]); SWAPV(e[0], e[2]); }
	if (l[2] > l[1]) { SWAPF(l[2], l[1]); SWAPV(e[2], e[1]); }
}

//-----------------------------------------------------------------------------
vec3f mat3fs::PrincDirection(int l)
{
	vec3f e[3];
	float lam[3];
	eigen(e, lam);
	return e[l] * lam[l];
}

//-----------------------------------------------------------------------------
double fractional_anisotropy(const mat3fs& m)
{
	vec3f e[3];
	float l[3];
	m.eigen(e, l);

	double la = (l[0] + l[1] + l[2]) / 3.0;
	double D = sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
	double fa = 0.0;
	if (D != 0) fa = sqrt(3.0 / 2.0) * sqrt((l[0] - la) * (l[0] - la) + (l[1] - la) * (l[1] - la) + (l[2] - la) * (l[2] - la)) / D;

	return fa;
}
