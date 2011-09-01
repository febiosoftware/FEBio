// mat3d.cpp: implementation of the mat3d class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "mat3d.h"

#define ROTATE(a, i, j, k, l) g=a[i][j]; h=a[k][l];a[i][j]=g-s*(h+g*tau); a[k][l] = h + s*(g - h*tau);

void mat3ds::eigen(double l[3], vec3d r[3])
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
