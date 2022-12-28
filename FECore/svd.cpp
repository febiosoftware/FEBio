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
#include "matrix.h"
#include <math.h>
using namespace std;

#define SQR(a) ((a)*(a))
#define FMAX(a, b) ((a)>(b)?(a):(b))
#define IMIN(a, b) ((a)<(b)?(a):(b))
#define SIGN(a, b) ((b) > 0.0 ? fabs(a) : -fabs(a))

void svbksb(matrix& u, vector<double>& w, matrix& v, vector<double>& b, vector<double>& x)
{
	int jj, j, i;
	double s;

	int m = u.rows();
	int n = u.columns();

	vector<double> tmp(n);
	for (j=0; j<n; ++j)
	{
		s = 0.0;
		if (w[j])
		{
			for (i=0; i<m; ++i) s += u[i][j]*b[i];
			s /= w[j];
		}
		tmp[j] = s;
	}
	for (j=0; j<n; ++j)
	{
		s = 0.0;
		for (jj=0; jj<n; ++jj) s += v[j][jj]*tmp[jj];
		x[j] = s;
	}
}

double pythag(double a, double b)
{
	double absa, absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0 + SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0 + SQR(absa/absb)));
}

void svdcmp(matrix& a, vector<double>& w, matrix& v)
{
	int flag, i, its, j, jj, k, l, nm;
	double anorm, c, f, g, h, s, scale, x, y, z;

	int m = a.rows();
	int n = a.columns();

	vector<double> rv1(n);
	g = scale = anorm = 0.0;
	for (i=0; i<n; ++i)
	{
		l=i+1;																// diff between c++/c NR
		rv1[i] = scale*g;
		g = s = scale = 0.0;
		if (i < m)
		{
			for (k=i; k<m; ++k) scale += fabs(a[k][i]);
			if (scale != 0.0)
			{
				for (k=i; k<m; ++k) 
				{
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f = a[i][i];
				g = -SIGN(sqrt(s), f);
				h = f*g-s;
				a[i][i] = f-g;
				for (j=l; j<n; ++j)										// diff between c++/c NR
				{
					for (s=0.0, k=i; k<m; k++) s += a[k][i]*a[k][j];
					f = s/h;
					for (k=i; k<m; ++k) a[k][j] += f*a[k][i];
				}
				for (k=i; k<m; k++) a[k][i] *= scale;
			}
		}
		w[i] = scale*g;
		g=s=scale = 0.0;
		if ((i+1<=m) && (i!=n-1))									// diff between c++/c NR
		{
			for (k=l; k<n; ++k) scale += fabs(a[i][k]);		// diff between c++/c NR
			if (scale != 0.0)
			{
				for (k=l; k<n; ++k)							// diff between c++/c NR
				{
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f = a[i][l];									// diff between c++/c NR
				g = -SIGN(sqrt(s), f);
				h=f*g-s;
				a[i][l] = f - g;
				for (k=l; k<n; k++) rv1[k] = a[i][k]/h;		// diff between c++/c NR
				for (j=l; j<m; ++j)							// diff between c++/c NR
				{
					for (s=0.0, k=l; k<n; k++) s += a[j][k]*a[i][k];		// diff between c++/c NR
					for (k=l; k<n; k++) a[j][k] += s*rv1[k];				// diff between c++/c NR
				}
				for (k=l; k<n; k++) a[i][k] *= scale;			// diff between c++/c NR
			}
		}
		anorm = FMAX(anorm, (fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n-1; i>=0; --i) 
	{
		if (i<n-1)
		{
			if (g != 0.0)
			{
				for (j=l; j<n; ++j)
					v[j][i] = (a[i][j]/a[i][l])/g;
				for (j=l; j<n; ++j)
				{
					for (s=0.0, k=l; k<n; ++k) s += a[i][k]*v[k][j];
					for (k=l; k<n; k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l; j<n; j++) v[i][j] = v[j][i] = 0.0;
		}
		v[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}
	for (i=IMIN(m,n)-1; i>=0; i--)
	{
		l=i+1;
		g=w[i];
		for (j=l; j<n; ++j) a[i][j] = 0.0;
		if (g != 0.0)
		{
			g = 1.0/g;
			for (j=l; j<n; ++j)
			{
				for (s=0.0, k=l; k<m; k++) s += a[k][i]*a[k][j];
				f = (s/a[i][i])*g;
				for (k=i; k<m; k++) a[k][j] += f*a[k][i];
			}
			for (j=i; j<m; j++) a[j][i] *= g;
		}
		else for (j=i; j<m; j++) a[j][i] = 0.0;
		++a[i][i];
	}
	for (k=n-1; k>=0; k--)
	{
		for (its=0; its<30; ++its)
		{
			flag = 1;
			for (l=k; l>=0; --l)
			{
				nm = l-1;
				if ((fabs(rv1[l])+anorm) == anorm)
				{
					flag = 0;
					break;
				}
				if ((fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag == 1)
			{
				c = 0.0;
				s = 1.0;
				for (i=l; i<=k; ++i)		// diff between c++/c NR
				{
					f = s*rv1[i];
					rv1[i] = c*rv1[i];
					if ((fabs(f)+anorm) == anorm) break;
					g = w[i];
					h = pythag(f, g);
					w[i] = h;
					h = 1.0/h;
					c = g*h;
					s= -f*h;
					for (j=0; j<m; j++)
					{
						y=a[j][nm];
						z=a[j][i];
						a[j][nm] = y*c+z*s;
						a[j][i] = z*c - y*s;
					}
				}
			}
			z = w[k];
			if (l == k)
			{
				if (z < 0.0)
				{
					w[k] = -z;
					for (j=0; j<n; j++) v[j][k] = -v[j][k];
				}
				break;
			}
			// TODO: do something drastic when its hits 29
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g = pythag(f, 1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g, f)))-h)) / x;
			c = s = 1.0;
			for (j=l; j<=nm; j++)
			{
				i = j+1;
				g = rv1[i];
				y = w[i];
				h = s*g;
				g = c*g;
				z = pythag(f, h);
				rv1[j] = z;
				c = f/z;
				s = h/z;
				f = x*c + g*s;
				g = g*c - x*s;
				h = y*s;
				y *= c;
				for (jj=0; jj<n; ++jj)
				{
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j] = x*c + z*s;
					v[jj][i] = z*c - x*s;
				}
				z = pythag(f, h);
				w[j] = z;
				if (z)
				{
					z = 1.0/z;
					c = f*z;
					s = h*z;
				}
				f = c*g + s*y;
				x = c*y - s*g;
				for (jj=0; jj<m; ++jj)
				{
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y*c + z*s;
					a[jj][i] = z*c - y*s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
}
