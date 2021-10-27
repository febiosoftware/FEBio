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
#include "tools.h"
#include <math.h>
#include <limits>
#include <assert.h>
#include <float.h>
#include "matrix.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

#define FMAX(a,b) ((a)>(b) ? (a) : (b))

//=============================================================================
// powell's method
// from Numerical Recipes in C, section 10.5, page 417-419
// modified for this application
//
void powell(double* p, double* xi, int n, double ftol, int* iter, double* fret, double(*fnc)(double[]))
{
	int i, j, ibig;
	double fp, fptt, del, t;

	const int ITMAX = 200;
	const double TINY = 1.0e-25;

	double* pt = new double[n];
	double* ptt = new double[n];
	double* xit = new double[n];

	*fret = (*fnc)(p);
	for (i = 0; i<n; i++) pt[i] = p[i];				// save the initial point
	for (*iter = 1;; ++(*iter))
	{
		fp = *fret;
		ibig = 0;
		del = 0.0;		// will be biggest function decrease

		for (i = 0; i<n; i++)	// in each iteration loop over all directions in the set
		{
			for (j = 0; j<n; j++) xit[j] = xi[i*n + j];	// copy the direction
			fptt = (*fret);
			linmin(p, xit, n, fret, fnc);		// minimize along it
			if (fptt - (*fret) > del)			// and record it as the largets decrease so far 	
			{
				del = fptt - (*fret);
				ibig = i;
			}
		}

		// termination criterion
		if (2.0*(fp - (*fret)) <= ftol*(fabs(fp) + fabs(*fret)) + TINY)
		{
			delete[] pt;
			delete[] ptt;
			delete[] xit;
			return;
		}

		// check we are not exceeding max number of iterations
		if (*iter == ITMAX)
		{
			printf("FATAL ERROR : Max iterations reached\n");
			exit(0);
		}

		for (j = 0; j<n; j++)
		{
			ptt[j] = 2.0*p[j] - pt[j];
			xit[j] = p[j] - pt[j];
			pt[j] = p[j];
		}

		fptt = (*fnc)(ptt);
		if (fptt < fp)
		{
			t = 2.0*(fp - 2.0*(*fret) + fptt)*SQR(fp - (*fret) - del) - del*SQR(fp - fptt);
			if (t < 0.0)
			{
				linmin(p, xit, n, fret, fnc);
				for (j = 0; j<n; j++)
				{
					xi[ibig *n + j] = xi[(n - 1)*n + j];
					xi[(n - 1)*n + j] = xit[j];
				}
			}
		}
	}
}

//--------------------------------------------------------------------------------
// line minimization routine
// from Numerical Recipes in C, section 10.5, page 419-420
// modified for this application
//
int ncom;
double *pcom, *xicom, *xt;
double(*nrfunc)(double[]);

double f1dim(double x)
{
	int j;
	double f;
	for (j = 0; j<ncom; j++) xt[j] = pcom[j] + x*xicom[j];
	f = (*nrfunc)(xt);
	return f;
}

void linmin(double* p, double* xi, int n, double* fret, double(*fnc)(double[]))
{
	int j;
	double xx, xmin, fx, fb, fa, bx, ax;
	const double TOL = 2.0e-4;

	ncom = n;
	pcom = new double[n];
	xicom = new double[n];
	xt = new double[n];

	nrfunc = fnc;
	for (j = 0; j<n; j++)
	{
		pcom[j] = p[j];
		xicom[j] = xi[j];
	}
	ax = 0.0;
	xx = 1.0e-4;
	mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, f1dim);
	*fret = brent(ax, xx, bx, f1dim, TOL, &xmin);
	for (j = 0; j<n; j++)
	{
		xi[j] *= xmin;
		p[j] += xi[j];
	}

	delete[] xt;
	delete[] pcom;
	delete[] xicom;
}


//--------------------------------------------------------------------------------
// routine to find 1D minimum
// from Numerical Recipes in C, section 10.3, page 404-405
// modified for this application
//
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b)>=0.0 ? (a) : (-(a)))

double brent(double ax, double bx, double cx, double(*f)(double), double tol, double* xmin)
{
	const int ITMAX = 100;
	const double CGOLD = 0.3819660;
	const double ZEPS = 1.0e-10;
	int iter;
	double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
	double e = 0.0;

	a = (ax < cx ? ax : cx);
	b = (ax > cx ? ax : cx);
	x = w = v = bx;
	fw = fv = fx = (*f)(x);
	for (iter = 1; iter <= ITMAX; iter++)
	{
		xm = 0.5*(a + b);
		tol2 = 2.0*(tol1 = tol*fabs(x) + ZEPS);
		if (fabs(x - xm) <= (tol2 - 0.5*(b - a)))
		{
			*xmin = x;
			return fx;
		}
		if (fabs(e) > tol1)
		{
			r = (x - w)*(fx - fv);
			q = (x - v)*(fx - fw);
			p = (x - v)*q - (x - w)*r;
			q = 2.0*(q - r);
			if (q > 0.0) p = -p;
			q = fabs(q);
			etemp = e;
			e = d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a - x) || p >= q*(b - x))
				d = CGOLD*(e = (x >= xm ? a - x : b - x));
			else
			{
				d = p / q;
				u = x + d;
				if (u - a < tol2 || b - u < tol2) d = SIGN(tol1, xm - x);
			}
		}
		else d = CGOLD*(e = (x >= xm ? a - x : b - x));
		u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
		fu = (*f)(u);

		if (fu <= fx)
		{
			if (u >= x) a = x; else b = x;
			SHFT(v, w, x, u);
			SHFT(fv, fw, fx, fu);
		}
		else
		{
			if (u<x) a = u; else b = u;
			if (fu <= fw || w == x)
			{
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			}
			else if (fu <= fv || v == x || v == w)
			{
				v = u;
				fv = fu;
			}
		}
	}

	fprintf(stderr, "ERROR : Too many iterations in brent routine\n");

	*xmin = x;
	return fx;
}


#define SHFT2(a, b, c) (a)=(b);(b)=(c);
#define SHFT(a, b, c, d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN2(a, b) ((b)>=0?fabs(a):(-fabs(a)))

void mnbrak(double* pa, double* pb, double* pc, double* pfa, double* pfb, double* pfc, double(*func)(double))
{
	const double GOLD = 1.618034;
	const double TINY = 1.0e-20;
	const double GLIMIT = 100;

	double& a = *pa;
	double& b = *pb;
	double& c = *pc;

	double& fa = *pfa;
	double& fb = *pfb;
	double& fc = *pfc;

	double ulim, u, r, q, fu, dum;

	fa = func(a);
	fb = func(b);
	if (fb>fa)
	{
		SHFT(dum, a, b, dum);
		SHFT(dum, fb, fa, dum);
	}

	c = b + GOLD*(b - a);
	fc = func(c);
	while (fb > fc)
	{
		r = (b - a)*(fb - fc);
		q = (b - c)*(fb - fa);
		u = b - ((b - c)*q - (b - a)*r) / (2.0*SIGN2(FMAX(fabs(q - r), TINY), q - r));

		ulim = b + GLIMIT*(c - b);

		if ((b - u)*(u - c) > 0)
		{
			fu = func(u);
			if (fu < fc)
			{
				a = b;
				b = u;
				fa = fb;
				fb = fu;
				return;
			}
			else if (fu > fb)
			{
				c = u;
				fc = fu;
				return;
			}

			u = c + GOLD*(c - b);
			fu = func(u);
		}
		else if ((c - u)*(u - ulim) > 0)
		{
			fu = func(u);
			if (fu < fc)
			{
				SHFT(b, c, u, c + GOLD*(c - b));
				SHFT(fb, fc, fu, func(u));
			}
		}
		else if ((u - ulim)*(ulim - c) >= 0)
		{
			u = ulim;
			fu = func(u);
		}
		else
		{
			u = c + GOLD*(c - b);
			fu = func(u);
		}

		SHFT(a, b, c, u);
		SHFT(fa, fb, fc, fu);
	}
}

double golden(double ax, double bx, double cx, double(*f)(double), double tol, double* xmin)
{
	const double R = 0.61803399;
	const double C = 1 - R;

	double f1, f2, x0, x1, x2, x3;

	x0 = ax;
	x3 = cx;

	if (fabs(cx - bx) > fabs(bx - ax))
	{
		x1 = bx;
		x2 = bx + C*(cx - bx);
	}
	else
	{
		x2 = bx;
		x1 = bx - C*(bx - ax);
	}
	f1 = f(x1);
	f2 = f(x2);

	while (fabs(x3 - x0) > tol*(fabs(x1) + fabs(x2)))
	{
		if (f2 < f1)
		{
			SHFT(x0, x1, x2, R*x1 + C*x3);
			SHFT2(f1, f2, f(x2));
		}
		else
		{
			SHFT(x3, x2, x1, R*x2 + C*x0);
			SHFT2(f2, f1, f(x1));
		}
	}

	if (f1 < f2)
	{
		*xmin = x1;
		return f1;
	}
	else
	{
		*xmin = x2;
		return f2;
	}
}

//-----------------------------------------------------------------------------
double zbrent(double func(double, void*), double x1, double x2, double tol, void* data)
{
	const int ITMAX = 100;
	const double EPS = DBL_EPSILON;
	int iter;
	double a=x1, b=x2, c=x2, d, e, min1, min2;
	double fa=func(a, data), fb = func(b, data), fc, p, q, r, s, tol1, xm;

	if ((fa > 0.0 && fb > 00) || (fa < 0.0 && fb < 0.0))
	{
		assert(false);
		return 0.0;
	}

	fc = fb;
	for (iter = 0; iter<ITMAX; iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c = a;
			fc = fa;
			e = d = b - a;
		}
		if (fabs(fc) < fabs(fb)) {
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}
		tol1 = 2.0*EPS*fabs(b) + 0.5*tol;
		xm = 0.5*(c - b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s = fb / fa;
			if (a == c) {
				p = 2.0*xm*s;
				q = 1.0 - s;
			} else {
				q = fa/fc;
				r = fb/fc;
				p = s*(2.0*xm*q*(q-r)-(b-a)*(r - 1.0));
				q = (q - 1.0)*(r - 1.0)*(s - 1.0);
			}
			if (p > 0.0) q = -q;
			p = fabs(p);
			min1 = 3.0*xm*q - fabs(tol1*q);
			min2 = fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e = d;
				d = p/q;
			} else {
				d = xm;
				e = d;
			}
		} else {
			d = xm;
			e = d;
		}
		a = b;
		fa = fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1, xm);
		fb = func(b, data);
	}
	assert(false);
	return 0.0;
}

//-----------------------------------------------------------------------------
bool zbrac(double f(double, void*), double& x1, double& x2, void* data)
{
	const int MAXTRY = 50;
	const double FACTOR = 1.6;

	if (x1 == x2)
	{
		assert(false);
		return false;
	}

	double f1 = f(x1, data);
	double f2 = f(x2, data);
	for (int j=0; j<MAXTRY; ++j)
	{
		if (f1*f2 < 0.0) return true;
		if (fabs(f1) < fabs(f2))
			f1 = f(x1 += FACTOR*(x1 - x2), data);
		else
			f2 = f(x2 += FACTOR*(x2 - x1), data);
	}
	return false;
}

//-----------------------------------------------------------------------------
void solve_3x3(double A[3][3], double b[3], double x[3])
{
	double D = A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] + A[1][0] * A[2][1] * A[0][2] \
		- A[1][1] * A[2][0] * A[0][2] - A[2][2] * A[1][0] * A[0][1] - A[0][0] * A[2][1] * A[1][2];

	assert(D != 0);

	double Ai[3][3];
	Ai[0][0] = A[1][1] * A[2][2] - A[2][1] * A[1][2];
	Ai[0][1] = A[2][1] * A[0][2] - A[0][1] * A[2][2];
	Ai[0][2] = A[0][1] * A[1][2] - A[1][1] * A[0][2];

	Ai[1][0] = A[2][0] * A[1][2] - A[1][0] * A[2][2];
	Ai[1][1] = A[0][0] * A[2][2] - A[2][0] * A[0][2];
	Ai[1][2] = A[1][0] * A[0][2] - A[0][0] * A[1][2];

	Ai[2][0] = A[1][0] * A[2][1] - A[2][0] * A[1][1];
	Ai[2][1] = A[2][0] * A[0][1] - A[0][0] * A[2][1];
	Ai[2][2] = A[0][0] * A[1][1] - A[0][1] * A[1][0];

	x[0] = (Ai[0][0] * b[0] + Ai[0][1] * b[1] + Ai[0][2] * b[2]) / D;
	x[1] = (Ai[1][0] * b[0] + Ai[1][1] * b[1] + Ai[1][2] * b[2]) / D;
	x[2] = (Ai[2][0] * b[0] + Ai[2][1] * b[1] + Ai[2][2] * b[2]) / D;


#ifdef _DEBUG
	double r[3];
	r[0] = b[0] - (A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2]);
	r[1] = b[1] - (A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2]);
	r[2] = b[2] - (A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2]);

	double nr = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
#endif
}

//=============================================================================

bool LinearRegression(const std::vector<std::pair<double, double> >& data, std::pair<double, double>& res)
{
	res.first = 0.0;
	res.second = 0.0;

	int n = (int)data.size();
	if (n == 0) return false;

	double mx = 0.0, my = 0.0;
	double sxx = 0.0, sxy = 0.0;
	for (int i = 0; i < n; ++i)
	{
		double xi = data[i].first;
		double yi = data[i].second;
		mx += xi;
		my += yi;

		sxx += xi * xi;
		sxy += xi * yi;
	}
	mx /= (double)n;
	my /= (double)n;
	sxx /= (double)n;
	sxy /= (double)n;

	double D = sxx - mx * mx;
	if (D == 0.0) return false;

	double a = (sxy - mx * my) / D;
	double b = my - a * mx;

	res.first = a;
	res.second = b;

	return true;
}

class Func
{
public:
	Func() {}
	virtual ~Func() {}
	virtual void setParams(const std::vector<double>& v) = 0;
	virtual double value(double x) = 0;
	virtual double derive1(double x, int n) = 0;
	virtual double derive2(double x, int n1, int n2) = 0;
};

class Quadratic : public Func
{
public:
	Quadratic() : m_a(0.0), m_b(0.0), m_c(0.0) {}
	void setParams(const std::vector<double>& v) override { m_a = v[0]; m_b = v[1]; m_c = v[2]; }
	double value(double x) override { return m_a * x * x + m_b * x + m_c; }
	double derive1(double x, int n) override
	{
		switch (n)
		{
		case 0: return x * x; break;
		case 1: return x; break;
		case 2: return 1; break;
		default:
			assert(false);
			return 0.0;
		}
	}

	double derive2(double x, int n1, int n2) override
	{
		return 0.0;
	}

private:
	double	m_a, m_b, m_c;
};

class Exponential : public Func
{
public:
	Exponential() : m_a(0.0), m_b(0.0) {}
	void setParams(const std::vector<double>& v) override { m_a = v[0]; m_b = v[1]; }
	double value(double x) override { return m_a * exp(x * m_b); }
	double derive1(double x, int n) override
	{
		switch (n)
		{
		case 0: return exp(x * m_b); break;
		case 1: return m_a * x * exp(x * m_b); break;
		default:
			assert(false);
			return 0.0;
		}
	}

	double derive2(double x, int n1, int n2) override
	{
		if ((n1 == 0) && (n2 == 0)) return 0;
		else if ((n1 == 0) && (n2 == 1)) return x * exp(x * m_b);
		else if ((n1 == 1) && (n2 == 0)) return x * exp(x * m_b);
		else if ((n1 == 1) && (n2 == 1)) return m_a * x * x * exp(x * m_b);
		else return 0.0;
	}

private:
	double	m_a, m_b;
};

bool NonlinearRegression(const std::vector<std::pair<double, double> >& data, std::vector<double>& res, int func)
{
	int MAX_ITER = 10;
	int niter = 0;

	int n = (int)data.size();
	int m = (int)res.size();

	Func* f = 0;
	switch (func)
	{
	case 1: f = new Quadratic; break;
	case 2: f = new Exponential; break;
	}
	if (f == 0) return false;

	std::vector<double> R(m, 0.0), da(m, 0.0);
	matrix K(m, m); K.zero();

	const double absTol = 1e-15;
	const double relTol = 1e-3;
	double norm0 = 0.0;
	do
	{
		f->setParams(res);

		// evaluate residual (and norm)
		double norm = 0.0;
		for (int i = 0; i < m; ++i)
		{
			R[i] = 0.0;
			for (int j = 0; j < n; ++j)
			{
				double xj = data[j].first;
				double yj = data[j].second;
				double fj = f->value(xj);
				double Dfi = f->derive1(xj, i);
				R[i] -= (fj - yj) * Dfi;
			}

			norm += R[i] * R[i];
		}
		norm = sqrt(norm / n);

		if (norm < absTol) break;

		if (niter == 0) norm0 = norm;
		else
		{
			double rel = norm / norm0;
			if (rel < relTol) break;
		}

		// evaluate Jacobian
		for (int i = 0; i < m; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				double Kij = 0.0;
				for (int k = 0; k < n; ++k)
				{
					double xk = data[k].first;
					double yk = data[k].second;
					double fk = f->value(xk);

					double Dfi = f->derive1(xk, i);
					double Dfj = f->derive1(xk, j);

					double Dfij = f->derive2(xk, i, j);

					Kij += Dfi * Dfj + (fk - yk) * Dfij;
				}

				K[i][j] = Kij;
			}
		}

		// solve linear system
		K.solve(da, R);

		for (int i = 0; i < m; ++i) res[i] += da[i];

		niter++;
	} while (niter < MAX_ITER);

	delete f;

	return (niter < MAX_ITER);
}
