#include "stdafx.h"
#include "FEPowellOptimizeMethod.h"
#include "FECore/log.h"

//-----------------------------------------------------------------------------
// forward declarations
void linmin(double* p, double* xi, int n, double* fret, double (*fnc)(double[]));
void powell(double* p, double* xi, int n, double ftol, int* iter, double* fret, double (*fnc)(double[]));
double brent(double ax, double bx, double cx, double (*f)(double), double tol, double* xmin);
void fecb(FEModel* pfem, void* pd);
void mnbrak(double* ax, double* bx, double* cx, double* fa, double* fb, double* fc, double (*fnc)(double));
double golden(double ax, double bx, double cx, double (*f)(double), double tol, double* xmin);

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

#define FMAX(a,b) ((a)>(b) ? (a) : (b))

FEPowellOptimizeMethod* FEPowellOptimizeMethod::m_pThis = 0;

//-----------------------------------------------------------------------------
// FEPowellOptimizeMethod
//-----------------------------------------------------------------------------

bool FEPowellOptimizeMethod::Solve(FEOptimizeData *pOpt)
{
	int i, j;

	m_pOpt = pOpt;
	FEOptimizeData& opt = *pOpt;

	int nvar = opt.Variables();

	// set the initial guess
	double* p = new double[nvar];
	for (i=0; i<nvar; ++i) p[i] = opt.Variable(i).m_val;

	// set the initial search directions
	double* xi = new double[nvar*nvar];
	for (i=0; i<nvar; ++i)
	{
		for (j=0; j<nvar; ++j) xi[i*nvar + j] = 0;
		xi[i*nvar + i] = 0.05*(1.0 + p[i]);
	}

	opt.m_niter = 0;

	// set the FEM callback function
	FEModel& fem = opt.GetFEM();
	fem.AddCallback(fecb, CB_MAJOR_ITERS, &opt);

	// don't plot anything
	fem.GetCurrentStep()->SetPlotLevel(FE_PLOT_NEVER);

	// don't forget to set this
	m_pThis = this;

	// call the powell routine
	int niter = 0;
	double fret = 0;
	powell(p, xi, nvar, 0.001, &niter, &fret, objfun);

	felog.SetMode(Logfile::FILE_AND_SCREEN);

	felog.printf("\nP A R A M E T E R   O P T I M I Z A T I O N   R E S U L T S\n\n");

	felog.printf("\tMajor iterations ....................... : %d\n\n", niter);
	felog.printf("\tMinor iterations ....................... : %d\n\n", opt.m_niter);

	felog.printf("\tVariables:\n\n");
	for (i=0; i<nvar; ++i)
	{
		OPT_VARIABLE& var = opt.Variable(i);
		felog.printf("\t\t%-15s : %.16lg\n", var.m_szname, p[i]);
	}

	// evaluate reaction forces at correct times
	OPT_OBJECTIVE& obj = opt.GetObjective();
	FELoadCurve& rlc = opt.ReactionLoad();
	FELoadCurve& olc = opt.GetLoadCurve(obj.m_nlc);

	felog.printf("\n\tFunction values:\n\n");
	felog.printf("               CURRENT        REQUIRED      DIFFERENCE\n");
	for (i=0; i<olc.Points(); ++i)
	{
		LOADPOINT& p = olc.LoadPoint(i);
		double f = rlc.Value(p.time);
		felog.printf("%5d: %15.10lg %15.10lg %15lg\n", i+1, f, p.value, fabs(f - p.value));
	}

	felog.printf("\n\tFinal objective value: %15lg\n\n", fret);

	// done
	delete [] p;
	delete [] xi;

	return true;
}

//-------------------------------------------------------------------
void fecb(FEModel* pmdl, void* pd)
{
	// get the optimizaton data
	FEOptimizeData& opt = *((FEOptimizeData*) pd);

	// get the FEM data
	FEModel& fem = opt.GetFEM();

	// get the current time value
	double time = fem.m_ftime;

	// evaluate the current reaction force value
	OPT_OBJECTIVE& obj = opt.GetObjective();
	double value = *(obj.m_pd);

	// add the data pair to the loadcurve
	FELoadCurve& lc = opt.ReactionLoad();
	lc.Add(time, value);
}

//------------------------------------------------------------------

double FEPowellOptimizeMethod::ObjFun(double *p)
{
	// get the optimization data
	FEOptimizeData& opt = *m_pOpt;

	// increase iterator counter
	opt.m_niter++;

	// get the FEM data
	FEModel& fem = opt.GetFEM();

	// reset reaction force data
	FELoadCurve& lc = opt.ReactionLoad();
	lc.Clear();

	felog.printf("\n----- Iteration: %d -----\n", opt.m_niter);

	// set the material parameters
	int nvar = opt.Variables();
	for (int i=0; i<nvar; ++i)
	{
		OPT_VARIABLE& var = opt.Variable(i);
		*(var.m_pd) = p[i];
		felog.printf("  %-15s: %.16lg\n", var.m_szname, p[i]);
	}

	// reset the FEM data
	fem.Reset();

	// suppress output
	felog.SetMode(Logfile::NEVER);

	double fobj = 0;

	if (fem.Solve() == false)
	{
		felog.printf("\n\n\nAAAAAAAAARRRRRRRRRGGGGGGGGHHHHHHHHHHH !!!!!!!!!!!!!\n\n\n\n");
		return 0;
	}
	else
	{
		felog.SetMode(Logfile::FILE_AND_SCREEN);

		// evaluate reaction forces at correct times
		OPT_OBJECTIVE& obj = opt.GetObjective();
		FELoadCurve& rlc = opt.ReactionLoad();
		FELoadCurve& olc = opt.GetLoadCurve(obj.m_nlc);

		felog.printf("\n\tFunction values:\n\n");
		felog.printf("               CURRENT        REQUIRED      DIFFERENCE\n");
		for (int i=0; i<olc.Points(); ++i)
		{
			LOADPOINT& p = olc.LoadPoint(i);

			double f = rlc.Value(p.time);
			fobj += (f - p.value)*(f - p.value);
			felog.printf("%5d: %15.10lg %15.10lg %15lg\n", i+1, f, p.value, fabs(f - p.value));
		}
	}

	felog.printf("\n objective function: %.16lg\n", fobj);

	return fobj;
}

//=============================================================================
// powell's method
// from Numerical Recipes in C, section 10.5, page 417-419
// modified for this application
//
void powell(double* p, double* xi, int n, double ftol, int* iter, double* fret, double (*fnc)(double[]))
{
	int i,j,ibig;
	double fp,fptt,del,t;

	const int ITMAX = 200;
	const double TINY = 1.0e-25;

	double* pt = new double[n];
	double* ptt = new double[n];	
	double* xit = new double[n];

	*fret = (*fnc)(p);
	for (i=0; i<n; i++) pt[i] = p[i];				// save the initial point
	for (*iter=1;;++(*iter))
	{
		fp = *fret;
		ibig = 0;
		del= 0.0;		// will be biggest function decrease

		for (i=0;i<n;i++)	// in each iteration loop over all directions in the set
		{
			for (j=0;j<n;j++) xit[j] = xi[i*n+j];	// copy the direction
			fptt = (*fret);
			linmin(p,xit,n,fret,fnc);		// minimize along it
			if (fptt-(*fret) > del)			// and record it as the largets decrease so far 	
			{
				del = fptt - (*fret);
				ibig = i;
			}
		}

		// termination criterion
		if (2.0*(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret)) + TINY)
		{
		  delete [] pt;
		  delete [] ptt;
		  delete [] xit;
		  return;
		}

		// check we are not exceeding max number of iterations
		if (*iter == ITMAX)
		{
		  printf("FATAL ERROR : Max iterations reached\n");
		  exit(0);
		}

		for (j=0;j<n;j++)
		{
			ptt[j] = 2.0*p[j] - pt[j];
			xit[j] = p[j] - pt[j];
			pt[j] = p[j];
		}

		fptt = (*fnc)(ptt);
		if (fptt < fp)
		{
			t = 2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
			if (t < 0.0)
			{
				linmin(p, xit, n, fret, fnc);
				for (j=0;j<n;j++)
				{
					xi[ibig *n + j] = xi [(n-1)*n + j];
					xi[(n-1)*n + j] = xit[j];
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
double (*nrfunc)(double[]);

double f1dim(double x)
{
	int j;
	double f;
	for (j=0;j<ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f = (*nrfunc)(xt);
	return f;
}

void linmin(double* p, double* xi, int n, double* fret, double (*fnc)(double[]))
{
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;
	const double TOL = 2.0e-4;

	ncom = n;
	pcom = new double[n];
	xicom = new double[n];
	xt = new double[n];

	nrfunc = fnc;
	for (j=0;j<n;j++)
	{
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax = 0.0;
	xx = 1.0e-4;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	*fret = brent(ax,xx,bx,f1dim,TOL,&xmin);
	for (j=0;j<n;j++)
	{
		xi[j] *= xmin;
		p[j] += xi[j];
	}

	delete [] xt;
	delete [] pcom;
	delete [] xicom;
}


//--------------------------------------------------------------------------------
// routine to find 1D minimum
// from Numerical Recipes in C, section 10.3, page 404-405
// modified for this application
//
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b)>=0.0 ? (a) : (-(a)))

double brent(double ax, double bx, double cx, double (*f)(double), double tol, double* xmin)
{
	const int ITMAX = 100;
	const double CGOLD = 0.3819660;
	const double ZEPS = 1.0e-10;
	int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;

	a = (ax < cx ? ax : cx);
	b = (ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x);
	for (iter=1;iter<=ITMAX;iter++)
	{
		xm = 0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2 - 0.5*(b-a)))
		{
			*xmin = x;
			return fx;
		}
		if (fabs(e) > tol1)
		{
			r = (x-w)*(fx-fv);
			q = (x-v)*(fx-fw);
			p = (x-v)*q - (x-w)*r;
			q = 2.0*(q-r);
			if (q > 0.0) p = -p;
			q = fabs(q);
			etemp = e;
			e=d;
			if (fabs(p)>=fabs(0.5*q*etemp)||p<=q*(a-x)||p>=q*(b-x))
				d = CGOLD*(e=(x >= xm ? a-x : b-x));
			else
			{
				d = p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2) d=SIGN(tol1,xm-x);
			}
		}
		else d = CGOLD*(e=(x >= xm ? a-x : b-x));
		u = (fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu = (*f)(u);

		if (fu <= fx)
		{
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u);
			SHFT(fv,fw,fx,fu);
		}
		else
		{
			if (u<x) a=u; else b=u;
			if (fu <= fw || w==x)
			{
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			}
			else if (fu <= fv || v==x || v==w)
			{
				v=u;
				fv=fu;
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

void mnbrak(double* pa, double* pb, double* pc, double* pfa, double* pfb, double* pfc, double (*func)(double))
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

	c = b+GOLD*(b - a);
	fc = func(c);
	while (fb > fc)
	{
		r = (b - a)*(fb - fc);
		q = (b - c)*(fb - fa);
		u = b - ((b - c)*q - (b - a)*r) / (2.0*SIGN2(FMAX(fabs(q-r), TINY), q-r));

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
		else if ((u-ulim)*(ulim - c) >= 0)
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

double golden(double ax, double bx, double cx, double (*f)(double), double tol, double* xmin)
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

