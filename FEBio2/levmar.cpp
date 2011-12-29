#include "stdafx.h"
#include "FEOptimizer.h"
#include "log.h"
#include "console.h"

FELMOptimizeMethod* FELMOptimizeMethod::m_pThis = 0;

// forward declarations
void mrqmin(vector<double>& x, 
			vector<double>& y, 
			vector<double>& sig, 
			vector<double>& a, 
			matrix& covar, 
			matrix& alpha, 
			double& chisq,
			void funcs(vector<double>& , vector<double>&, vector<double>&, matrix&),
			double& alamda);

void mrqcof(vector<double>& x, 
			vector<double>& y, 
			vector<double>& sig, 
			vector<double>& a, 
			matrix& alpha,
			vector<double>& beta,
			double& chisq,
			void funcs(vector<double>& , vector<double>&, vector<double>&, matrix&));

extern void fecb(FEModel* pfem, void* pd);

//-----------------------------------------------------------------------------
FELMOptimizeMethod::FELMOptimizeMethod()
{
	m_objtol = 0.001;
	m_fdiff  = 0.001;
}

//-----------------------------------------------------------------------------
bool FELMOptimizeMethod::Solve(FEOptimizeData *pOpt)
{
	m_pOpt = pOpt;
	FEOptimizeData& opt = *pOpt;

	int i;

	// set the variables
	int ma = opt.Variables();
	vector<double> a(ma);
	for (i=0; i<ma; ++i)
	{
		OPT_VARIABLE& var = opt.Variable(i);
		a[i] = var.m_val;
	}

	// set the data
	OPT_OBJECTIVE& obj = opt.GetObjective();
	FELoadCurve& lc = opt.GetLoadCurve(obj.m_nlc);
	int ndata = lc.Points();
	vector<double> x(ndata), y(ndata);
	for (i=0; i<ndata; ++i) 
	{
		x[i] = lc.LoadPoint(i).time;
		y[i] = lc.LoadPoint(i).value;
	}
	m_y0 = y;

	// set the sigma's
	// for now we set them all to 1
	vector<double> sig(ndata);
	for (i=0; i<ndata; ++i) sig[i] = 1;

	// allocate matrices
	matrix covar(ma, ma), alpha(ma, ma);

	// set the FEM callback function
	FEM& fem = opt.GetFEM();
	fem.AddCallback(fecb, &opt);

	// don't plot anything
	fem.m_pStep->m_nplot = FE_PLOT_NEVER;

	// set the this pointer
	m_pThis = this;

	opt.m_niter = 0;

	// return value
	double fret = 0.0;

	clog.SetMode(Logfile::FILE_AND_SCREEN);

	int niter = 1;

	try
	{
		// do the first call with lamda to intialize the minimization
		double alamda = -1.0;
		clog.printf("\n----- Major Iteration: %d -----\n", 0);
		mrqmin(x, y, sig, a, covar, alpha, fret, objfun, alamda);

		// repeat until converged
		double fprev = fret, lam1 = alamda;
		bool bconv = false;
		int NMAX = 100;
		do
		{
			clog.printf("\n----- Major Iteration: %d -----\n", niter);
			mrqmin(x, y, sig, a, covar, alpha, fret, objfun, alamda);

			if (alamda < lam1)
			{
				if (niter != 1)
				{
					double df = (fprev - fret)/(fprev + fret + 1);
					if ( df < m_objtol) bconv = true;
					clog.printf("objective value: %lg (diff = %lg)\n\n", fret, df);
				}
			}
			else clog.printf("\n objective value: %lg\n\n", fret);

			fprev = fret;
			lam1 = alamda;

			++niter;
		}
		while ((bconv == false) && (niter < NMAX));

		// do final call with lamda = 0
		alamda = 0.0;
		mrqmin(x, y, sig, a, covar, alpha, fret, objfun, alamda);

	}
	catch (FEErrorTermination)
	{
		clog.printbox("F A T A L   E R R O R", "FEBio error terminated. Parameter optimization cannot continue.");
		return false;
	}

	clog.SetMode(Logfile::FILE_AND_SCREEN);

	clog.printf("\nP A R A M E T E R   O P T I M I Z A T I O N   R E S U L T S\n\n");

	clog.printf("\tMajor iterations ....................... : %d\n\n", niter);
	clog.printf("\tMinor iterations ....................... : %d\n\n", opt.m_niter);

	clog.printf("\tVariables:\n\n");
	for (i=0; i<ma; ++i)
	{
		OPT_VARIABLE& var = opt.Variable(i);
		clog.printf("\t\t%-15s : %.16lg\n", var.m_szname, a[i]);
	}

	// print reaction forces
	clog.printf("\n\tFunction values:\n\n");
	clog.printf("               CURRENT        REQUIRED      DIFFERENCE\n");
	for (i=0; i<ndata; ++i)
	{
		clog.printf("%5d: %15.10lg %15.10lg %15lg\n", i+1, m_yopt[i], y[i], fabs(m_yopt[i] - y[i]));
	}

	clog.printf("\n\tFinal objective value: %15lg\n\n", fret);

	return true;
}

//-----------------------------------------------------------------------------
void FELMOptimizeMethod::ObjFun(vector<double>& x, vector<double>& a, vector<double>& y, matrix& dyda)
{
	// get the optimization data
	FEOptimizeData& opt = *m_pOpt;

	// evaluate at a
	bool bret = FESolve(x, a, y);
	if (bret == false) throw FEErrorTermination();
	
	m_yopt = y;

	// now calculate the derivatives using forward differences
	int ndata = x.size();
	vector<double> a1(a);
	vector<double> y1(ndata);
	int ma = a.size();
	for (int i=0; i<ma; ++i)
	{
		double b = opt.Variable(i).m_sf;

		a1[i] = a1[i] + m_fdiff*(fabs(b) + fabs(a[i]));
		assert(a1[i] != a[i]);

		FESolve(x, a1, y1);
		for (int j=0; j<ndata; ++j) dyda[j][i] = (y1[j] - y[j])/(a1[i] - a[i]);
		a1[i] = a[i];
	}
}

//-----------------------------------------------------------------------------
bool FELMOptimizeMethod::FESolve(vector<double> &x, vector<double> &a, vector<double> &y)
{
	// get the optimization data
	FEOptimizeData& opt = *m_pOpt;

	// increase iterator counter
	opt.m_niter++;

	// get the FEM data
	FEM& fem = opt.GetFEM();

	// reset reaction force data
	FELoadCurve& lc = opt.ReactionLoad();
	lc.Clear();

	// set the material parameters
	int nvar = opt.Variables();
	for (int i=0; i<nvar; ++i)
	{
		OPT_VARIABLE& var = opt.Variable(i);
		*(var.m_pd) = a[i];
	}

	// reset the FEM data
	fem.Reset();

	clog.SetMode(Logfile::FILE_AND_SCREEN);
	clog.printf("\n----- Iteration: %d -----\n", opt.m_niter);
	for (int i=0; i<nvar; ++i) 
	{
		OPT_VARIABLE& var = opt.Variable(i);
		clog.printf("%-15s = %lg\n", var.m_szname, a[i]);
	}

	// solve the FE problem
	clog.SetMode(Logfile::NEVER);
	Console* pwnd = Console::GetHandle();
	pwnd->Deactivate();

	bool bret = fem.Solve();

	clog.SetMode(Logfile::FILE_AND_SCREEN);
	if (bret)
	{
		FELoadCurve& rlc = opt.ReactionLoad();
		int ndata = x.size();
		clog.printf("               CURRENT        REQUIRED      DIFFERENCE\n");
		for (int i=0; i<ndata; ++i) 
		{
			y[i] = rlc.Value(x[i]);
			clog.printf("%5d: %15.10lg %15.10lg %15lg\n", i+1, y[i], m_y0[i], fabs(y[i] - m_y0[i]));
		}
	}

	return bret;
}

//-----------------------------------------------------------------------------
void mrqmin(vector<double>& x, 
			vector<double>& y, 
			vector<double>& sig, 
			vector<double>& a, 
			matrix& covar, 
			matrix& alpha, 
			double& chisq,
			void funcs(vector<double>& , vector<double>&, vector<double>&, matrix&),
			double& alamda)
{
	static double ochisq;
	int j, k, l;

	int ma = a.size();
	static vector<double> oneda(ma);
	static vector<double> atry(ma), beta(ma), da(ma);
	if (alamda < 0)
	{
		alamda = 0.001;
		mrqcof(x, y, sig, a, alpha, beta, chisq, funcs);
		ochisq = chisq;
		for (j=0; j<ma; j++) atry[j] = a[j];
	}
	matrix temp(ma, ma);
	for (j=0; j<ma; j++)
	{
		for (k=0; k<ma; k++) covar[j][k] = alpha[j][k];
		covar[j][j] = alpha[j][j]*(1.0 + alamda);
		for (k=0; k<ma; k++) temp[j][k] = covar[j][k];
		oneda[j] = beta[j];
	}
	matrix tempi = temp.inverse();
	oneda = tempi*oneda;
	for (j=0; j<ma; j++)
	{
		for (k=0; k<ma; k++) covar[j][k] = tempi[j][k];
		da[j] = oneda[j];
	}

	if (alamda == 0.0) return;

	for (j=0, l=0; l<ma; l++) atry[l] = a[l] + da[j++];
	mrqcof(x, y, sig, atry, covar, da, chisq, funcs);
	if (chisq < ochisq)
	{
		alamda *= 0.1;
		ochisq = chisq;
		for (j=0; j<ma; j++)
		{
			for (k=0; k<ma; k++) alpha[j][k] = covar[j][k];
			beta[j] = da[j];
		}
		for (l=0; l<ma; l++) a[l] = atry[l];
	}
	else
	{
		alamda *= 10.0;
		chisq = ochisq;
	}
}

//-----------------------------------------------------------------------------
void mrqcof(vector<double>& x, 
			vector<double>& y, 
			vector<double>& sig, 
			vector<double>& a, 
			matrix& alpha,
			vector<double>& beta,
			double& chisq,
			void funcs(vector<double>& , vector<double>&, vector<double>&, matrix&))
{
	int i, j, k, l, m;
	double wt, sig2i, dy;

	int ndata = x.size();
	int ma = a.size();
	for (j=0; j<ma; j++)
	{
		for (k=0; k<=j; k++) alpha[j][k] = 0.0;
		beta[j] = 0.0;
	}

	vector<double> ymod(ndata);
	matrix dyda(ndata, ma);
	funcs(x, a, ymod, dyda);
	
	chisq = 0.0;
	for (i=0; i<ndata; i++)
	{
		sig2i = 1.0 / (sig[i]*sig[i]);
		dy = y[i] - ymod[i];
		for (j=0, l=0; l<ma; l++)
		{
			wt = dyda[i][l]*sig2i;
			for (k=0, m=0; m<l+1; m++) alpha[j][k++] += wt*dyda[i][m];
			beta[j++] += dy*wt;
		}
		chisq += dy*dy*sig2i;
	}
	for (j=1; j<ma; j++)
		for (k=0; k<j; k++) alpha[k][j] = alpha[j][k];
}
