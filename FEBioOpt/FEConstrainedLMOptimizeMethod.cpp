#include "stdafx.h"
#include "FEConstrainedLMOptimizeMethod.h"
#include "FECore/log.h"

#ifdef HAVE_LEVMAR
#include "levmar.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEConstrainedLMOptimizeMethod, FEOptimizeMethod)
	ADD_PARAMETER(m_objtol, FE_PARAM_DOUBLE, "obj_tol"     );
	ADD_PARAMETER(m_tau   , FE_PARAM_DOUBLE, "tau"         );
	ADD_PARAMETER(m_fdiff , FE_PARAM_DOUBLE, "f_diff_scale");
	ADD_PARAMETER(m_nmax  , FE_PARAM_INT   , "max_iter"    );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEConstrainedLMOptimizeMethod* FEConstrainedLMOptimizeMethod::m_pThis = 0;

//-----------------------------------------------------------------------------
extern bool fecb(FEModel* pfem, unsigned int, void* pd);

//-----------------------------------------------------------------------------
void clevmar_cb(double *p, double *hx, int m, int n, void *adata)
{
	FEConstrainedLMOptimizeMethod* pLM = (FEConstrainedLMOptimizeMethod*) adata;

	// get the optimization data
	FEOptimizeData& opt = *pLM->GetOptimizeData();

	// set the variables
	vector<double> a(m);
	for (int i=0; i<m; ++i) a[i] = p[i];

	// set the data
	OPT_OBJECTIVE& obj = opt.GetObjective();
	FELoadCurve& lc = opt.GetLoadCurve(obj.m_nlc);
	vector<double> x(n), y(n);
	for (int i=0; i<n; ++i) x[i] = lc.LoadPoint(i).time;
	
	// evaluate at a
	if (pLM->FESolve(x, a, y) == false) throw FEErrorTermination();

	// store the measurement vector
	for (int i=0; i<n; ++i) hx[i] = y[i];

	// store the last calculated values
	pLM->m_yopt = y;
}

//-----------------------------------------------------------------------------
FEConstrainedLMOptimizeMethod::FEConstrainedLMOptimizeMethod()
{
	m_tau = 1e-3;
	m_objtol = 0.001;
	m_fdiff  = 0.001;
	m_nmax   = 100;
    m_loglevel = Logfile::NEVER;
}

//-----------------------------------------------------------------------------
bool FEConstrainedLMOptimizeMethod::Solve(FEOptimizeData *pOpt)
{
	m_pOpt = pOpt;
	FEOptimizeData& opt = *pOpt;

	// set the variables
	int ma = opt.Variables();
	vector<double> a(ma);
	for (int i=0; i<ma; ++i)
	{
		OPT_VARIABLE& var = opt.Variable(i);
		a[i] = var.m_val;
	}

	// set the data
	OPT_OBJECTIVE& obj = opt.GetObjective();
	FELoadCurve& lc = opt.GetLoadCurve(obj.m_nlc);
	int ndata = lc.Points();
	vector<double> x(ndata), y(ndata);
	for (int i=0; i<ndata; ++i) 
	{
		x[i] = lc.LoadPoint(i).time;
		y[i] = lc.LoadPoint(i).value;
	}
	m_y0 = y;

	// set the sigma's
	// for now we set them all to 1
	vector<double> sig(ndata);
	for (int i=0; i<ndata; ++i) sig[i] = 1;

	// allocate matrices
	matrix covar(ma, ma), alpha(ma, ma);

	// set the FEM callback function
	FEModel& fem = opt.GetFEM();
	fem.AddCallback(fecb, CB_MAJOR_ITERS | CB_INIT, &opt);

	// don't plot anything
	fem.GetCurrentStep()->SetPlotLevel(FE_PLOT_NEVER);

	// set the this pointer
	m_pThis = this;

	opt.m_niter = 0;

	// return value
	double fret = 0.0;

	felog.SetMode(Logfile::FILE_AND_SCREEN);

	int niter = 1;

	try
	{
		double* p = new double[ma];
		for (int i=0; i<ma; ++i) p[i] = a[i];

		double* lb = new double[ma];
		for (int i=0; i<ma; ++i) lb[i] = opt.Variable(i).m_min;

		double* ub = new double[ma];
		for (int i=0; i<ma; ++i) ub[i] = opt.Variable(i).m_max;

		double* q = new double[ndata];
		for (int i=0; i<ndata; ++i) q[i] = y[i];

		const double tol = m_objtol;
		double opts[5] = {m_tau, tol, tol, tol, m_fdiff};

		int itmax = m_nmax;
		if (opt.Constraints() > 0)
		{
			int NC = opt.Constraints();
			double* A = new double[NC*ma];
			double* b = new double[NC];
			for (int i=0; i<NC; ++i)
			{
				OPT_LIN_CONSTRAINT& con = opt.Constraint(i);
				for (int j=0; j<ma; ++j) A[i*ma + j] = con.a[j];
				b[i] = con.b;
			}

			int ret = dlevmar_blec_dif(clevmar_cb, p, q, ma, ndata, lb, ub, A, b, NC, 0, itmax, opts, 0, 0, 0, (void*) this);

			delete [] b;
			delete [] A;
		}
		else
		{
			int ret = dlevmar_bc_dif(clevmar_cb, p, q, ma, ndata, lb, ub, 0, itmax, opts, 0, 0, 0, (void*) this);
		}

		for (int i=0; i<ma; ++i) a[i] = p[i];

		// store the optimal values
		fret = 0.0;
		for (int i=0; i<ndata; ++i)
		{
			double dy = m_yopt[i] - m_y0[i];
			fret += dy*dy;
		}

		delete [] q;
		delete [] ub;
		delete [] lb;
		delete [] p;
	}
	catch (FEErrorTermination)
	{
		felog.printbox("F A T A L   E R R O R", "FEBio error terminated. Parameter optimization cannot continue.");
		return false;
	}

	felog.SetMode(Logfile::FILE_AND_SCREEN);

	felog.printf("\nP A R A M E T E R   O P T I M I Z A T I O N   R E S U L T S\n\n");

	// print reaction forces
	felog.printf("\n\tFunction values:\n\n");
	felog.printf("               CURRENT        REQUIRED      DIFFERENCE\n");
	for (int i=0; i<ndata; ++i)
	{
		felog.printf("%5d: %15.10lg %15.10lg %15lg\n", i+1, m_yopt[i], y[i], fabs(m_yopt[i] - y[i]));
	}
    
	felog.printf("\n\tFinal objective value: %15lg\n\n", fret);
    
	felog.printf("\tMajor iterations ....................... : %d\n\n", niter);
	felog.printf("\tMinor iterations ....................... : %d\n\n", opt.m_niter);

	felog.printf("\tVariables:\n\n");
	for (int i=0; i<ma; ++i)
	{
		OPT_VARIABLE& var = opt.Variable(i);
		felog.printf("\t\t%-15s : %.16lg\n", var.m_szname, a[i]);
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEConstrainedLMOptimizeMethod::ObjFun(vector<double>& x, vector<double>& a, vector<double>& y, matrix& dyda)
{
	// get the optimization data
	FEOptimizeData& opt = *m_pOpt;

	// poor man's box constraints
	int ma = a.size();
	vector<int> dir(ma,1);	// forward difference by default
	for (int i=0; i<opt.Variables(); ++i)
	{
		OPT_VARIABLE& var = opt.Variable(i);
		if (a[i] < var.m_min) {
			a[i] = var.m_min;
		} else if (a[i] >= var.m_max) {
			a[i] = var.m_max;
			dir[i] = -1;	// use backward difference
		}
	}
	
	// evaluate at a
	if (FESolve(x, a, y) == false) throw FEErrorTermination();
	
	m_yopt = y;

	// now calculate the derivatives using forward differences
	int ndata = x.size();
	vector<double> a1(a);
	vector<double> y1(ndata);
	for (int i=0; i<ma; ++i)
	{
		double b = opt.Variable(i).m_sf;

		a1[i] = a1[i] + dir[i]*m_fdiff*(b + fabs(a[i]));

		if (FESolve(x, a1, y1) == false) throw FEErrorTermination();
		for (int j=0; j<ndata; ++j) dyda[j][i] = (y1[j] - y[j])/(a1[i] - a[i]);
		a1[i] = a[i];
	}
}

//-----------------------------------------------------------------------------
bool FEConstrainedLMOptimizeMethod::FESolve(vector<double> &x, vector<double> &a, vector<double> &y)
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

	// set the material parameters
	int nvar = opt.Variables();
	for (int i=0; i<nvar; ++i)
	{
		OPT_VARIABLE& var = opt.Variable(i);
		*(var.m_pd) = a[i];
	}

	// reset the FEM data
	fem.Reset();

	felog.SetMode(Logfile::FILE_AND_SCREEN);
	felog.printf("\n----- Iteration: %d -----\n", opt.m_niter);
	for (int i=0; i<nvar; ++i) 
	{
		OPT_VARIABLE& var = opt.Variable(i);
		felog.printf("%-15s = %lg\n", var.m_szname, a[i]);
	}

	// solve the FE problem
	felog.SetMode((Logfile::MODE)m_loglevel);

	bool bret = m_pOpt->RunTask();

	felog.SetMode(Logfile::FILE_AND_SCREEN);
	if (bret)
	{
		double chisq = 0.0;
		FELoadCurve& rlc = opt.ReactionLoad();
		int ndata = x.size();
		if (m_print_level == PRINT_VERBOSE) felog.printf("               CURRENT        REQUIRED      DIFFERENCE\n");
		for (int i=0; i<ndata; ++i) 
		{
			y[i] = rlc.Value(x[i]);
			double dy = (m_y0[i] - y[i]);
			chisq += dy*dy;
			if (m_print_level == PRINT_VERBOSE) felog.printf("%5d: %15.10lg %15.10lg %15lg\n", i+1, y[i], m_y0[i], fabs(y[i] - m_y0[i]));
		}
		felog.printf("objective value: %lg\n", chisq);
	}


	return bret;
}
#endif
