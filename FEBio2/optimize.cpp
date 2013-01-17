// optimize.cpp : implementation of a parameter optimization routine
// 
//////////////////////////////////////////////////////////////////////////
#include "stdafx.h"
#include "FEBioLib/FEOptimizer.h"
#include "FECore/log.h"

#ifdef NAGLIB

#include <nag.h>
#include <nage04.h>
#include <time.h>
#include "FEBioXML/XMLReader.h"
#include "FEBioLib/Timer.h"
#include "console.h"

//-----------------------------------------------------------------------------
//! Parameter optimization data

//! This class stores all the variables that are needed to do a parameter
//! optimization.

struct FE_OPTIMIZE
{
	enum {MAX_VARS   =   8};	//!< maximum number of variables to optimize
	enum {MAX_FUNCS  = 128};	//!< maximum number of subfunctions

	int		nprint;				//!< print flag

	int		nvar;					//!< number of variables to optimize
	double*	pvar[MAX_VARS];			//!< pointers to variables to optimize
	char	szvar[MAX_VARS][80];	//!< name of variables

	double* pobj;				//!< variable that generates objective function data

	double	c0[MAX_VARS];		//!< initial values of variables
	double	cf[MAX_VARS];		//!< final values of variables
	double	bl[MAX_VARS];		//!< lower bounds of variables
	double	bu[MAX_VARS];		//!< upper bounds of variables

	int		nfunc;				//!< number of subfunctions
	double	x[MAX_FUNCS];		//!< "time" values for the subfunctions
	double	y[MAX_FUNCS];		//!< y-values for the subfunctions

	double	diff;				//!< finite difference interval
	double	prec;				//!< expected precision of objective function
	double	lstol;				//!< line search tolerance

	FEModel*	pfem;				//!< pointer to FEM object

	FELoadCurve	lc;				//!< loadcurve containing the reaction forces

	Nag_E04_Opt	op;				//!< optional parameters

	//! constructor
	FE_OPTIMIZE()
	{
		nag_opt_init(&op);
		op.obj_deriv      = Nag_FALSE;
		op.con_deriv      = Nag_FALSE;
		op.verify_grad    = Nag_NoCheck;
	}
};

//-----------------------------------------------------------------------------
// forward declarations

bool ReadScript(FEModel& fem, XMLReader& xml, XMLTag& tag, FE_OPTIMIZE& opt);
void objfun(Integer m, Integer n, double x[], double f[],
			double fjac[], Integer tdfjac, Nag_Comm* comm);
void fecb(FEModel*, void*);

//-----------------------------------------------------------------------------
//! Perform a parameter optimization.
//! This function calls the NAG routine to do the actual parameter optimization.

bool optimize(FEModel& fem, const char* szfile)
{
	// try to open the file
	XMLReader xml;
	if (xml.Open(szfile) == false) return false;

	// find the root element
	XMLTag tag;
	if (xml.FindTag("Paramopt", tag) == false) return false;

	// check the version number
	if (strcmp(tag.m_szatv[0], "2.0") != 0) 
	{
		fprintf(stderr, "Invalid script version.");
		return false;
	}

	// Read the input file
	FE_OPTIMIZE opt;
	try
	{
		if (!ReadScript(fem, xml, tag, opt))
		{
			fprintf(stderr, "FATAL ERROR: Failed reading script file.\n\n");
			return false;
		}
	}
	catch (InvalidVariableName e)
	{
		fprintf(stderr, "FATAL ERROR: the variable %s is not recognized.\n\n", e.szname);
		return false;
	}
	catch (NothingToOptimize)
	{
		fprintf(stderr, "FATAL ERROR: there is nothing to optimize.\n\n");
		return false;
	}
	catch (XMLReader::InvalidTag e)
	{
		fprintf(stderr, "FATAL ERROR: unrecognized tag \"%s\" (line %d)\n", e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::InvalidValue e)
	{
		fprintf(stderr, "FATAL ERROR: The element %s has an invalid value.\n\n", e.tag.Name(), e.tag.szvalue());
		return false;
	}
	catch (...)
	{
		fprintf(stderr, "FATAL ERROR: an exception occured in the optimize routine.\n\n");
		return false;
	}

	// data variables for nag routine
	double objf;
	double fjac[FE_OPTIMIZE::MAX_FUNCS*FE_OPTIMIZE::MAX_VARS];
	double yf[FE_OPTIMIZE::MAX_FUNCS];

	opt.pfem = &fem;

	// Nag error handling
	NagError fail;
	fail.print = Nag_TRUE;
	fail.handler = 0;

	// NAG communication structure
	Nag_Comm comm;
	comm.p = (void*) &opt;

	// keep track of time
	Timer timer;

	// call the NAG routine
	timer.start();
	{
		e04unc(
			opt.nfunc,				// nr of subfunctions
			opt.nvar,				// nr of parameters
			0, 0, 0, 0,				// (zeroes, not used)
			opt.bl,					// lower bounds of parameters
			opt.bu,					// upper bounds of parameters
			opt.y,					// values of y vector
			objfun,					// pointer to the objective function
			NULLFN,					// pointer to the nonlinear constraint function (not used)
			opt.cf,					// initial (input) / final (output) optimized values
			&objf,					// final objective function value
			yf,						// final subfunction values
			fjac,					// final subfunction derivatives
			FE_OPTIMIZE::MAX_VARS,	// second dimension of fjac
			&opt.op,				// optional parameters
			&comm,					// communication information
			&fail);					// fail return code
	}
	timer.stop();

	if ((opt.nprint) || (fail.code == NE_NOERROR))
	{
		int i;
		printf("\nO P T I M I Z A T I O N   S U M M A R Y\n\n");
		printf(" var        initial           final           lower           upper\n");
		for (i=0; i<opt.nvar; ++i)
			printf(" V%d %15lg %15lg %15lg %15lg\n", i+1, opt.c0[i], opt.cf[i], opt.bl[i], opt.bu[i]);

		printf("\n           Time         y-value           final           diff\n");
		for (i=0; i<opt.nfunc; ++i)
		{
			double diff = opt.y[i] - yf[i];
			printf("%15lg %15lg %15lg%15lg\n", opt.x[i], opt.y[i], yf[i], diff);
		}

		char sztime[64];
		timer.time_str(sztime);
		printf("\n time elapsed: %s\n", sztime);
	}

	if (fail.code == NE_NOERROR)
	{
		printf("\nN O R M A L   T E R M I N A T I O N\n\n");
	}
	else
	{
		printf("\nE R R O R   T E R M I N A T I O N\n\n");
	}

  return true;
}

//-----------------------------------------------------------------------------

void objfun(Integer m, 
			Integer n, 
			double x[], 
			double f[],
			double fjac[], 
			Integer tdfjac, 
			Nag_Comm* comm)
{
	int i;

	// get the optimization data
	FE_OPTIMIZE& opt = *((FE_OPTIMIZE*) comm->p);

	// get the FE Model data
	FEModel& fem = *(opt.pfem);

	// reset reaction force data
	opt.lc.Clear();

	// set the callback function
	fem.AddCallback(fecb, &opt);

	// set the material parameters
	for (i=0; i<n; ++i) *(opt.pvar[i]) = x[i];

	// reset the FE data
	fem.Reset();

	// suppress output to the screen
	clog.SetMode(Logfile::FILE_ONLY);
	Console* pwnd = Console::GetHandle();
	pwnd->Deactivate();

	// solve the problem
	if (opt.nprint > 0)
	{
		static int niter = 1;

		fprintf(stderr, "\nSolving FE problem (iteration %d):\n", niter++);
		for (i=0; i<n; ++i)
			fprintf(stderr, "\t%s: %lg\n", opt.szvar[i], x[i]);
	}
	if (fem.Solve() == false)
	{
		comm->flag = -1;
	}
	else
	{
		// evaluate reaction forces at correct times
		double fobj = 0;
		for (i=0; i<m; ++i)
		{
			f[i] = opt.lc.Value(opt.x[i]);
			fobj += (f[i] - opt.y[i])*(f[i] - opt.y[i]);
		}

		fobj *= 0.5;

		if (opt.nprint > 0)
		{
			fprintf(stderr, "\tObjective function: %lg\n", fobj);
		}
	}
}

//-----------------------------------------------------------------------------

void fecb(FEModel* pfem, void* pd)
{
	// get the optimizaton data
	FE_OPTIMIZE& opt = *((FE_OPTIMIZE*) pd);

	// get the FE Model data
	FEModel& fem = *(opt.pfem);

	// get the current time value
	double time = fem.m_ftime;

	// evaluate the current reaction force value
	double value = *opt.pobj;

	// add the data pair to the loadcurve
	opt.lc.Add(time, value);
}

//-----------------------------------------------------------------------------

bool ReadScript(FEModel& fem, XMLReader& xml, XMLTag& tag, FE_OPTIMIZE& opt)
{
	// initialize optimization data
	opt.nvar  = 0;
	opt.nfunc = 0;
	opt.diff = 0.001;
	opt.prec = 0.0001;
	opt.lstol = 0.9;
	opt.pobj = 0;
	opt.nprint = 0;

	double* pd;
	const char* sz;
	double d[3];

	// read the optimization data
	++tag;
	do
	{
		if (tag == "model")
		{
			// get the input file
			char szfile[256];
			tag.value(szfile);

			// read input data
			if (fem.Input(szfile) == false) return false;

			// initialize data 
			if (fem.Init() == false) return false;
		}
		else if (tag == "variables")
		{
			// read the parameters
			++tag;
			do
			{
				if (tag == "var")
				{
					if (opt.nvar >= FE_OPTIMIZE::MAX_VARS) throw XMLReader::InvalidTag(tag);

					// get the variable name
					sz = tag.AttributeValue("name");
					if (sz == 0) throw InvalidVariableName("[Unknown]");

					strcpy(opt.szvar[ opt.nvar ], sz);

					// find the variable
					pd = fem.FindParameter(sz);
					if (pd == 0) throw InvalidVariableName(sz);
					opt.pvar[ opt.nvar ] = pd;

					// set initial values and bounds
					tag.value(d, 3);
					opt.c0[ opt.nvar ] = opt.cf[ opt.nvar ] = d[0];
					opt.bl[ opt.nvar ] = d[1];
					opt.bu[ opt.nvar ] = d[2];

					// increase the variable counter
					opt.nvar++;
				}
				else throw XMLReader::InvalidTag(tag);

				++tag;
			}
			while (!tag.isend());
		}
		else if (tag == "sub_functions")
		{
			// read the subfunction values
			++tag;
			do
			{
				if (tag == "func")
				{
					if (opt.nfunc >= FE_OPTIMIZE::MAX_FUNCS) throw XMLReader::InvalidTag(tag);

					// get the function ordinate and coordinate value
					tag.value(d, 2);

					opt.x[ opt.nfunc ] = d[0];
					opt.y[ opt.nfunc ] = d[1];

					// increase function counter
					opt.nfunc++;
				}
				else throw XMLReader::InvalidTag(tag);				

				++tag;
			}
			while (!tag.isend());
		}
		else if (tag == "objective")
		{
			// get the variable name
			char szval[256];
			tag.value(szval);

			// find the variable
			pd = fem.FindParameter(szval);
			if (pd == 0) throw InvalidVariableName(szval);
			opt.pobj = pd;
		}
		else if (tag == "print") tag.value(opt.nprint);
		else if (tag == "Nag_options")
		{
			// read the NAG options
			++tag;
			do
			{
				if (tag == "list")
				{
					if (strcmp(tag.szvalue(), "TRUE" ) == 0) opt.op.list = Nag_TRUE;
					else if (strcmp(tag.szvalue(), "FALSE") == 0) opt.op.list = Nag_FALSE;
					else throw XMLReader::InvalidValue(tag);
				}
				else if (tag == "print_level")
				{
					if      (strcmp(tag.szvalue(), "Nag_NoPrint"        ) == 0) opt.op.print_level = Nag_NoPrint;
					else if (strcmp(tag.szvalue(), "Nag_Soln"           ) == 0) opt.op.print_level = Nag_Soln;
					else if (strcmp(tag.szvalue(), "Nag_Iter"           ) == 0) opt.op.print_level = Nag_Iter;
					else if (strcmp(tag.szvalue(), "Nag_Iter_Long"      ) == 0) opt.op.print_level = Nag_Iter_Long;
					else if (strcmp(tag.szvalue(), "Nag_Soln_Iter"      ) == 0) opt.op.print_level = Nag_Soln_Iter;
					else if (strcmp(tag.szvalue(), "Nag_Soln_Iter_Long" ) == 0) opt.op.print_level = Nag_Soln_Iter_Long;
					else if (strcmp(tag.szvalue(), "Nag_Soln_Iter_Const") == 0) opt.op.print_level = Nag_Soln_Iter_Const;
					else if (strcmp(tag.szvalue(), "Nag_Soln_Iter_Full" ) == 0) opt.op.print_level = Nag_Soln_Iter_Full;
					else throw XMLReader::InvalidValue(tag);
				}
				else if (tag == "minor_print_level")
				{
					if      (strcmp(tag.szvalue(), "Nag_NoPrint"        ) == 0) opt.op.minor_print_level = Nag_NoPrint;
					else if (strcmp(tag.szvalue(), "Nag_Soln"           ) == 0) opt.op.minor_print_level = Nag_Soln;
					else if (strcmp(tag.szvalue(), "Nag_Iter"           ) == 0) opt.op.minor_print_level = Nag_Iter;
					else if (strcmp(tag.szvalue(), "Nag_Iter_Long"      ) == 0) opt.op.minor_print_level = Nag_Iter_Long;
					else if (strcmp(tag.szvalue(), "Nag_Soln_Iter"      ) == 0) opt.op.minor_print_level = Nag_Soln_Iter;
					else if (strcmp(tag.szvalue(), "Nag_Soln_Iter_Long" ) == 0) opt.op.minor_print_level = Nag_Soln_Iter_Long;
					else if (strcmp(tag.szvalue(), "Nag_Soln_Iter_Const") == 0) opt.op.minor_print_level = Nag_Soln_Iter_Const;
					else if (strcmp(tag.szvalue(), "Nag_Soln_Iter_Full" ) == 0) opt.op.minor_print_level = Nag_Soln_Iter_Full;
					else throw XMLReader::InvalidValue(tag);
				}
				else if (tag == "outfile"       ) strcpy(opt.op.outfile, tag.szvalue());
				else if (tag == "f_diff_int"    ) tag.value(opt.op.f_diff_int);
				else if (tag == "c_diff_int"    ) tag.value(opt.op.c_diff_int);
				else if (tag == "max_iter"      ) tag.value(opt.op.max_iter);
				else if (tag == "minor_max_iter") tag.value(opt.op.minor_max_iter);
				else if (tag == "f_prec"        ) tag.value(opt.op.f_prec);
				else if (tag == "optim_tol"     ) tag.value(opt.op.optim_tol);
				else if (tag == "linesearch_tol") tag.value(opt.op.linesearch_tol);
				else if (tag == "step_limit"    ) tag.value(opt.op.step_limit);
				else if (tag == "crash_tol"     ) tag.value(opt.op.crash_tol);
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());

	// make sure we have something to optimize
	if ((opt.nfunc == 0) || (opt.nvar == 0) || (opt.pobj == 0)) throw NothingToOptimize();

	return true;	
}


#else

//-----------------------------------------------------------------------------
//! If the NAG libraries are not defined this function simply returns false

bool optimize(FEModel& fem, const char* szfile)
{
	clog.printbox("P A R A M E T E R   O P T I M I Z A T I O N   M O D U L E", "version 0.1");

	// create an optimizer object
	FEOptimizeData opt(fem);

	// read the data from the xml input file
	if (opt.Input(szfile) == false) return false;

	// do initialization
	if (opt.Init() == false) return false;

	// solve the problem
	bool bret = opt.Solve();

	if (bret)
		clog.printf(" N O R M A L   T E R M I N A T I O N\n\n");
	else 
		clog.printf(" E R R O R   T E R M I N A T I O N\n\n");

	return bret;
}

#endif // NAGLIB
