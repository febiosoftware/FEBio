#pragma once

#include "FEBioXML/XMLReader.h"
#include "FECore/FEModel.h"
#include "NumCore/vector.h"
#include <string.h>

//-----------------------------------------------------------------------------
// IO exceptions

//! the variable name is not recognized
class InvalidVariableName
{
public:
	InvalidVariableName(const char* sz) { strcpy(szname, sz); }
	char szname[256];
};

//! there is nothing to optimize
class NothingToOptimize{};

//! FEBio error terminated during the optimization
class FEErrorTermination{};

//-----------------------------------------------------------------------------
class FEOptimizeData;

//=============================================================================
//! Class that reads the optimization input file
class FEOptimizeInput
{
public:
	bool Input(const char* szfile, FEOptimizeData* pOpt);

protected:
	bool ParseModel     (XMLTag& tag, FEOptimizeData& opt);
	bool ParseOptions   (XMLTag& tag, FEOptimizeData& opt);
	bool ParseObjective (XMLTag& tag, FEOptimizeData& opt);
	bool ParseParameters(XMLTag& tag, FEOptimizeData& opt);
	bool ParseLoadData  (XMLTag& tag, FEOptimizeData& opt);
};

//=============================================================================
//! optimization method - this class does the actual work
class FEOptimizeMethod
{
public:
	virtual bool Solve(FEOptimizeData* pOpt) = 0;
};

//----------------------------------------------------------------------------
//! Optimization method using Powell's method
class FEPowellOptimizeMethod : public FEOptimizeMethod
{
public:
	bool Solve(FEOptimizeData* pOpt);

protected:
	double ObjFun(double* p);

	FEOptimizeData*	m_pOpt;
	
	static FEPowellOptimizeMethod*	m_pThis;
	static double objfun(double* p) { return (m_pThis)->ObjFun(p); }
};

//----------------------------------------------------------------------------
//! Optimization method using Levenberg-Marquardt method
class FELMOptimizeMethod : public FEOptimizeMethod
{
public:
	FELMOptimizeMethod();
	bool Solve(FEOptimizeData* pOpt);

protected:
	FEOptimizeData* m_pOpt;

	void ObjFun(vector<double>& x, vector<double>& a, vector<double>& y, matrix& dyda);

	bool FESolve(vector<double>& x, vector<double>& a, vector<double>& y);

	static FELMOptimizeMethod* m_pThis;
	static void objfun(vector<double>& x, vector<double>& a, vector<double>& y, matrix& dyda) { return m_pThis->ObjFun(x, a, y, dyda); }

public:
	double	m_objtol;	// objective tolerance
	double	m_fdiff;	// forward difference step size

protected:
	vector<double>	m_yopt;	// optimal y-values
	vector<double>	m_y0;	// initial (target) y-values
};

//----------------------------------------------------------------------------
//! optimization method using the NAG library
class FENAGOptimizeMethod : public FEOptimizeMethod
{
public:
	bool Solve(FEOptimizeData* pOpt);
};

//=============================================================================
struct OPT_VARIABLE
{
	char	m_szname[128];	//!< variable name
	double*	m_pd;			//!< pointer to variable data
	double	m_val;			//!< value
	double	m_min, m_max;	//!< variable bounds
	double	m_sf;			//!< variable scale factor
};

//=============================================================================
struct OPT_OBJECTIVE
{
	char	m_szname[128];	//!< name of objective
	double*	m_pd;			//!< pointer to variable data
	int		m_nlc;			//!< load curve
};

//=============================================================================
//! optimization analyses
//! 
class FEOptimizeData
{
public:
	//! constructor
	FEOptimizeData(FEModel& fem);
	~FEOptimizeData(void);

	//! input function
	bool Input(const char* sz);

	//! Initialize data
	bool Init();

	//! solver the problem
	bool Solve();

	//! return the FE Model
	FEModel& GetFEM() { return m_fem; }

	//! add a loadcurve
	void AddLoadCurve(FELoadCurve* plc) { m_LC.push_back(plc); }

	//! set the objective function
	void SetObjective(OPT_OBJECTIVE o) { m_obj = o; }

	//! add a variable to optimize
	void AddVariable(OPT_VARIABLE& var) { m_Var.push_back(var); }

	int Variables() { return m_Var.size(); }

	OPT_VARIABLE& Variable(int n) { return m_Var[n]; }

	OPT_OBJECTIVE& GetObjective() { return m_obj; }

	FELoadCurve& ReactionLoad() { return m_rf; }

	FELoadCurve& GetLoadCurve(int n) { return *m_LC[n]; }

	void SetSolver(FEOptimizeMethod* po) { m_pSolver = po; }

public:
	int	m_niter;	// nr of minor iterations (i.e. FE solves)

protected:
	FEModel&	m_fem;

	OPT_OBJECTIVE	m_obj;		//!< the objective function

	FEOptimizeMethod*	m_pSolver;

	FELoadCurve	m_rf;	// reaction force data

	vector<FELoadCurve*>	m_LC;
	vector<OPT_VARIABLE>	m_Var;
};
