#pragma once
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
#define QN_BFGS		0
#define QN_BFGS2	1
#define QN_BROYDEN	2

//-----------------------------------------------------------------------------
class LinearSolver;

//-----------------------------------------------------------------------------
//! A Base class for newton-type solution strategies
class FENewtonStrategy
{
public:
	FENewtonStrategy();
	virtual ~FENewtonStrategy();

public:
	//! Data allocation and initialization
	virtual void Init(int neq, LinearSolver* pls) = 0;

	//! perform a Newton udpate
	virtual bool Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1) = 0;

	//! solve the equations
	virtual void SolveEquations(vector<double>& x, vector<double>& b) = 0;

public:
	int		m_maxups;		//!< max nr of QN iters permitted between stiffness reformations
	double	m_cmax;			//!< maximum value for the condition number
	int		m_nups;			//!< nr of stiffness updates
};
