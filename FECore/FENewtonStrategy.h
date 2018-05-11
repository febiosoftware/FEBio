#pragma once
#include <vector>
#include "LinearSolver.h"
using namespace std;

//-----------------------------------------------------------------------------
class FENewtonSolver;

//-----------------------------------------------------------------------------
//! A Base class for newton-type solution strategies
class FENewtonStrategy
{
public:
	FENewtonStrategy(FENewtonSolver* pns);
	virtual ~FENewtonStrategy();

public:
	//! Data allocation and initialization
	virtual void Init(int neq, LinearSolver* pls) = 0;

	//! initialize the linear system
	virtual SparseMatrix* CreateSparseMatrix(Matrix_Type mtype);

	//! Presolve update
	virtual void PreSolveUpdate() {}

	//! perform a Newton udpate
	virtual bool Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1) = 0;

	//! solve the equations
	virtual void SolveEquations(vector<double>& x, vector<double>& b) = 0;

	//! reform the stiffness matrix
	virtual bool ReformStiffness();

public:
	int		m_maxups;		//!< max nr of QN iters permitted between stiffness reformations
	int		m_max_buf_size;	//!< max buffer size for update vector storage
	bool	m_cycle_buffer;	//!< recycle the buffer when updates is larger than buffer size
	double	m_cmax;			//!< maximum value for the condition number
	int		m_nups;			//!< nr of stiffness updates

protected:
	FENewtonSolver*	m_pns;
};
