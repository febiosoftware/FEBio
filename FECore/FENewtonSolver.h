#pragma once
#include "FESolver.h"
#include "FENewtonStrategy.h"
#include "FETypes.h"

//-----------------------------------------------------------------------------
// forward declarations
class FEModel;
class FEGlobalMatrix;

//-----------------------------------------------------------------------------
//! This class defines the base class for Newton-type solvers. 
//! The class implements the basic logic behind a newton-solver but defers some
//! of the logic, especially relating to the update of the stiffness matrix to
//! the FENewtonStrategy class. 
//! \todo there is some common functionality with the FELinearSolver. Perhaps
//! I can make the FELinearSolver a base class?
//! \todo Perhaps I can introduce a base class for linear search algorithms
//! so that the line search strategy can be customized as well.
class FENewtonSolver : public FESolver
{
public:
	//! constructor
	FENewtonSolver(FEModel* pfem);

	//! destrcutor
	~FENewtonSolver();

	//! Set the solution strategy
	void SetSolutionStrategy(FENewtonStrategy* pstrategy);

public: // overloaded from FESolver

	//! Initialization
	bool Init();

	//! Initialize linear equation system
	bool InitEquations();

	//! Clean up
	void Clean();

	//! serialization
	void Serialize(DumpFile& ar);

	//! Solve an analysis step
	bool SolveStep(double time);

public:
	//! return the stiffness matrix
	FEGlobalMatrix& GetStiffnessMatrix();

	//! reform the stiffness matrix
    bool ReformStiffness(const FETimePoint& tp);

    //! recalculates the shape of the stiffness matrix
    bool CreateStiffness(bool breset);
    
protected:
	//! Performs a linesearch
	double LineSearch(double s);

	//! Do a Quasi-Newton step
	//! This is called from SolveStep and must be implemented by derived classes.
	virtual bool Quasin(double time) = 0;

    //! calculates the global stiffness matrix (needs to be overwritten by derived classes)
    virtual bool StiffnessMatrix(const FETimePoint& tp) = 0;

public:
	// line search options
	double	m_LSmin;		//!< minimum line search step
	double	m_LStol;		//!< line search tolerance
	int		m_LSiter;		//!< max nr of line search iterations

	// solver parameters
	int					m_maxref;		//!< max nr of reformations per time step
	FENewtonStrategy*	m_pbfgs;		//!< class handling the specific stiffness update logic

	// linear solver data
	LinearSolver*		m_plinsolve;	//!< the linear solver
	FEGlobalMatrix*		m_pK;			//!< global stiffness matrix
	int					m_neq;			//!< number of equations
    bool				m_breshape;		//!< Matrix reshape flag

	// data used by Quasin
	vector<double> m_R0;	//!< residual at iteration i-1
	vector<double> m_R1;	//!< residual at iteration i
	vector<double> m_ui;	//!< displacement increment vector

	DECLARE_PARAMETER_LIST();
};
