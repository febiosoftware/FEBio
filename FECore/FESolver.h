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



#pragma once
#include "FECoreBase.h"
#include "Timer.h"
#include "matrix.h"
#include "vector.h"
#include "FEDofList.h"
#include "FETimeInfo.h"

//-----------------------------------------------------------------------------
// Scheme for assigning equation numbers
// STAGGERED: | a0, b0, a1, b1, ..., an, bn |
// BLOCK    : | a0, a1, ..., an, b0, b1, ..., bn |
enum EQUATION_SCHEME
{
	STAGGERED,
	BLOCK
};

//-----------------------------------------------------------------------------
enum EQUATION_ORDER
{
	NORMAL_ORDER,
	REVERSE_ORDER,
	FEBIO2_ORDER
};

//-----------------------------------------------------------------------------
// Solution variable
class FESolutionVariable
{
public:
	FESolutionVariable(const char* szname, FEDofList* dofs = nullptr, int order = 2)
	{
		m_szname = szname;
		m_dofs = dofs;
		m_order = order;
	}

public:
	FEDofList*	m_dofs;		// the dof list
	int			m_order;	// the order of interpolation (0 = constant, 1 = linear, 2 = default)
	const char* m_szname;	// name of solution variable
};

//-----------------------------------------------------------------------------
// structure identifying nodal dof info
struct FECORE_API FENodalDofInfo
{
	int		m_eq = -1;		// equation number
	int		m_node = -1;		// 0-based index into mesh!
	int		m_dof = -1;		// index into nodal m_ID array
	const char* szdof = nullptr;
};

//-----------------------------------------------------------------------------
class FEModel;
class FEGlobalMatrix;
class LinearSolver;
class FEGlobalVector;

//-----------------------------------------------------------------------------
//! This is the base class for all FE solvers.

//! A class derived from FESolver implements a solver for a specific type
//! of physics problem. It takes the FEModel in its constructor and implements
//! the SolveStep function to solve the FE problem.
class FECORE_API FESolver : public FECoreBase
{
	FECORE_SUPER_CLASS(FESOLVER_ID)
	FECORE_BASE_CLASS(FESolver)

public:
	//! constructor
	FESolver(FEModel* fem);

	//! destructor
	virtual ~FESolver();

public:
	//! Data serialization
	void Serialize(DumpStream& ar) override;

	//! This is called by FEAnalaysis::Deactivate
	virtual void Clean();

	//! rewind the solver (This is called when the time step fails and needs to retry)
	virtual void Rewind() {}

	//! called during model reset
	virtual void Reset();

	// Initialize linear equation system
	virtual bool InitEquations();

	// New equation initialization procedure
	// TODO: work in progress
	virtual bool InitEquations2();

	//! add equations
	void AddEquations(int neq, int partition = 0);

	//! initialize the step (This is called before SolveStep)
	virtual bool InitStep(double time);

	//! Solve an analysis step
	virtual bool SolveStep() = 0;

	//! Update the state of the model
	virtual void Update(std::vector<double>& u);

	//! Do the augmentations
	virtual bool Augment();

	//! Calculates concentrated nodal loads
//	virtual void NodalLoads(FEGlobalVector& R, const FETimeInfo& tp);

public:
	//! Set the equation allocation scheme
	void SetEquationScheme(int scheme);

	//! set the linear system partitions
	void SetPartitions(const vector<int>& part);

	//! Get the size of a partition
	int GetPartitionSize(int partition);

	//! get the current stiffness matrix
	virtual FEGlobalMatrix* GetStiffnessMatrix();

	//! get the current load vector
	virtual std::vector<double> GetLoadVector();

	// get the linear solver
	virtual LinearSolver* GetLinearSolver();

	//! Matrix symmetry flag
	int MatrixSymmetryFlag() const;

	//! get matrix type
	Matrix_Type MatrixType() const;

	//! build the matrix profile
	virtual void BuildMatrixProfile(FEGlobalMatrix& G, bool breset);

	// see if the dofs in the dof list are active in this solver
	bool HasActiveDofs(const FEDofList& dof);

	// get the active dof map (returns nr of functions)
	int GetActiveDofMap(vector<int>& dofMap);

	// return the node (mesh index) from an equation number
	FENodalDofInfo GetDOFInfoFromEquation(int ieq);

public:
	// extract the (square) norm of a solution vector
	double ExtractSolutionNorm(const vector<double>& v, const FEDofList& dofs) const;

	// return the solution vector
	virtual std::vector<double> GetSolutionVector() const;

public: //TODO Move these parameters elsewhere
	bool				m_bwopt;	    //!< bandwidth optimization flag
	int					m_msymm;		//!< matrix symmetry flag for linear solver allocation
	int					m_eq_scheme;	//!< equation number scheme (used in InitEquations)
	int					m_eq_order;		//!< normal or reverse ordering
	int					m_neq;			//!< number of equations
	std::vector<int>	m_part;			//!< partitions of linear system
	std::vector<int>	m_dofMap;		//!< array stores for each equation the corresponding dof index

	// counters
	int		m_nrhs;			//!< nr of right hand side evalutations
	int		m_niter;		//!< nr of quasi-newton iterations
	int		m_nref;			//!< nr of stiffness retormations
	int		m_ntotref;		//!< nr of total stiffness reformations

	// augmentation
	int		m_naug;			//!< nr of augmentations
	bool	m_baugment;		//!< do augmentations flag

protected:
	// list of solution variables
	vector<FESolutionVariable>	m_Var;

	DECLARE_FECORE_CLASS();
};
