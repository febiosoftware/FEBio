#pragma once
#include "FELinearConstraint.h"
#include "table.h"

class FEGlobalMatrix;
class matrix;

//-----------------------------------------------------------------------------
// This class helps manage all the linear constraints
class FECORE_API FELinearConstraintManager
{
public:
	FELinearConstraintManager(FEModel* fem);

	// Clear all constraints
	void Clear();

	// copy data
	void CopyFrom(const FELinearConstraintManager& lcm);

	// serialize linear constraints
	void Serialize(DumpStream& ar);

	// build the matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& G);

	// add the linear constraint
	void AddLinearConstraint(FELinearConstraint& lc);

	// return number of linear constraints
	int LinearConstraints() const;

	// return linear constraint
	const FELinearConstraint& LinearConstraint(int i) const;

	// return linear constraint
	FELinearConstraint& LinearConstraint(int i);

	//! remove a linear constraint
	void RemoveLinearConstraint(int i);

public:
	// one-time initialization
	bool Initialize();

	// activation
	bool Activate();

	// assemble element residual into global residual
	void AssembleResidual(vector<double>& R, vector<int>& en, vector<int>& elm, vector<double>& fe);

	// assemble element matrix into (reduced) global matrix
	void AssembleStiffness(FEGlobalMatrix& K, vector<double>& R, vector<double>& ui, vector<int>& en, vector<int>& elm, matrix& ke);

	// called before the first reformation for each time step
	void PrepStep();

	// update nodal variables
	void Update();

protected:
	void InitTable();

private:
	FEModel* m_fem;
	vector<FELinearConstraint>	m_LinC;		//!< linear constraints data
	table<int>					m_LCT;		//!< linear constraint table
	vector<double>				m_up;		//!< the inhomogenous component of the linear constraint
};
