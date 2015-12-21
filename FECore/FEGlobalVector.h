#pragma once
#include <vector>
using namespace std;

class FEModel;

//-----------------------------------------------------------------------------
//! This class represents a global system array. It provides functions to assemble
//! local (element) vectors into this array
class FEGlobalVector
{
public:
	//! constructor
	FEGlobalVector(FEModel& fem, vector<double>& R, vector<double>& Fr);

	//! destructor
	virtual ~FEGlobalVector();

	//! Assemble the element vector into this global vector
	virtual void Assemble(vector<int>& en, vector<int>& elm, vector<double>& fe);

	//! Assemble into this global vector
	virtual void Assemble(vector<int>& lm, vector<double>& fe);
    
	//! access operator
	double& operator [] (int i) { return m_R[i]; }

	//! Get the FE model
	FEModel& GetFEModel() { return m_fem; }

protected:
	FEModel&			m_fem;	//!< model
	vector<double>&		m_R;	//!< residual
	vector<double>&		m_Fr;	//!< nodal reaction forces \todo I want to remove this
};
