#pragma once
#include "FEBioMech/FESolidSolver2.h"
#include "FECore/FEElementTraits.h"

//-----------------------------------------------------------------------------
// This class adds additional functionality to the FESolidSolver to solve
// biphasic problems. 
class FEBiphasicSolver : public FESolidSolver2
{
public:
	//! constructor
	FEBiphasicSolver(FEModel* pfem);
	virtual ~FEBiphasicSolver() {}

	//! Initialize data structures
	bool Init();

	//! Initialize linear equation system
	bool InitEquations();

	//! prepares the data for the first QN iteration
	virtual void PrepStep(double time);

	//! Performs a Newton-Raphson iteration
	bool Quasin(double time);

	//! serialize data to/from dump file
	void Serialize(DumpFile& ar);

public:
	//! update contact
	virtual void UpdateContact();

	//! update kinematics
	virtual void UpdateKinematics(vector<double>& ui);

	//! Update poroelastic data
	void UpdatePoro(vector<double>& ui);

public:

	//! Calculates concentrated nodal forces (overridden from FESolidSolver2)
	//! (This function is called from FESolidSolver2::PrepStep)
	virtual void NodalForces(vector<double>& F, const FETimePoint& tp);

	//! Calculates residual (overridden from FESolidSolver2)
	virtual bool Residual(vector<double>& R);

	//! calculates the global stiffness matrix (overridden from FESolidSolver2)
	virtual bool StiffnessMatrix(const FETimePoint& tp);

protected:
	void GetDisplacementData(vector<double>& di, vector<double>& ui);
	void GetPressureData(vector<double>& pi, vector<double>& ui);

public:
	// additional convergence norms
	double	m_Ptol;			//!< pressure tolerance

	// equation numbers
	int		m_ndeq;				//!< number of equations related to displacement dofs
	int		m_npeq;				//!< number of equations related to pressure dofs
	vector<int>		m_nceq;	//!< number of equations related to concentration dofs

	// poro data
	vector<double>	m_di;	//!< displacement increment vector
	vector<double>	m_Di;	//!< total displacement vector for iteration
	vector<double>	m_pi;	//!< pressure increment vector
	vector<double>	m_Pi;	//!< Total pressure vector for iteration

protected:
	int	m_dofP;	//!< pressure dof index

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
