#pragma once

#include "FECore/FENewtonSolver.h"
#include "FECore/FETypes.h"
#include "FECore/FEGlobalVector.h"

//-----------------------------------------------------------------------------
//! The FEFluidSolver class solves fluid mechanics problems
//! It can deal with quasi-static and dynamic problems
//!
class FEFluidSolver : public FENewtonSolver
{
public:
    //! constructor
    FEFluidSolver(FEModel* pfem);
    
    //! destructor
    ~FEFluidSolver();
    
    //! serialize data to/from dump file
    void Serialize(DumpStream& ar);
    
    //! Initializes data structures
    bool Init();
    
public:
    //! assemble the element residual into the global residual
    //! \todo This was implemented for nodal forces
    void AssembleResidual(int node, int dof, double f, vector<double>& R);
    
    //! adjust the residual matrix for prescribed velocities
    void AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke);
    
    //! assemble global stiffness matrix \todo this is only used by rigid joints
    void AssembleStiffness(vector<int>& elm, matrix& ke);
    
public:
    //{ --- evaluation and update ---
    //! Perform an update
    void Update(vector<double>& ui);
    
    //! Evaluate system, i.e. calculate residual
    void Evaluate(vector<double>& R) { Residual(R); }
    //}
    
    //{ --- Solution functions ---
    
    //! prepares the data for the first QN iteration
    void PrepStep(const FETimeInfo& timeInfo);
    
    //! Performs a Newton-Raphson iteration
    bool Quasin(double time);
    
    //! update nodal positions, velocities, accelerations, etc.
    void UpdateKinematics(vector<double>& ui);
    
    //! Update Stresses
    void UpdateStresses();
    
    //! Lagrangian augmentation
    bool Augment();
    
    //{ --- Stiffness matrix routines ---
    
    //! calculates the global stiffness matrix
    bool StiffnessMatrix(const FETimeInfo& tp);
    
    //! calculates stiffness contributon of nonlinear constraints
    void NonLinearConstraintStiffness(const FETimeInfo& tp);
    
    //{ --- Residual routines ---
    
    //! Calculates concentrated nodal forces
    void NodalForces(vector<double>& F, const FETimeInfo& tp);
    
    //! Calculates residual
    bool Residual(vector<double>& R);
    
    //! Calculate nonlinear constraint forces
    void NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp);
    
public:
    // convergence tolerances
    double	m_Rtol;			//!< residual tolerance
    double	m_Vtol;			//!< velocity tolerance
    double	m_Etol;			//!< energy tolerance
    double	m_Rmin;			//!< min residual value
    
    // strategy parameters
    bool	m_bdivreform;	//!< reform when diverging
    bool	m_bdoreforms;	//!< do reformations
	bool	m_breformtimestep;	//!< reform at start of time step

public:
    vector<double> m_Fn;	//!< concentrated nodal force vector
    vector<double> m_Fr;	//!< nodal reaction forces
    vector<double> m_Vi;	//!< Total velocity vector for iteration
    vector<double> m_Vt;	//!< Total velocity vector at time t (incl all previous timesteps)
    vector<double> m_Fd;	//!< residual correction due to prescribed velocities

public:
    bool		m_baugment;		//!< augmentation flag
    
protected:
	int		m_dofVX;
	int		m_dofVY;
	int		m_dofVZ;
	int		m_dofE;
    
    // declare the parameter list
    DECLARE_PARAMETER_LIST();
};
