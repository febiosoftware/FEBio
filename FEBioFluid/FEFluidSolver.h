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
    void Serialize(DumpStream& ar) override;
    
    //! Initializes data structures
    bool Init() override;
    
    //! initialize the step
    bool InitStep(double time) override;
    
    //! Initialize linear equation system
    bool InitEquations() override;
    
public:
    //! assemble the element residual into the global residual
    //! \todo This was implemented for nodal forces
    void AssembleResidual(int node, int dof, double f, vector<double>& R);
    
    //! adjust the residual matrix for prescribed velocities
    void AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke) override;
    
    //! assemble global stiffness matrix \todo this is only used by rigid joints
    void AssembleStiffness(vector<int>& elm, matrix& ke) override;
    
public:
    //{ --- evaluation and update ---
    //! Perform an update
    void Update(vector<double>& ui) override;
    
    //}
    
    //{ --- Solution functions ---
    
    //! prepares the data for the first QN iteration
    void PrepStep(const FETimeInfo& timeInfo);
    
    //! Performs a Newton-Raphson iteration
    bool Quasin(double time) override;
    
    //! update nodal positions, velocities, accelerations, etc.
    void UpdateKinematics(vector<double>& ui);
    
    //! Update Stresses
    void UpdateStresses();
    
    //! Lagrangian augmentation
    bool Augment();
    
    //{ --- Stiffness matrix routines ---
    
    //! calculates the global stiffness matrix
    bool StiffnessMatrix(const FETimeInfo& tp) override;
    
    //! calculates stiffness contributon of nonlinear constraints
    void NonLinearConstraintStiffness(const FETimeInfo& tp);
    
    //{ --- Residual routines ---
    
    //! Calculates concentrated nodal forces
    void NodalForces(vector<double>& F, const FETimeInfo& tp);
    
    //! Calculates residual
    bool Residual(vector<double>& R) override;
    
    //! Calculate nonlinear constraint forces
    void NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp);
    
protected:
    void GetVelocityData(vector<double>& vi, vector<double>& ui);
    void GetDilatationData(vector<double>& ei, vector<double>& ui);
    
public:
    // convergence tolerances
    double	m_Rtol;			//!< residual tolerance
    double	m_Vtol;			//!< velocity tolerance
    double	m_Ftol;			//!< dilatation tolerance
    double	m_Etol;			//!< energy tolerance
    double	m_Rmin;			//!< min residual value
    
    // strategy parameters
    bool	m_bdivreform;	//!< reform when diverging
    bool	m_bdoreforms;	//!< do reformations
	bool	m_breformtimestep;	//!< reform at start of time step

public:
    // equation numbers
    int		m_nveq;				//!< number of equations related to velocity dofs
    int		m_ndeq;				//!< number of equations related to dilatation dofs
    
public:
    vector<double> m_Fn;	//!< concentrated nodal force vector
    vector<double> m_Ui;	//!< Total DOF vector for iteration
    vector<double> m_Ut;	//!< Total DOF vector at time t (incl all previous timesteps)
    vector<double> m_Fr;	//!< nodal reaction forces
    vector<double> m_vi;	//!< velocity increment vector
    vector<double> m_Vi;	//!< Total velocity vector for iteration
    vector<double> m_Fd;	//!< residual correction due to prescribed velocities
    vector<double> m_di;	//!< dilatation increment vector
    vector<double> m_Di;	//!< Total dilatation vector for iteration
    
    // generalized alpha method
    double  m_rhoi;         //!< rho infinity
    double  m_alphaf;       //!< alpha step for Y={v,e}
    double  m_alpham;       //!< alpha step for Ydot={∂v/∂t,∂e/∂t}
    double  m_gammaf;       //!< gamma
    int     m_pred;         //!< predictor method

public:
    bool		m_baugment;		//!< augmentation flag
    
protected:
    int     m_dofX;
    int     m_dofY;
    int     m_dofZ;
    
	int		m_dofWX;
	int		m_dofWY;
	int		m_dofWZ;
	int		m_dofEF;
    
    int		m_dofWXP;
    int		m_dofWYP;
    int		m_dofWZP;
    int     m_dofEFP;
    
    int		m_dofAWX;
    int		m_dofAWY;
    int		m_dofAWZ;
    int     m_dofAEF;
    
    int		m_dofAWXP;
    int		m_dofAWYP;
    int		m_dofAWZP;
    int     m_dofAEFP;
    
    // declare the parameter list
    DECLARE_PARAMETER_LIST();
};
