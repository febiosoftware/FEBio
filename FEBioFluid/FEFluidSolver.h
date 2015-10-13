#pragma once

#include "FECore/FENewtonSolver2.h"
#include "FECore/FETypes.h"

//-----------------------------------------------------------------------------
//! The FEFluidSolver class solves fluid mechanics problems
//! It can deal with quasi-static and dynamic problems
//!
class FEFluidSolver : public FENewtonSolver2
{
public:
    //! constructor
    FEFluidSolver(FEModel* pfem);
    
    //! destructor
    ~FEFluidSolver();
    
    //! serialize data to/from dump file
    void Serialize(DumpFile& ar);
    
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
    void PrepStep(double time);
    
    //! Performs a Newton-Raphson iteration
    bool Quasin(double time);
    
    //! update nodal positions, velocities, accelerations, etc.
    void UpdateKinematics(vector<double>& ui);
    
    //! Update Stresses
    void UpdateStresses();
    
    //{ --- Stiffness matrix routines ---
    
    //! calculates the global stiffness matrix
    bool StiffnessMatrix(const FETimePoint& tp);
    
    //{ --- Residual routines ---
    
    //! Calculates concentrated nodal forces
    void NodalForces(vector<double>& F, const FETimePoint& tp);
    
    //! Calculates residual
    bool Residual(vector<double>& R);
    
public:
    // convergence tolerances
    double	m_Rtol;			//!< residual tolerance
    double	m_Vtol;			//!< velocity tolerance
    double	m_Rmin;			//!< min residual value
    
    // strategy parameters
    bool	m_bdivreform;	//!< reform when diverging
    bool	m_bdoreforms;	//!< do reformations
    
public:
    vector<double> m_Fn;	//!< concentrated nodal force vector
    vector<double> m_Fr;	//!< nodal reaction forces
    vector<double> m_Vi;	//!< Total velocity vector for iteration
    vector<double> m_Vt;	//!< Total velocity vector at time t (incl all previous timesteps)
    vector<double> m_Fd;	//!< residual correction due to prescribed velocities
    
    // declare the parameter list
    DECLARE_PARAMETER_LIST();
};
