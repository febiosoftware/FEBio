//
//  FEFluidVSolver.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 1/30/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#ifndef FEFluidVSolver_hpp
#define FEFluidVSolver_hpp

#include <FECore/FENewtonSolver.h>
#include <FECore/FETimeInfo.h>
#include <FECore/FEGlobalVector.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! The FEFluidVSolver class solves fluid mechanics problems
//! It can deal with quasi-static and dynamic problems
//!
class FEBIOFLUID_API FEFluidVSolver : public FENewtonSolver
{
public:
    //! constructor
    FEFluidVSolver(FEModel* pfem);
    
    //! destructor
    ~FEFluidVSolver();
    
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
    
    //! update nodal positions, velocities, accelerations, etc.
    void UpdateKinematics(vector<double>& ui);
    //}
    
    //{ --- Solution functions ---
    
    //! prepares the data for the first QN iteration
    void PrepStep() override;
    
    //! Performs a Newton-Raphson iteration
    bool Quasin() override;
    
    //! Lagrangian augmentation
    bool Augment() override;
    
    //{ --- Stiffness matrix routines ---
    
    //! calculates the global stiffness matrix
    bool StiffnessMatrix() override;
    
    //! contact stiffness
    void ContactStiffness();
    
    //! calculates stiffness contributon of nonlinear constraints
    void NonLinearConstraintStiffness(const FETimeInfo& tp);
    
    //{ --- Residual routines ---
    
    //! Calculates concentrated nodal forces
    void NodalForces(vector<double>& F, const FETimeInfo& tp);
    
    //! Calculate the contact forces
    void ContactForces(FEGlobalVector& R);
    
    //! Calculates residual
    bool Residual(vector<double>& R) override;
    
    //! Calculate nonlinear constraint forces
    void NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp);
    
public:
    // convergence tolerances
    double    m_Rtol;            //!< residual tolerance
    double    m_Vtol;            //!< velocity tolerance
    double    m_Etol;            //!< energy tolerance
    double    m_Rmin;            //!< min residual value
    double    m_Rmax;            //!< max residual value
    
public:
    // equation numbers
    int        m_nveq;                //!< number of equations related to velocity dofs
    
public:
    vector<double> m_Fn;    //!< concentrated nodal force vector
    vector<double> m_Ui;    //!< Total DOF vector for iteration
    vector<double> m_Ut;    //!< Total DOF vector at time t (incl all previous timesteps)
    vector<double> m_Fr;    //!< nodal reaction forces
    
    // generalized alpha method
    double  m_rhoi;         //!< rho infinity
    double  m_alphaf;       //!< alpha step for Y={v,e}
    double  m_alpham;       //!< alpha step for Ydot={∂v/∂t,∂e/∂t}
    double  m_gammaf;       //!< gamma
    int     m_pred;         //!< predictor method
    
protected:
    int     m_dofX;
    int     m_dofY;
    int     m_dofZ;
    
    int     m_dofWX;
    int     m_dofWY;
    int     m_dofWZ;
    
    int     m_dofWXP;
    int     m_dofWYP;
    int     m_dofWZP;
    
    int     m_dofAWX;
    int     m_dofAWY;
    int     m_dofAWZ;
    
    int     m_dofAWXP;
    int     m_dofAWYP;
    int     m_dofAWZP;
    
    // declare the parameter list
    DECLARE_FECORE_CLASS();
};

#endif /* FEFluidVSolver_hpp */
