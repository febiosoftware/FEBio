/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2022 University of Utah, The Trustees of Columbia University in
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
#include <FECore/FENewtonSolver.h>
#include <FECore/FETimeInfo.h>
#include <FECore/FEGlobalVector.h>
#include <FECore/FEDofList.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
class FELinearSystem;

//-----------------------------------------------------------------------------
//! The FEPolarFluidSolver class solves polar fluid problems
//! It can deal with quasi-static and dynamic problems
//!
class FEBIOFLUID_API FEPolarFluidSolver : public FENewtonSolver
{
public:
    //! constructor
    FEPolarFluidSolver(FEModel* pfem);
    
    //! destructor
    ~FEPolarFluidSolver();
    
    //! serialize data to/from dump file
    void Serialize(DumpStream& ar) override;
    
    //! Initializes data structures
    bool Init() override;
    
    //! initialize the step
    bool InitStep(double time) override;
    
    //! Initialize linear equation system
    bool InitEquations() override;
    bool InitEquations2() override;

    //! Generate warnings if needed
    void SolverWarnings();
    
public:
    //{ --- evaluation and update ---
    //! Perform an update
    void Update(vector<double>& ui) override;
    void Update2(const vector<double>& ui) override;

    //! update nodal positions, velocities, accelerations, etc.
    void UpdateKinematics(vector<double>& ui);
    
    void UpdateConstraints();
    
    //! update DOF increments
    virtual void UpdateIncrements(vector<double>& Ui, vector<double>& ui, bool emap);
    //}
    
    //{ --- Solution functions ---
    
    //! prepares the data for the first QN iteration
    void PrepStep() override;
    
    //! Performs a Newton-Raphson iteration
    bool Quasin() override;
    
    //{ --- Stiffness matrix routines ---
    
    //! calculates the global stiffness matrix
    bool StiffnessMatrix(FELinearSystem& LS) override;

    //! contact stiffness
    void ContactStiffness(FELinearSystem& LS);
    
    //! calculates stiffness contributon of nonlinear constraints
    void NonLinearConstraintStiffness(FELinearSystem& LS, const FETimeInfo& tp);
    
    //{ --- Residual routines ---
    
    //! Calculate the contact forces
    void ContactForces(FEGlobalVector& R);
    
    //! Calculates residual
    bool Residual(vector<double>& R) override;
    
    //! Calculate nonlinear constraint forces
    void NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp);
    
protected:
    void GetVelocityData(vector<double>& vi, vector<double>& ui);
    void GetAngularVelocityData(vector<double>& xi, vector<double>& ui);
    void GetDilatationData(vector<double>& ei, vector<double>& ui);
    
public:
    // convergence tolerances
    double  m_Vtol;     //!< velocity tolerance
    double  m_Gtol;     //!< angular velocity tolerance
    double  m_Ftol;     //!< dilatation tolerance
    double  m_minJf;    //!< minimum allowable compression ratio
    
public:
    // equation numbers
    int m_nveq;     //!< number of equations related to velocity dofs
    int m_ngeq;     //!< number of equations related to angular velocity dofs
    int m_nfeq;     //!< number of equations related to dilatation dofs
    
public:
    vector<double> m_Fn;    //!< concentrated nodal force vector
    vector<double> m_Fr;    //!< nodal reaction forces
    vector<double> m_vi;    //!< velocity increment vector
    vector<double> m_Vi;    //!< Total velocity vector for iteration
    vector<double> m_gi;    //!< angular velocity increment vector
    vector<double> m_Gi;    //!< Total angular velocity vector for iteration
    vector<double> m_fi;    //!< dilatation increment vector
    vector<double> m_Fi;    //!< Total dilatation vector for iteration
    
    // generalized alpha method
    double  m_rhoi;         //!< spectral radius (rho infinity)
    double  m_alphaf;       //!< alpha step for Y={v,e}
    double  m_alpham;       //!< alpha step for Ydot={∂v/∂t,∂e/∂t}
    double  m_beta;         //!< beta
    double  m_gamma;        //!< gamma
    int     m_pred;         //!< predictor method
    int     m_order;        //!< generalized-alpha integration order
    
protected:
    FEDofList   m_dofW;     // fluid velocity
    FEDofList   m_dofAW;    // material time derivative of fluid velocity
    FEDofList   m_dofG;     // fluid angular velocity
    FEDofList   m_dofAG;    // material time derivative of fluid angular velocity
    FEDofList   m_dofEF;    // fluid dilatation
    int         m_dofAEF;   // material time derivative of fluid dilatation
    
    // declare the parameter list
    DECLARE_FECORE_CLASS();
};
