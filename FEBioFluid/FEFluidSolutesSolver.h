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

#include <FECore/FENewtonSolver.h>
#include <FECore/FETimeInfo.h>
#include <FECore/FEGlobalVector.h>
#include <FECore/FEDofList.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
class FELinearSystem;

//-----------------------------------------------------------------------------
//! The FEFluidSolutesSolver class solves fluid mechanics problems
//! with mass transport (solutes)
//! It can deal with quasi-static and dynamic problems
//!
class FEBIOFLUID_API FEFluidSolutesSolver : public FENewtonSolver
{
	enum SOLVE_STRATEGY {
		SOLVE_COUPLED,			// monolithic solution approach
		SOLVE_SEQUENTIAL		// first solve velocity, then solutes
	};

public:
    //! constructor
    FEFluidSolutesSolver(FEModel* pfem);
    
    //! destructor
    ~FEFluidSolutesSolver();
    
    //! Initializes data structures
    bool Init() override;
    
    //! initialize the step
    bool InitStep(double time) override;
    
    //! Initialize linear equation system
    bool InitEquations() override;

	//! preferred matrix type should be unsymmetric.
	Matrix_Type PreferredMatrixType() const override { return REAL_UNSYMMETRIC; };
    
public:
    //{ --- evaluation and update ---
    //! Perform an update
    void Update(vector<double>& ui) override;
    
    //! update nodal positions, velocities, accelerations, etc.
    void UpdateKinematics(vector<double>& ui);
    
    //! used by JFNK
    void Update2(const vector<double>& ui) override;
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
    
    //! Serialization
    void Serialize(DumpStream& ar) override;
    
protected:
    void GetVelocityData(vector<double>& vi, vector<double>& ui);
    void GetDilatationData(vector<double>& ei, vector<double>& ui);
    void GetConcentrationData(vector<double>& ci, vector<double>& ui, const int sol);

public:
    // convergence tolerances
    double  m_Vtol;     //!< velocity tolerance
    double  m_Ftol;     //!< dilatation tolerance
    double  m_Ctol;     //!< concentration tolerance
    double  m_minJf;    //!< minimum allowable compression ratio
    bool    m_forcePositive;    //!< force conentrations to remain positive
    double  m_cmin;     //!< threshold for detecting sudden drop in concentration
    double  m_cmax;     //!< threshold for detecting sudden increase in concentration
    int     m_cnum;     //!< minimum number of points at which C drops suddenly

public:
    // equation numbers
    int         m_nveq;     //!< number of equations related to velocity dofs
    int         m_ndeq;     //!< number of equations related to dilatation dofs
	int			m_nseq;		//!< total number of concentration dofs
    vector<int> m_nceq;     //!< number of equations related to concentration dofs

public:
    vector<double> m_Fr;    //!< nodal reaction forces
    vector<double> m_vi;    //!< velocity increment vector
    vector<double> m_Vi;    //!< Total velocity vector for iteration
    vector<double> m_di;    //!< dilatation increment vector
    vector<double> m_Di;    //!< Total dilatation vector for iteration
    
    // solute data
    vector< vector<double> >    m_ci;    //!< concentration increment vector
    vector< vector<double> >    m_Ci;    //!< Total concentration vector for iteration

    // generalized alpha method
    double  m_rhoi;         //!< rho infinity
    double  m_alphaf;       //!< alpha step for Y={v,e}
    double  m_alpham;       //!< alpha step for Ydot={∂v/∂t,∂e/∂t}
    double  m_gammaf;       //!< gamma
    int     m_pred;         //!< predictor method
    
protected:
    FEDofList       m_dofW;
    FEDofList       m_dofAW;
    FEDofList       m_dofEF;
    int             m_dofAEF;
    FEDofList       m_dofC;
    int             m_dofAC;

	int				m_solve_strategy;
    
private:
    bool            m_sudden_C_change;
    
    // declare the parameter list
    DECLARE_FECORE_CLASS();
};
