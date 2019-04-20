/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! The FEFluidSolver class solves fluid mechanics problems
//! It can deal with quasi-static and dynamic problems
//!
class FEBIOFLUID_API FEFluidSolver : public FENewtonSolver
{
public:
    //! constructor
    FEFluidSolver(FEModel* pfem);
    
    //! destructor
    ~FEFluidSolver();
    
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
    void AssembleStiffness(vector<int>& en, vector<int>& elmi, vector<int>& elmj, matrix& ke) override;
    
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
    void NodalForces(FEGlobalVector& R, const FETimeInfo& tp);
    
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
    
public:
    // convergence tolerances
    double	m_Rtol;			//!< residual tolerance
    double	m_Vtol;			//!< velocity tolerance
    double	m_Ftol;			//!< dilatation tolerance
    double	m_Etol;			//!< energy tolerance
    double	m_Rmin;			//!< min residual value
	double	m_Rmax;			//!< max residual value
    double  m_minJf;        //!< minimum allowable compression ratio

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
    vector<double> m_di;	//!< dilatation increment vector
    vector<double> m_Di;	//!< Total dilatation vector for iteration
    
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
    DECLARE_FECORE_CLASS();
};
