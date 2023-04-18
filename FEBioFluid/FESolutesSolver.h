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

class FESolutesSolver : public FENewtonSolver
{
public:
    //! constructor
	FESolutesSolver(FEModel* pfem);
    
    //! destructor
    ~FESolutesSolver();
    
    //! Initializes data structures
    bool Init() override;
    
    //! initialize the step
    bool InitStep(double time) override;
    
    //! Initialize linear equation system
    bool InitEquations() override;
	bool InitEquations2() override;

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
    
    //{ --- Residual routines ---
    
    //! Calculates residual
    bool Residual(vector<double>& R) override;

public:
    //! Serialization
    void Serialize(DumpStream& ar) override;
    
protected:
    void GetConcentrationData(vector<double>& ci, vector<double>& ui, const int sol);

public:
    // convergence tolerances
    double  m_Ctol;     //!< concentration tolerance
    bool    m_forcePositive;    //!< force conentrations to remain positive

public:
    // equation numbers
    vector<int> m_nceq;     //!< number of equations related to concentration dofs
	int			m_nCeq;		//!< total number of concentration dofs

public:
    vector<double> m_Fr;    //!< nodal reaction forces
    
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
    FEDofList	m_dofC;
    FEDofList	m_dofAC;
    
    // declare the parameter list
    DECLARE_FECORE_CLASS();
};
