//
//  BFGSSolver2.h
//  FECore
//
//  Created by Gerard Ateshian on 10/12/15.
//  Copyright © 2015 febio.org. All rights reserved.
//

#pragma once

#include "matrix.h"
#include "vector.h"
#include "LinearSolver.h"

class FESolver;

//-----------------------------------------------------------------------------
//! The BFGSSolver solves a nonlinear system of equations using the BFGS method.
//! It depends on the NonLinearSystem to evaluate the function and its jacobian.

class BFGSSolver2
{
public:
    //! constructor
    BFGSSolver2();
    
    //! New initialization method
    void Init(int neq, FESolver* pNLS, LinearSolver* pls);
    
    //! perform a BFGS udpate
    bool Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1);
    
    //! solve the equations
    void SolveEquations(vector<double>& x, vector<double>& b);
    
public:
    int		m_maxups;		//!< max nr of QN iters permitted between stiffness reformations
    double	m_cmax;			//!< maximum value for the condition number
    
public:
    // keep a pointer to the linear solver
    LinearSolver*	m_plinsolve;	//!< pointer to linear solver
    
    // the non-linear system to solve
    FESolver*		m_pNLS;		//!< pointer to nonlinear system to solve
    int				m_neq;		//!< number of equations
    
    // counters
    int		m_nups;			//!< nr of stiffness updates
    
    // BFGS update vectors
    vector<double>	m_dx;       //!< temp vectors for calculating BFGS update vectors
    vector<double>	m_f;        //!< temp vectors for calculating BFGS update vectors
};
/*{
public:
    //! constructor
    BFGSSolver();
    
    //! New initialization method
    void Init(int neq, FESolver* pNLS, LinearSolver* pls);
    
    //! perform a BFGS udpate
    bool Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1);
    
    //! solve the equations
    void SolveEquations(vector<double>& x, vector<double>& b);
    
public:
    int		m_maxups;		//!< max nr of QN iters permitted between stiffness reformations
    double	m_cmax;			//!< maximum value for the condition number
    
public:
    // keep a pointer to the linear solver
    LinearSolver*	m_plinsolve;	//!< pointer to linear solver
    
    // the non-linear system to solve
    FESolver*		m_pNLS;		//!< pointer to nonlinear system to solve
    int				m_neq;		//!< number of equations
    
    // counters
    int		m_nups;			//!< nr of stiffness updates
    
    // BFGS update vectors
    matrix			m_r;		//!< BFGS update vector
    matrix			m_delta;	//!< BFGS update vector
    matrix          m_q;
    vector<double>	m_rho;      //!< temp vectors for calculating BFGS update vectors
};*/
