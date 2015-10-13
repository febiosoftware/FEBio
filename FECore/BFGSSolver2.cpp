//
//  BFGSSolver2.cpp
//  FECore
//
//  Created by Gerard Ateshian on 10/12/15.
//  Copyright © 2015 febio.org. All rights reserved.
//

#include "BFGSSolver2.h"
#include "stdafx.h"
#include "FESolver.h"

//-----------------------------------------------------------------------------
// BFGSSolver2
//-----------------------------------------------------------------------------

BFGSSolver2::BFGSSolver2()
{
    m_maxups = 10;
    m_cmax   = 1e5;
    
    // pointer to linear solver
    m_plinsolve = 0;
    
    // pointer to non-linear system
    m_pNLS = 0;
}

//-----------------------------------------------------------------------------
// Initialization method
void BFGSSolver2::Init(int neq, FESolver* pNLS, LinearSolver* pls)
{
    // allocate storage for BFGS update vectors
    m_dx.resize(neq);
    m_f.resize(neq);
    
    m_neq = neq;
    m_nups = 0;
    
    m_plinsolve = pls;
    
    assert(pNLS);
    m_pNLS = pNLS;
}
/*{
    // allocate storage for BFGS update vectors
    m_r.resize(m_maxups, neq);
    m_delta.resize(m_maxups, neq);
    m_q.resize(m_maxups, neq);
    
    m_rho.resize(neq);
    
    m_neq = neq;
    m_nups = 0;
    
    m_plinsolve = pls;
    
    assert(pNLS);
    m_pNLS = pNLS;
}*/

//-----------------------------------------------------------------------------
//! This function performs a BFGS stiffness update.
//! The last line search step is input to this function.
//! This function returns false if
//! the condition number is too big. A too big condition number might
//! be an indication of an ill-conditioned matrix and the update should
//! not be performed.

bool BFGSSolver2::Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1)
{
    int i;
    
    int neq = (int)ui.size();
    vector<double> df(neq);
    
    for (i=0; i<neq; ++i)
    {
        m_dx[i] = s*ui[i];
        m_f[i] = s*R0[i];
//        df[i] = R1[i] - R0[i];
        df[i] = R1[i] - m_f[i];
    }
    
    double c = sqrt(fabs((m_dx*df)/(m_dx*m_f)));
    
    // make sure c is less than the the maximum.
    if (c >  m_cmax) return false;
    
    // increment update counter
    ++m_nups;
    
    return true;
}
/*{
    int i, j;
    
    int neq = (int)ui.size();
    
    for (j=0; j<m_nups; ++j)
    {
        double dq = 0;
        for (i=0; i<neq; ++i) dq += m_delta[j][i]*m_q[j][i];
        for (i=0; i<neq; ++i)
            m_q[j+i][i] = m_q[j][i] + m_rho[j]*dq*(m_delta[j][i]-m_r[j][i]);
    }
    
    m_rho[m_nups] = 0;
    for (i=0; i<neq; ++i)
    {
        m_r[m_nups][i] = m_q[m_nups][i] + ui[i];
        m_delta[m_nups][i] = s*ui[i];
        m_rho[m_nups] += m_r[m_nups][i]*m_delta[m_nups][i];
    }
    m_rho[m_nups] = 1./m_rho[m_nups];
    
    // increment update counter
    ++m_nups;
    
    return true;
}*/

//-----------------------------------------------------------------------------
// This function solves a system of equations using the BFGS update vectors
// The variable m_nups keeps track of how many updates have been made so far.
//

void BFGSSolver2::SolveEquations(vector<double>& x, vector<double>& b)
{
    int i;
    
    // get the nr of equations
    int neq = (int)x.size();
    
    // make sure we need to do work
    if (neq==0) return;
    
    if (m_nups == 0)
    {
        // create temporary storage
        static vector<double> tmp;
        tmp = b;
        
        // perform a backsubstitution
        m_plinsolve->BackSolve(x, tmp);
    }
    else
    {
        vector<double> df(neq);
        for (i=0; i<neq; ++i) {
            df[i] = b[i] - m_f[i];
        }
        
        double r = (m_dx*b)/(m_dx*df);
        
        for (i=0; i<neq; ++i) {
            x[i] = -r*m_dx[i];
        }
    }
}
/*{
    int i;
    
    // get the nr of equations
    int neq = (int)x.size();
    
    // make sure we need to do work
    if (neq==0) return;
    
    if (m_nups == 0)
    {
        // create temporary storage
        static vector<double> tmp;
        tmp = b;
        
        // perform a backsubstitution
        m_plinsolve->BackSolve(x, tmp);
    }
    else
    {
        double dq = 0;
        for (i=0; i<neq; ++i) dq += m_delta[m_nups][i]*m_q[m_nups][i];
        
        for (i=0; i<neq; ++i)
            x[i] = -m_q[m_nups][i] - m_rho[m_nups]*dq*(m_delta[m_nups][i]-m_r[m_nups][i]);
    }
}*/
