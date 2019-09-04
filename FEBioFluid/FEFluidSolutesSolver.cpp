/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2019 University of Utah, The Trustees of Columbia University in
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



#include "stdafx.h"
#include "FEFluidSolutesSolver.h"
#include "FEFluidDomain.h"
#include "FEFluidResidualVector.h"
#include "FECore/FEModel.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include "NumCore/NumCore.h"
#include <assert.h>
#include "FECore/FEGlobalMatrix.h"
#include "FECore/sys.h"
#include <FEBioMech/FEBodyForce.h>
#include <FECore/FEBoundaryCondition.h>
#include <FECore/FENodalLoad.h>
#include <FECore/FESurfaceLoad.h>
#include "FEFluidResistanceBC.h"
#include "FEBackFlowStabilization.h"
#include "FEFluidNormalVelocity.h"
#include "FEFluidVelocity.h"
#include "FEFluidRotationalVelocity.h"
#include "FETiedFluidInterface.h"
#include <FECore/FEModelLoad.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FELinearConstraintManager.h>
#include <FECore/FELinearSystem.h>
#include "FEBioFluidSolutes.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_FECORE_CLASS(FEFluidSolutesSolver, FENewtonSolver)
ADD_PARAMETER(m_Vtol , "vtol"        );
ADD_PARAMETER(m_Ftol , "ftol"        );
ADD_PARAMETER(m_Ctol , "ctol"        );
ADD_PARAMETER(m_rhoi , "rhoi"        );
ADD_PARAMETER(m_pred , "predictor"   );
ADD_PARAMETER(m_minJf, "min_volume_ratio");
ADD_PARAMETER(m_forcePositive, "force_positive_concentrations");
ADD_PARAMETER(m_solve_strategy, "solve_strategy");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! FEFluidSolutesSolver Construction
//
FEFluidSolutesSolver::FEFluidSolutesSolver(FEModel* pfem) : FENewtonSolver(pfem), m_dofW(pfem), m_dofAW(pfem)
{
    // default values
    m_Rtol = 0.001;
    m_Etol = 0.01;
    m_Vtol = 0.001;
    m_Ftol = 0.001;
    m_Ctol = 0.01;
    m_Rmin = 1.0e-20;
    m_Rmax = 0;     // not used if zero
    m_minJf = 0;    // not used if zero
    
    m_nveq = 0;
    m_ndeq = 0;
    m_niter = 0;
    
    // assume non-symmetric stiffness
    m_msymm = REAL_UNSYMMETRIC;
    
    m_forcePositive = true;    // force all concentrations to remain positive

	m_solve_strategy = SOLVE_COUPLED;
    
    m_rhoi = 0;
    m_pred = 0;
    
    // Preferred strategy is Broyden's method
    SetDefaultStrategy(QN_BROYDEN);
    
    // turn off checking for a zero diagonal
    CheckZeroDiagonal(false);
    
    // Allocate degrees of freedom
    DOFS& dofs = pfem->GetDOFS();
    int varD = dofs.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::DISPLACEMENT), VAR_VEC3);
    dofs.SetDOFName(varD, 0, "x");
    dofs.SetDOFName(varD, 1, "y");
    dofs.SetDOFName(varD, 2, "z");
    
    int nW = dofs.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::RELATIVE_FLUID_VELOCITY), VAR_VEC3);
    dofs.SetDOFName(nW, 0, "wx");
    dofs.SetDOFName(nW, 1, "wy");
    dofs.SetDOFName(nW, 2, "wz");
    
    int nE = dofs.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_DILATATION), VAR_SCALAR);
    dofs.SetDOFName(nE, 0, "ef");
    
    int nAW = dofs.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::RELATIVE_FLUID_ACCELERATION), VAR_VEC3);
    dofs.SetDOFName(nAW, 0, "awx");
    dofs.SetDOFName(nAW, 1, "awy");
    dofs.SetDOFName(nAW, 2, "awz");
    
    int nAE = dofs.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_DILATATION_TDERIV), VAR_SCALAR);
    dofs.SetDOFName(nAE, 0, "aef");
    
    int varC = dofs.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION), VAR_ARRAY);
    int varAC = dofs.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION_TDERIV), VAR_ARRAY);

    // get the dof indices
    m_dofW.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::RELATIVE_FLUID_VELOCITY));
    m_dofAW.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::RELATIVE_FLUID_ACCELERATION));
    m_dofEF  = pfem->GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_DILATATION), 0);
    m_dofAEF = pfem->GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_DILATATION_TDERIV), 0);
    m_dofC = m_dofAC = -1;
}

//-----------------------------------------------------------------------------
FEFluidSolutesSolver::~FEFluidSolutesSolver()
{
    
}

//-----------------------------------------------------------------------------
//! Allocates and initializes the data structures used by the FEFluidSolutesSolver
//
bool FEFluidSolutesSolver::Init()
{
    // initialize base class
    if (FENewtonSolver::Init() == false) return false;
    
    // check parameters
    if (m_Vtol <  0.0) { feLogError("vtol must be nonnegative."); return false; }
    if (m_Ftol <  0.0) { feLogError("dtol must be nonnegative."); return false; }
    if (m_Ctol <  0.0) { feLogError("ctol must be nonnegative."); return false; }
    if (m_Etol <  0.0) { feLogError("etol must be nonnegative."); return false; }
    if (m_Rtol <  0.0) { feLogError("rtol must be nonnegative."); return false; }
    
    if (m_rhoi == -1) {
        m_alphaf = m_alpham = m_gammaf = 1.0;
    }
    else if ((m_rhoi >= 0) && (m_rhoi <= 1)) {
        m_alphaf = 1.0/(1+m_rhoi);
        m_alpham = (3-m_rhoi)/(1+m_rhoi)/2;
        m_gammaf = 0.5 + m_alpham - m_alphaf;
    }
    else { feLogError("rhoi must be -1 or between 0 and 1."); return false; }
    
    // allocate vectors
    int neq = m_neq;
    m_Fn.assign(neq, 0);
    m_Fr.assign(neq, 0);
    m_Ui.assign(neq, 0);
    m_Ut.assign(neq, 0);
    m_vi.assign(m_nveq,0);
    m_Vi.assign(m_nveq,0);
    m_di.assign(m_ndeq,0);
    m_Di.assign(m_ndeq,0);
    
    // get number of DOFS
    FEModel& fem = *GetFEModel();
    DOFS& fedofs = fem.GetDOFS();
    int MAX_CDOFS = fedofs.GetVariableSize(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION));
    // allocate concentration-vectors
    m_ci.assign(MAX_CDOFS,vector<double>(0,0));
    m_Ci.assign(MAX_CDOFS,vector<double>(0,0));
    for (int i=0; i<MAX_CDOFS; ++i) {
        m_ci[i].assign(m_nceq[i], 0);
        m_Ci[i].assign(m_nceq[i], 0);
    }
    vector<int> dofs;
    for (int j=0; j<(int)m_nceq.size(); ++j) {
        if (m_nceq[j]) {
            dofs.push_back(m_dofC + j);
        }
    }

    // we need to fill the total DOF vector m_Ut
    // TODO: I need to find an easier way to do this
    FEMesh& mesh = fem.GetMesh();
    gather(m_Ut, mesh, m_dofW[0]);
    gather(m_Ut, mesh, m_dofW[1]);
    gather(m_Ut, mesh, m_dofW[2]);
    gather(m_Ut, mesh, m_dofEF);
    gather(m_Ut, mesh, dofs);

    // set flag for transient or steady-state analyses
    for (int i = 0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        if (fem.GetCurrentStep()->m_nanalysis == FE_STEADY_STATE)
            dom.SetSteadyStateAnalysis();
        else
            dom.SetTransientAnalysis();
    }
    
    return true;
}

//-----------------------------------------------------------------------------
//! Initialize equations
bool FEFluidSolutesSolver::InitEquations()
{
    // base class initialization
    FENewtonSolver::InitEquations();
    
    // determined the nr of velocity and dilatation equations
    FEModel& fem = *GetFEModel();
    m_dofC = fem.GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION), 0);
    m_dofAC = fem.GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION_TDERIV), 0);
    FEMesh& mesh = fem.GetMesh();
    m_nveq = m_ndeq = 0;
    
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& n = mesh.Node(i);
        if (n.m_ID[m_dofW[0]] != -1) m_nveq++;
        if (n.m_ID[m_dofW[1]] != -1) m_nveq++;
        if (n.m_ID[m_dofW[2]] != -1) m_nveq++;
        if (n.m_ID[m_dofEF ] != -1) m_ndeq++;
    }
    
    // determine the nr of concentration equations
    DOFS& fedofs = fem.GetDOFS();
    int MAX_CDOFS = fedofs.GetVariableSize(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION));
    m_nceq.assign(MAX_CDOFS, 0);
    
    // get number of DOFS
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& n = mesh.Node(i);
        for (int j=0; j<MAX_CDOFS; ++j) {
            if (n.m_ID[m_dofC+j] != -1) m_nceq[j]++;
        }
    }

	// get the total concentration equations
	m_nseq = 0;
	for (int i = 0; i < MAX_CDOFS; ++i) m_nseq += m_nceq[i];

	// check that we are using a block scheme for sequential solves
	if ((m_solve_strategy == SOLVE_SEQUENTIAL) && (m_eq_scheme != EQUATION_SCHEME::BLOCK))
	{
		feLogWarning("You need a block solver when using the sequential solve strategy.");
		return false;
	}

	if (m_eq_scheme == EQUATION_SCHEME::BLOCK)
	{
		// repartition the equations so that we only have two partitions,
		// one for the fluid-dilataion, and one for the solutes. 

		// fluid equations is all the rest
		int nfeq = m_neq - m_nseq;

		// create the new partitions
		// Note that this assumes that the solute equations are always last!
		vector<int> p = { nfeq, m_nseq };
		SetPartitions(p);
	}
    
    return true;
}

//-----------------------------------------------------------------------------
void FEFluidSolutesSolver::GetVelocityData(vector<double> &vi, vector<double> &ui)
{
    FEModel& fem = *GetFEModel();
    int N = fem.GetMesh().Nodes(), nid, m = 0;
    zero(vi);
    for (int i=0; i<N; ++i)
    {
        FENode& n = fem.GetMesh().Node(i);
        nid = n.m_ID[m_dofW[0]];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            vi[m++] = ui[nid];
            assert(m <= (int) vi.size());
        }
        nid = n.m_ID[m_dofW[1]];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            vi[m++] = ui[nid];
            assert(m <= (int) vi.size());
        }
        nid = n.m_ID[m_dofW[2]];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            vi[m++] = ui[nid];
            assert(m <= (int) vi.size());
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolutesSolver::GetDilatationData(vector<double> &ei, vector<double> &ui)
{
    FEModel& fem = *GetFEModel();
    int N = fem.GetMesh().Nodes(), nid, m = 0;
    zero(ei);
    for (int i=0; i<N; ++i)
    {
        FENode& n = fem.GetMesh().Node(i);
        nid = n.m_ID[m_dofEF];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            ei[m++] = ui[nid];
            assert(m <= (int) ei.size());
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolutesSolver::GetConcentrationData(vector<double> &ci, vector<double> &ui, const int sol)
{
    FEModel& fem = *GetFEModel();
    int N = fem.GetMesh().Nodes(), nid, m = 0;
    zero(ci);
    for (int i=0; i<N; ++i)
    {
        FENode& n = fem.GetMesh().Node(i);
        nid = n.m_ID[m_dofC+sol];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            ci[m++] = ui[nid];
            assert(m <= (int) ci.size());
        }
    }
}

//-----------------------------------------------------------------------------
//! Update the kinematics of the model, such as nodal positions, velocities,
//! accelerations, etc.
void FEFluidSolutesSolver::UpdateKinematics(vector<double>& ui)
{
    // get the mesh
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    
    // update nodes
    vector<double> U(m_Ut.size());
    for (size_t i=0; i<m_Ut.size(); ++i) U[i] = ui[i] + m_Ui[i] + m_Ut[i];
    
    scatter(U, mesh, m_dofW[0]);
    scatter(U, mesh, m_dofW[1]);
    scatter(U, mesh, m_dofW[2]);
    scatter(U, mesh, m_dofEF);
    
    // get number of DOFS
    DOFS& fedofs = fem.GetDOFS();
    int MAX_CDOFS = fedofs.GetVariableSize(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION));
    
    // update solute data
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        
        // update nodal concentration
        for (int j=0; j<MAX_CDOFS; ++j) {
            int n = node.m_ID[m_dofC+j];
            // Force the concentrations to remain positive
            if (n >= 0) {
                double ct = 0 + m_Ut[n] + m_Ui[n] + ui[n];
                if ((ct < 0.0) && m_forcePositive) ct = 0.0;
                node.set(m_dofC + j, ct);
            }
        }
    }
    
    // force dilatations to remain greater than -1
    if (m_minJf > 0) {
        const int NN = mesh.Nodes();
        for (int i=0; i<NN; ++i)
        {
            FENode& node = mesh.Node(i);
            if (node.get(m_dofEF) <= -1.0)
                node.set(m_dofEF, m_minJf - 1.0);
        }
    }
    
    // make sure the prescribed velocities are fullfilled
    int nvel = fem.BoundaryConditions();
    for (int i=0; i<nvel; ++i)
    {
        FEBoundaryCondition& bc = *fem.BoundaryCondition(i);
        if (bc.IsActive()) bc.Update();
    }
    
    // prescribe DOFs for specialized surface loads
    int nsl = fem.SurfaceLoads();
    for (int i=0; i<nsl; ++i)
    {
        FESurfaceLoad& psl = *fem.SurfaceLoad(i);
        if (psl.IsActive()) psl.Update();
    }
    
    // enforce the linear constraints
    // TODO: do we really have to do this? Shouldn't the algorithm
    // already guarantee that the linear constraints are satisfied?
    FELinearConstraintManager& LCM = fem.GetLinearConstraintManager();
    if (LCM.LinearConstraints() > 0)
    {
        LCM.Update();
    }
    
    // update time derivatives of velocity and dilatation
    // for dynamic simulations
    FEAnalysis* pstep = fem.GetCurrentStep();
    if (pstep->m_nanalysis == FE_DYNAMIC)
    {
        int N = mesh.Nodes();
        double dt = fem.GetTime().timeIncrement;
        double cgi = 1 - 1.0/m_gammaf;
        for (int i=0; i<N; ++i)
        {
            FENode& n = mesh.Node(i);
            
            // velocity time derivative
            vec3d vft = n.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2]);
            vec3d vfp = n.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2]);
            vec3d aft = n.get_vec3d(m_dofAW[0], m_dofAW[1], m_dofAW[2]);
            vec3d afp = n.get_vec3d_prev(m_dofAW[0], m_dofAW[1], m_dofAW[2]);
            aft = afp*cgi + (vft - vfp)/(m_gammaf*dt);
            n.set_vec3d(m_dofAW[0], m_dofAW[1], m_dofAW[2], aft);
            
            // dilatation time derivative
            double eft = n.get(m_dofEF);
            double efp = n.get_prev(m_dofEF);
            double aefp = n.get_prev(m_dofAEF);
            double aeft = aefp*cgi + (eft - efp)/(m_gammaf*dt);
            n.set(m_dofAEF, aeft);
            
            // concentration time derivative
            // update nodal concentration
            for (int j=0; j<MAX_CDOFS; ++j) {
                int k = n.m_ID[m_dofC+j];
                // Force the concentrations to remain positive
                if (k >= 0) {
                    double ct = n.get(m_dofC+j);
                    double cp = n.get_prev(m_dofC+j);
                    double acp = n.get_prev(m_dofAC);
                    double act = acp*cgi + (ct - cp)/(m_gammaf*dt);
                    n.set(m_dofC + j, ct);
                    n.set(m_dofAC + j, act);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolutesSolver::Update2(const vector<double>& ui)
{
    // get the mesh
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    
    // update nodes
    vector<double> U(m_Ut.size());
    for (size_t i = 0; i<m_Ut.size(); ++i) U[i] = ui[i] + m_Ui[i] + m_Ut[i];
    
    scatter(U, mesh, m_dofW[0]);
    scatter(U, mesh, m_dofW[1]);
    scatter(U, mesh, m_dofW[2]);
    scatter(U, mesh, m_dofEF);
    
    // get number of DOFS
    DOFS& fedofs = fem.GetDOFS();
    int MAX_CDOFS = fedofs.GetVariableSize(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION));
    
    // update solute data
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        
        // update nodal concentration
        for (int j=0; j<MAX_CDOFS; ++j) {
            int n = node.m_ID[m_dofC+j];
            // Force the concentrations to remain positive
            if (n >= 0) {
                double ct = 0 + m_Ut[n] + m_Ui[n] + ui[n];
                if ((ct < 0.0) && m_forcePositive) ct = 0.0;
                node.set(m_dofC + j, ct);
            }
        }
    }
    
    // Update the prescribed nodes
    for (int i = 0; i<mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        if (node.m_rid == -1)
        {
            vec3d dv(0, 0, 0);
            for (int j = 0; j < node.m_ID.size(); ++j)
            {
                int nj = -node.m_ID[j] - 2; if (nj >= 0) node.set(j, node.get(j) + ui[nj]);
            }
        }
    }
    
    // update model state
    GetFEModel()->Update();
}

//-----------------------------------------------------------------------------
//! Updates the current state of the model
void FEFluidSolutesSolver::Update(vector<double>& ui)
{
    // update kinematics
    UpdateKinematics(ui);
    
    // update model state
    GetFEModel()->Update();
}

//-----------------------------------------------------------------------------
bool FEFluidSolutesSolver::InitStep(double time)
{
    FEModel& fem = *GetFEModel();
    
    // set time integration parameters
    FETimeInfo& tp = fem.GetTime();
    tp.alphaf = m_alphaf;
    tp.alpham = m_alpham;
    tp.gamma = m_gammaf;
    
    // evaluate load curve values at current (or intermediate) time
    double t = tp.currentTime;
    double dt = tp.timeIncrement;
    double ta = (t > 0) ? t - (1-m_alphaf)*dt : m_alphaf*dt;
    
    return FESolver::InitStep(ta);
}

//-----------------------------------------------------------------------------
//! Prepares the data for the first BFGS-iteration.
void FEFluidSolutesSolver::PrepStep()
{
    FEModel& fem = *GetFEModel();
    // get number of DOFS
    DOFS& fedofs = fem.GetDOFS();
    int MAX_CDOFS = fedofs.GetVariableSize(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION));

    const FETimeInfo& tp = fem.GetTime();
    double dt = tp.timeIncrement;
    
    // zero total DOFs
    zero(m_Ui);
    zero(m_Vi);
    zero(m_Di);
    for (int j=0; j<(int)m_nceq.size(); ++j) if (m_nceq[j]) zero(m_Ci[j]);
    
    // store previous mesh state
    // we need them for strain and acceleration calculations
    FEMesh& mesh = fem.GetMesh();
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& ni = mesh.Node(i);
        ni.m_rp = ni.m_rt = ni.m_r0;
        ni.UpdateValues();
        
        switch (m_pred) {
            case 0:
            {
                // initial guess at start of new time step (default)
                vec3d afp = ni.get_vec3d_prev(m_dofAW[0], m_dofAW[1], m_dofAW[2]);
                ni.set_vec3d(m_dofAW[0], m_dofAW[1], m_dofAW[2], afp*(m_gammaf-1)/m_gammaf);
                ni.set(m_dofAEF, ni.get_prev(m_dofAEF)*(m_gammaf-1)/m_gammaf);
                // update nodal concentration
                for (int j=0; j<MAX_CDOFS; ++j)
                    ni.set(m_dofAC+j, ni.get_prev(m_dofAC+j)*(m_gammaf-1)/m_gammaf);
            }
                break;
                
            case 1:
            {
                // initial guess at start of new time step (Zero Ydot)
                ni.set_vec3d(m_dofAW[0], m_dofAW[1], m_dofAW[2], vec3d(0,0,0));
                ni.set(m_dofAEF, 0);
                for (int j=0; j<MAX_CDOFS; ++j)
                    ni.set(m_dofAC+j, 0);

                vec3d vfp = ni.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2]);
                vec3d afp = ni.get_vec3d_prev(m_dofAW[0], m_dofAW[1], m_dofAW[2]);
                ni.set_vec3d(m_dofW[0], m_dofW[1], m_dofW[2], vfp + afp*dt*(1-m_gammaf)*m_alphaf);
                ni.set(m_dofEF, ni.get_prev(m_dofEF) + ni.get_prev(m_dofAEF)*dt*(1-m_gammaf)*m_alphaf);
                for (int j=0; j<MAX_CDOFS; ++j)
                    ni.set(m_dofC+j, ni.get_prev(m_dofC+j) + ni.get_prev(m_dofAC+j)*dt*(1-m_gammaf)*m_alphaf);
            }
                break;
                
            case 2:
            {
                // initial guess at start of new time step (Same Ydot)
                vec3d afp = ni.get_vec3d_prev(m_dofAW[0], m_dofAW[1], m_dofAW[2]);
                ni.set_vec3d(m_dofAW[0], m_dofAW[1], m_dofAW[2], afp);
                ni.set(m_dofAEF, ni.get_prev(m_dofAEF));
                for (int j=0; j<MAX_CDOFS; ++j)
                    ni.set(m_dofAC+j, ni.get_prev(m_dofAC+j));

                vec3d vfp = ni.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2]);
                ni.set_vec3d(m_dofW[0], m_dofW[1], m_dofW[2], vfp + afp*dt);
                ni.set(m_dofEF, ni.get_prev(m_dofEF) + ni.get_prev(m_dofAEF)*dt);
                for (int j=0; j<MAX_CDOFS; ++j)
                    ni.set(m_dofC+j, ni.get_prev(m_dofC+j) + ni.get_prev(m_dofAC+j)*dt);
            }
                break;
                
            default:
                break;
        }
    }
    
    // apply concentrated nodal forces
    // since these forces do not depend on the geometry
    // we can do this once outside the NR loop.
    vector<double> dummy(m_neq, 0.0);
    zero(m_Fn);
    FEGlobalVector Fn(*GetFEModel(), m_Fn, dummy);
    NodalLoads(Fn, tp);
    
    // apply prescribed velocities
    // we save the prescribed velocity increments in the ui vector
    vector<double>& ui = m_ui;
    zero(ui);
    int nbc = fem.BoundaryConditions();
    for (int i=0; i<nbc; ++i)
    {
        FEBoundaryCondition& bc = *fem.BoundaryCondition(i);
        if (bc.IsActive()) bc.PrepStep(ui);
    }
    
    // apply prescribed DOFs for specialized surface loads
    int nsl = fem.SurfaceLoads();
    for (int i=0; i<nsl; ++i)
    {
        FESurfaceLoad& psl = *fem.SurfaceLoad(i);
        if (psl.IsActive()) psl.Update();
    }
    
    // intialize material point data
    // NOTE: do this before the stresses are updated
    // TODO: does it matter if the stresses are updated before
    //       the material point data is initialized
    // update domain data
    for (int i=0; i<mesh.Domains(); ++i) mesh.Domain(i).PreSolveUpdate(tp);
    
    // update stresses
    fem.Update();
    
    // see if we need to do contact augmentations
    m_baugment = false;
    for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
    {
        FEContactInterface& ci = dynamic_cast<FEContactInterface&>(*fem.SurfacePairConstraint(i));
        if (ci.IsActive() && ci.m_blaugon) m_baugment = true;
    }
    
    // see if we have to do nonlinear constraint augmentations
    if (fem.NonlinearConstraints() != 0) m_baugment = true;
}

//-----------------------------------------------------------------------------
bool FEFluidSolutesSolver::Quasin()
{
    FEModel& fem = *GetFEModel();
    // convergence norms
    double    normR1;        // residual norm
    double    normE1;        // energy norm
    double    normV;        // velocity norm
    double    normv;        // velocity increment norm
    double    normRi = 0;    // initial residual norm
    double    normVi = 0;    // initial velocity norm
    double    normEi = 0; // initial energy norm
    double    normEm = 0;    // max energy norm
    double    normDi = 0;    // initial dilatation norm
    double    normD;        // current dilatation norm
    double    normd;        // incremement dilatation norm
    
    // get number of DOFS
    DOFS& fedofs = fem.GetDOFS();
    int MAX_CDOFS = fedofs.GetVariableSize(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION));
    
    // solute convergence data
    vector<double>    normCi(MAX_CDOFS);    // initial concentration norm
    vector<double>    normC(MAX_CDOFS);    // current concentration norm
    vector<double>    normc(MAX_CDOFS);    // incremement concentration norm
    
    // Get the current step
    FEAnalysis* pstep = fem.GetCurrentStep();
    
    // prepare for the first iteration
    const FETimeInfo& tp = fem.GetTime();
    PrepStep();
    
    // Init QN method
    if (QNInit() == false) return false;

	// this flag indicates whether the velocity has converged for a sequential solve
	// (This is not used for a coupled solve.)
	bool vel_converged = false;
    
    // loop until converged or when max nr of reformations reached
    bool bconv = false; // convergence flag
    do
    {
        feLog(" %d\n", m_niter+1);
        
        // assume we'll converge.
        bconv = true;

		// for sequential solve, we set one of the residual components to zero
		if (m_solve_strategy == SOLVE_SEQUENTIAL)
		{
			int veq = m_neq - m_nseq;
			if (vel_converged == false)
			{
				// zero the solute residual
				for (int i = veq; i < m_neq; ++i) m_R0[i] = 0.0;
			}
			else
			{
				// zero the velocity residual
				for (int i = 0; i < veq; ++i) m_R0[i] = 0.0;
			}
		}
        
        // solve the equations (returns line search; solution stored in m_ui)
        double s = QNSolve();

		// for sequential solve, we set one of the residual components to zero
		if (m_solve_strategy == SOLVE_SEQUENTIAL)
		{
			int veq = m_neq - m_nseq;
			if (vel_converged == false)
			{
				// zero the solute residual
				for (int i = veq; i < m_neq; ++i) m_R1[i] = 0.0;

				// zero the solute solution
				for (int i = veq; i < m_neq; ++i) m_ui[i] = 0.0;
			}
			else
			{
				// zero the velocity residual
				for (int i = 0; i < veq; ++i) m_R1[i] = 0.0;
			}
		}

        // extract the velocity and dilatation increments
        GetVelocityData(m_vi, m_ui);
        GetDilatationData(m_di, m_ui);
        
        // set initial convergence norms
        if (m_niter == 0)
        {
            normRi = fabs(m_R0*m_R0);
            normEi = fabs(m_ui*m_R0);
            normVi = fabs(m_vi*m_vi);
            normDi = fabs(m_di*m_di);
            normEm = normEi;
        }
        
        // calculate norms
        // update all degrees of freedom
        for (int i=0; i<m_neq; ++i) m_Ui[i] += s*m_ui[i];
        
        // update velocities
        for (int i = 0; i<m_nveq; ++i) m_Vi[i] += s*m_vi[i];
        
        // update dilatations
        for (int i = 0; i<m_ndeq; ++i) m_Di[i] += s*m_di[i];
        
        // calculate the norms
        normR1 = m_R1*m_R1;
        normv  = (m_vi*m_vi)*(s*s);
        normV  = m_Vi*m_Vi;
        normd  = (m_di*m_di)*(s*s);
        normD  = m_Di*m_Di;
        normE1 = s*fabs(m_ui*m_R1);
        
        // check for nans
        if (ISNAN(normR1)) throw NANDetected();
        
        // check residual norm
        if ((m_Rtol > 0) && (normR1 > m_Rtol*normRi)) bconv = false;
        
        // check velocity norm
        if ((m_Vtol > 0) && (normv  > (m_Vtol*m_Vtol)*normV )) bconv = false;
        
        // check dilatation norm
        if ((m_Ftol > 0) && (normd  > (m_Ftol*m_Ftol)*normD )) bconv = false;
        
        // check energy norm
        if ((m_Etol > 0) && (normE1 > m_Etol*normEi)) bconv = false;
        
        // check linestep size
        if ((m_lineSearch->m_LStol > 0) && (s < m_lineSearch->m_LSmin)) bconv = false;
        
        // check energy divergence
        if (normE1 > normEm) bconv = false;
        
        // check solute convergence
        // extract the concentration increments
        for (int j = 0; j<(int)m_nceq.size(); ++j) {
            if (m_nceq[j]) {
                GetConcentrationData(m_ci[j], m_ui,j);
                
                // set initial norm
                if (m_niter == 0)
                    normCi[j] = fabs(m_ci[j]*m_ci[j]);
                
                // update total concentration
                for (int i = 0; i<m_nceq[j]; ++i) m_Ci[j][i] += s*m_ci[j][i];
                
                // calculate norms
                normC[j] = m_Ci[j]*m_Ci[j];
                normc[j] = (m_ci[j]*m_ci[j])*(s*s);
                
            }
        }
        
        // check convergence
        if (m_Ctol > 0) {
            for (int j = 0; j<(int)m_nceq.size(); ++j)
                if (m_nceq[j]) bconv = bconv && (normc[j] <= (m_Ctol*m_Ctol)*normC[j]);
        }
        
        // print convergence summary
        feLog(" Nonlinear solution status: time= %lg\n", tp.currentTime);
        feLog("\tstiffness updates             = %d\n", m_qnstrategy->m_nups);
        feLog("\tright hand side evaluations   = %d\n", m_nrhs);
        feLog("\tstiffness matrix reformations = %d\n", m_nref);
        if (m_lineSearch->m_LStol > 0) feLog("\tstep from line search         = %lf\n", s);
        feLog("\tconvergence norms :     INITIAL         CURRENT         REQUIRED\n");
        feLog("\t   residual         %15le %15le %15le \n", normRi, normR1, m_Rtol*normRi);
        feLog("\t   energy           %15le %15le %15le \n", normEi, normE1, m_Etol*normEi);
        feLog("\t   velocity         %15le %15le %15le \n", normVi, normv ,(m_Vtol*m_Vtol)*normV );
        feLog("\t   dilatation       %15le %15le %15le \n", normDi, normd ,(m_Ftol*m_Ftol)*normD );
        for (int j = 0; j<(int)m_nceq.size(); ++j) {
            if (m_nceq[j])
                feLog("\t solute %d concentration %15le %15le %15le\n", j+1, normCi[j], normc[j] ,(m_Ctol*m_Ctol)*normC[j] );
        }

        // see if we may have a small residual
        if ((bconv == false) && (normR1 < m_Rmin))
        {
            // check for almost zero-residual on the first iteration
            // this might be an indication that there is no force on the system
            feLogWarning("No force acting on the system.");
            bconv = true;
        }
        
        // see if we have exceeded the max residual
        if ((bconv == false) && (m_Rmax > 0) && (normR1 >= m_Rmax))
        {
            // doesn't look like we're getting anywhere, so let's retry the time step
            throw MaxResidualError();
        }
        
        // check if we have converged.
        // If not, calculate the BFGS update vectors
        if (bconv == false)
        {
            if (s < m_lineSearch->m_LSmin)
            {
                // check for zero linestep size
                feLogWarning("Zero linestep size. Stiffness matrix will now be reformed");
                QNForceReform(true);
            }
            else if ((normE1 > normEm) && m_bdivreform)
            {
                // check for diverging
                feLogWarning("Problem is diverging. Stiffness matrix will now be reformed");
                normEm = normE1;
                normEi = normE1;
                normRi = normR1;
                normVi = normv;
                normDi = normd;
                for (int j = 0; j<(int)m_nceq.size(); ++j)
                    if (m_nceq[j]) normCi[j] = normc[j];
                QNForceReform(true);
            }
            
            // Do the QN update (This may also do a stiffness reformation if necessary)
            bool bret = QNUpdate();
            
            // something went wrong with the update, so we'll need to break
            if (bret == false) break;
        }
        else if (m_baugment)
        {
            // Do augmentations
            bconv = DoAugmentations();
        }

		if (bconv && (m_solve_strategy == SOLVE_SEQUENTIAL))
		{
			if (vel_converged == false)
			{
				vel_converged = true;
				bconv = false;
				m_qnstrategy->m_nups = 0;
				m_niter = -1;
				Residual(m_R0);
				feLog("\n*** Velocity converged. Now solving for solutes.\n");
			}
		}
        
        // increase iteration number
        m_niter++;
        
        // do minor iterations callbacks
        fem.DoCallback(CB_MINOR_ITERS);
    }
    while (bconv == false);
    
    // if converged we update the total velocities
    if (bconv)
    {
        m_Ut += m_Ui;
        zero(m_Ui);
    }
    
    return bconv;
}

//-----------------------------------------------------------------------------
//! Calculates global stiffness matrix.

bool FEFluidSolutesSolver::StiffnessMatrix(FELinearSystem& LS)
{
    FEModel& fem = *GetFEModel();
    
    const FETimeInfo& tp = fem.GetTime();
    
    // get the mesh
    FEMesh& mesh = fem.GetMesh();
    
    // calculate the stiffness matrix for each domain
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        dom.StiffnessMatrix(LS, tp);
    }
    
    // calculate the body force stiffness matrix for each domain
    int NBL = fem.BodyLoads();
    for (int j = 0; j<NBL; ++j)
    {
        FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.GetBodyLoad(j));
        if (pbf && pbf->IsActive())
        {
            for (int i = 0; i<pbf->Domains(); ++i)
            {
                FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(*pbf->Domain(i));
                dom.BodyForceStiffness(LS, tp, *pbf);
            }
        }
    }
    
    // calculate contact stiffness
    ContactStiffness(LS);
    
    // calculate stiffness matrix due to surface loads
    int nsl = fem.SurfaceLoads();
    for (int i=0; i<nsl; ++i)
    {
        FESurfaceLoad* psl = fem.SurfaceLoad(i);
        if (psl->IsActive()) psl->StiffnessMatrix(LS, tp);
    }
    
    // Add mass matrix
    // loop over all domains
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        dom.MassMatrix(LS, tp);
    }
    
    // calculate nonlinear constraint stiffness
    // note that this is the contribution of the
    // constrainst enforced with augmented lagrangian
    NonLinearConstraintStiffness(LS, tp);
    
    return true;
}

//-----------------------------------------------------------------------------
//! Calculate the stiffness contribution due to nonlinear constraints
void FEFluidSolutesSolver::NonLinearConstraintStiffness(FELinearSystem& LS, const FETimeInfo& tp)
{
    FEModel& fem = *GetFEModel();
    
    int N = fem.NonlinearConstraints();
    for (int i=0; i<N; ++i)
    {
        FENLConstraint* plc = fem.NonlinearConstraint(i);
        if (plc->IsActive()) plc->StiffnessMatrix(LS, tp);
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the contact stiffness matrix

void FEFluidSolutesSolver::ContactStiffness(FELinearSystem& LS)
{
    FEModel& fem = *GetFEModel();
    
    const FETimeInfo& tp = fem.GetTime();
    for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
    {
        FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));
        if (pci->IsActive()) pci->StiffnessMatrix(LS, tp);
    }
}

//-----------------------------------------------------------------------------
//! Calculates the contact forces
void FEFluidSolutesSolver::ContactForces(FEGlobalVector& R)
{
    FEModel& fem = *GetFEModel();
    
    const FETimeInfo& tp = fem.GetTime();
    for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
    {
        FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));
        if (pci->IsActive()) pci->Residual(R, tp);
    }
}

//-----------------------------------------------------------------------------
//! calculates the residual vector
//! Note that the concentrated nodal forces are not calculated here.
//! This is because they do not depend on the geometry
//! so we only calculate them once (in Quasin) and then add them here.

bool FEFluidSolutesSolver::Residual(vector<double>& R)
{
    FEModel& fem = *GetFEModel();
    
    // get the time information
    const FETimeInfo& tp = fem.GetTime();
    
    // initialize residual with concentrated nodal loads
    R = m_Fn;
    
    // zero nodal reaction forces
    zero(m_Fr);
    
    // setup the global vector
    FEFluidResidualVector RHS(fem, R, m_Fr);
    
    // get the mesh
    FEMesh& mesh = fem.GetMesh();
    
    // calculate the internal (stress) forces
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        dom.InternalForces(RHS, tp);
    }
    
    // calculate the body forces
    for (int j = 0; j<fem.BodyLoads(); ++j)
    {
        FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.GetBodyLoad(j));
        if (pbf && pbf->IsActive())
        {
            for (int i = 0; i<pbf->Domains(); ++i)
            {
                FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(*pbf->Domain(i));
                dom.BodyForce(RHS, tp, *pbf);
            }
        }
    }
    
    // calculate inertial forces
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        dom.InertialForces(RHS, tp);
    }
    
    // calculate forces due to surface loads
    int nsl = fem.SurfaceLoads();
    for (int i=0; i<nsl; ++i)
    {
        FESurfaceLoad* psl = fem.SurfaceLoad(i);
        if (psl->IsActive()) psl->Residual(RHS, tp);
    }
    
    // calculate contact forces
    ContactForces(RHS);
    
    // calculate nonlinear constraint forces
    // note that these are the linear constraints
    // enforced using the augmented lagrangian
    NonLinearConstraintForces(RHS, tp);
    
    // add model loads
    int NML = fem.ModelLoads();
    for (int i=0; i<NML; ++i)
    {
        FEModelLoad& mli = *fem.ModelLoad(i);
        if (mli.IsActive())
        {
            mli.Residual(RHS, tp);
        }
    }
    
    // set the nodal reaction forces
    // TODO: Is this a good place to do this?
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        node.m_Fr = vec3d(0,0,0);
        
        int n;
        if ((n = -node.m_ID[m_dofW[0]]-2) >= 0) node.m_Fr.x = -m_Fr[n];
        if ((n = -node.m_ID[m_dofW[1]]-2) >= 0) node.m_Fr.y = -m_Fr[n];
        if ((n = -node.m_ID[m_dofW[2]]-2) >= 0) node.m_Fr.z = -m_Fr[n];
    }
    
    // increase RHS counter
    m_nrhs++;
    
    return true;
}

//-----------------------------------------------------------------------------
//! calculate the nonlinear constraint forces
void FEFluidSolutesSolver::NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp)
{
    FEModel& fem = *GetFEModel();
    
    int N = fem.NonlinearConstraints();
    for (int i=0; i<N; ++i)
    {
        FENLConstraint* plc = fem.NonlinearConstraint(i);
        if (plc->IsActive()) plc->Residual(R, tp);
    }
}

//-----------------------------------------------------------------------------
//! Serialization
void FEFluidSolutesSolver::Serialize(DumpStream& ar)
{
    FENewtonSolver::Serialize(ar);
    if (ar.IsShallow()) return;
    ar & m_nveq & m_ndeq;
    ar & m_alphaf & m_alpham;
    ar & m_gammaf;
    ar & m_pred;
    
    ar & m_Fn & m_Fr & m_Ui &m_Ut;
    ar & m_Vi & m_vi;
    ar & m_Di & m_di;

    ar & m_ci & m_Ci;
}
