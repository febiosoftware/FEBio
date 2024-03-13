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



#include "stdafx.h"
#include "FEBioThermoFluid.h"
#include "FEFluidHeatSupply.h"
#include "FEFluidResistanceBC.h"
#include "FEBackFlowStabilization.h"
#include "FEFluidNormalVelocity.h"
#include "FEFluidVelocity.h"
#include "FEFluidRotationalVelocity.h"
#include "FETiedFluidInterface.h"
#include "FEThermoFluidSolver.h"
#include "FEThermoFluidDomain3D.h"
#include "FEFluidDomain.h"
#include <assert.h>
#include "FEFluidResidualVector.h"
#include <FEBioMech/FEResidualVector.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <FECore/DOFS.h>
#include <FECore/FEGlobalMatrix.h>
#include <FECore/sys.h>
#include <FEBioMech/FEBodyForce.h>
#include <FECore/FEBoundaryCondition.h>
#include <FECore/FENodalLoad.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEModelLoad.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FELinearConstraintManager.h>
#include <FECore/FELinearSystem.h>
#include <FECore/FENLConstraint.h>
#include "FEThermoFluidAnalysis.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_FECORE_CLASS(FEThermoFluidSolver, FENewtonSolver)
    ADD_PARAMETER(m_Vtol , "vtol"        );
    ADD_PARAMETER(m_Ftol , "ftol"        );
    ADD_PARAMETER(m_Ttol , "ttol"        );
    ADD_PARAMETER(m_Etol, FE_RANGE_GREATER_OR_EQUAL(0.0), "etol");
    ADD_PARAMETER(m_Rtol, FE_RANGE_GREATER_OR_EQUAL(0.0), "rtol");
    ADD_PARAMETER(m_rhoi , "rhoi"        );
    ADD_PARAMETER(m_pred , "predictor"   );
    ADD_PARAMETER(m_minJf, "min_volume_ratio");
    ADD_PARAMETER(m_minT , "min_abs_temperature");
    ADD_PARAMETER(m_solve_strategy, "solve_strategy")->setEnums("coupled\0sequential\0");
    ADD_PARAMETER(m_Tmin , "min_T_drop");
    ADD_PARAMETER(m_Tmax , "min_T_rise");
    ADD_PARAMETER(m_Tnum , "min_T_num");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! FEThermoFluidSolver Construction
//
FEThermoFluidSolver::FEThermoFluidSolver(FEModel* pfem) : FENewtonSolver(pfem), m_dofW(pfem), m_dofAW(pfem), m_dofEF(pfem), m_dofT(pfem)
{
    // default values
    m_Rtol = 0.001;
    m_Etol = 0.01;
    m_Vtol = 0.001;
    m_Ftol = 0.001;
    m_Ttol = 0.001;
    m_Rmin = 1.0e-20;
    m_Rmax = 0;     // not used if zero
    m_minJf = 0;    // not used if zero
    m_minT = 0;     // not used if zero
    m_Tmin = 0;     // not used if zero
    m_Tmax = 0;     // not used if zero
    m_Tnum = 1;

    m_nveq = 0;
    m_ndeq = 0;
    m_nteq = 0;
    m_niter = 0;
    
    // assume non-symmetric stiffness
    m_msymm = REAL_UNSYMMETRIC;

    m_solve_strategy = SOLVE_COUPLED;
    
    m_rhoi = 0;
    m_pred = 0;
    
    m_sudden_T_change = false;
    
    // Preferred strategy is Broyden's method
    SetDefaultStrategy(QN_BROYDEN);

    // turn off checking for a zero diagonal
    CheckZeroDiagonal(false);

    // get the dof indices
    // TODO: Can this be done in Init, since  there is no error checking
    if (pfem)
    {
        m_dofW.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::RELATIVE_FLUID_VELOCITY));
        m_dofAW.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::RELATIVE_FLUID_ACCELERATION));
        m_dofEF.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::FLUID_DILATATION));
        m_dofAEF = pfem->GetDOFIndex(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::FLUID_DILATATION_TDERIV), 0);
        m_dofT.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::TEMPERATURE));
        m_dofAT = pfem->GetDOFIndex(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::TEMPERATURE_TDERIV), 0);
    }
}

//-----------------------------------------------------------------------------
FEThermoFluidSolver::~FEThermoFluidSolver()
{

}

//-----------------------------------------------------------------------------
//! Allocates and initializes the data structures used by the FEThermoFluidSolver
//
bool FEThermoFluidSolver::Init()
{
    // initialize base class
    if (FENewtonSolver::Init() == false) return false;

    // check parameters
    if (m_Vtol <  0.0) { feLogError("vtol must be nonnegative."); return false; }
    if (m_Ftol <  0.0) { feLogError("dtol must be nonnegative."); return false; }
    if (m_Ttol <  0.0) { feLogError("ttol must be nonnegative."); return false; }
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
    m_Fr.assign(neq, 0);
    m_Ui.assign(neq, 0);
    m_Ut.assign(neq, 0);
    m_vi.assign(m_nveq,0);
    m_Vi.assign(m_nveq,0);
    m_di.assign(m_ndeq,0);
    m_Di.assign(m_ndeq,0);
    m_ti.assign(m_nteq,0);
    m_Ti.assign(m_nteq,0);

    // we need to fill the total DOF vector m_Ut
    // TODO: I need to find an easier way to do this
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    gather(m_Ut, mesh, m_dofW[0]);
    gather(m_Ut, mesh, m_dofW[1]);
    gather(m_Ut, mesh, m_dofW[2]);
    gather(m_Ut, mesh, m_dofEF[0]);
    gather(m_Ut, mesh, m_dofT[0]);

    // set flag for transient or steady-state analyses
    FEAnalysis* pstep = fem.GetCurrentStep();
    for (int i = 0; i<mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        if (dom.IsActive()) {
            FEFluidDomain* fdom = dynamic_cast<FEFluidDomain*>(&dom);
            if (pstep->m_nanalysis == FEThermoFluidAnalysis::STEADY_STATE)
                fdom->SetSteadyStateAnalysis();
            else
                fdom->SetTransientAnalysis();
        }
    }

    return true;
}

//-----------------------------------------------------------------------------
//! Initialize equations
bool FEThermoFluidSolver::InitEquations()
{
    // Add the solution variables
    AddSolutionVariable(&m_dofW , 1, "velocity"   , m_Vtol);
    AddSolutionVariable(&m_dofEF, 1, "dilatation" , m_Ftol);
    AddSolutionVariable(&m_dofT , 1, "temperature", m_Ttol);

    // base class initialization
    if (FENewtonSolver::InitEquations() == false) return false;
    
    // determined the nr of velocity and dilatation equations
    FEMesh& mesh = GetFEModel()->GetMesh();
    m_nveq = m_ndeq = m_nteq = 0;
    
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& n = mesh.Node(i);
        if (n.m_ID[m_dofW[0]] != -1) m_nveq++;
        if (n.m_ID[m_dofW[1]] != -1) m_nveq++;
        if (n.m_ID[m_dofW[2]] != -1) m_nveq++;
        if (n.m_ID[m_dofEF[0]] != -1) m_ndeq++;
        if (n.m_ID[m_dofT[0]] != -1) m_nteq++;
    }

    // check that we are using a block scheme for sequential solves
    if ((m_solve_strategy == SOLVE_SEQUENTIAL) && (m_eq_scheme != EQUATION_SCHEME::BLOCK))
    {
        feLogWarning("You need a block solver when using the sequential solve strategy.");
        return false;
    }
    
    // Next, we add any Lagrange Multipliers
    FEModel& fem = *GetFEModel();
    for (int i = 0; i < fem.NonlinearConstraints(); ++i)
    {
        FENLConstraint* lmc = fem.NonlinearConstraint(i);
        if (lmc->IsActive())
        {
            m_neq += lmc->InitEquations(m_neq);
        }
    }
    for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
    {
        FESurfacePairConstraint* spc = fem.SurfacePairConstraint(i);
        if (spc->IsActive())
        {
            m_neq += spc->InitEquations(m_neq);
        }
    }
    
    if (m_eq_scheme == EQUATION_SCHEME::BLOCK)
    {
        // repartition the equations so that we only have two partitions,
        // one for the fluid-dilatation, and one for the temperature.
        
        // fluid equations is all the rest
        int nfeq = m_neq - m_nteq;
        
        // create the new partitions
        // Note that this assumes that the temperature equations are always last!
        vector<int> p = { nfeq, m_nteq };
        SetPartitions(p);
    }
    
    return true;
}

//-----------------------------------------------------------------------------
//! Initialize equations
bool FEThermoFluidSolver::InitEquations2()
{
    // Add the solution variables
    AddSolutionVariable(&m_dofW , 1, "velocity"   , m_Vtol);
    AddSolutionVariable(&m_dofEF, 1, "dilatation" , m_Ftol);
    AddSolutionVariable(&m_dofT , 1, "temperature", m_Ttol);

    // base class initialization
    if (FENewtonSolver::InitEquations2() == false) return false;
    
    // determined the nr of velocity and dilatation equations
    FEMesh& mesh = GetFEModel()->GetMesh();
    m_nveq = m_ndeq = m_nteq = 0;
    
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& n = mesh.Node(i);
        if (n.m_ID[m_dofW[0]] != -1) m_nveq++;
        if (n.m_ID[m_dofW[1]] != -1) m_nveq++;
        if (n.m_ID[m_dofW[2]] != -1) m_nveq++;
        if (n.m_ID[m_dofEF[0]] != -1) m_ndeq++;
        if (n.m_ID[m_dofT[0]] != -1) m_nteq++;
    }

    // Next, we add any Lagrange Multipliers
    FEModel& fem = *GetFEModel();
    for (int i = 0; i < fem.NonlinearConstraints(); ++i)
    {
        FENLConstraint* lmc = fem.NonlinearConstraint(i);
        if (lmc->IsActive())
        {
            m_neq += lmc->InitEquations(m_neq);
        }
    }
    for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
    {
        FESurfacePairConstraint* spc = fem.SurfacePairConstraint(i);
        if (spc->IsActive())
        {
            m_neq += spc->InitEquations(m_neq);
        }
    }
    
    return true;
}

//-----------------------------------------------------------------------------
void FEThermoFluidSolver::GetVelocityData(vector<double> &vi, vector<double> &ui)
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
void FEThermoFluidSolver::GetDilatationData(vector<double> &ei, vector<double> &ui)
{
    FEModel& fem = *GetFEModel();
    int N = fem.GetMesh().Nodes(), nid, m = 0;
    zero(ei);
    for (int i=0; i<N; ++i)
    {
        FENode& n = fem.GetMesh().Node(i);
        nid = n.m_ID[m_dofEF[0]];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            ei[m++] = ui[nid];
            assert(m <= (int) ei.size());
        }
    }
}

//-----------------------------------------------------------------------------
void FEThermoFluidSolver::GetTemperatureData(vector<double> &ti, vector<double> &ui)
{
    FEModel& fem = *GetFEModel();
    int N = fem.GetMesh().Nodes(), nid, m = 0;
    zero(ti);
    for (int i=0; i<N; ++i)
    {
        FENode& n = fem.GetMesh().Node(i);
        nid = n.m_ID[m_dofT[0]];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            ti[m++] = ui[nid];
            assert(m <= (int) ti.size());
        }
    }
}

//-----------------------------------------------------------------------------
//! Update the kinematics of the model, such as nodal positions, velocities,
//! accelerations, etc.
void FEThermoFluidSolver::UpdateKinematics(vector<double>& ui)
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
    scatter(U, mesh, m_dofEF[0]);
//    scatter(U, mesh, m_dofT[0]);

    // update temperature data
    int nssd = 0, nssr = 0;
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        
        // update nodal temperature
        int n = node.m_ID[m_dofT[0]];
        // Force the temperature to remain positive
        if (n >= 0) {
            double Tt = 0 + m_Ut[n] + m_Ui[n] + ui[n];
            double Tp = node.get_prev(m_dofT[0]);
            if ((m_Tmin > 0) && (node.get_bc(m_dofT[0]) == DOF_OPEN) && (Tp - Tt >= m_Tmin))
                nssd++;
            if ((m_Tmax > 0) && (node.get_bc(m_dofT[0]) == DOF_OPEN) && (Tt - Tp >= m_Tmax))
                nssr++;
            node.set(m_dofT[0], Tt);
        }
    }
    
    if (nssd >= m_Tnum) m_sudden_T_change = true;
    if (nssr >= m_Tnum) m_sudden_T_change = true;
    
    // force dilatations to remain greater than -1
    if (m_minJf > 0) {
        const int NN = mesh.Nodes();
        for (int i=0; i<NN; ++i)
        {
            FENode& node = mesh.Node(i);
            if (node.get(m_dofEF[0]) <= -1.0)
                node.set(m_dofEF[0], m_minJf - 1.0);
        }
    }

    // force absolute temperature to remain greater than 0
    double Tr = fem.GetGlobalConstant("T");
    if (m_minT > 0) {
        const int NN = mesh.Nodes();
        for (int i=0; i<NN; ++i)
        {
            FENode& node = mesh.Node(i);
            if (node.get(m_dofT[0]) <= -Tr)
                node.set(m_dofT[0], m_minT - Tr);
        }
    }
    
    // make sure the prescribed velocities are fulfilled
    int nvel = fem.BoundaryConditions();
    for (int i=0; i<nvel; ++i)
    {
        FEBoundaryCondition& bc = *fem.BoundaryCondition(i);
        if (bc.IsActive() && HasActiveDofs(bc.GetDofList())) bc.Update();
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
    if (pstep->m_nanalysis == FEThermoFluidAnalysis::DYNAMIC)
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
            double eft = n.get(m_dofEF[0]);
            double efp = n.get_prev(m_dofEF[0]);
            double aefp = n.get_prev(m_dofAEF);
            double aeft = aefp*cgi + (eft - efp)/(m_gammaf*dt);
            n.set(m_dofAEF, aeft);
            
            // temperature time derivative
            double Tt = n.get(m_dofT[0]);
            double Tp = n.get_prev(m_dofT[0]);
            double aTp = n.get_prev(m_dofAT);
            double aTt = aTp*cgi + (Tt - Tp)/(m_gammaf*dt);
            n.set(m_dofAT, aTt);
        }
    }
    // update nonlinear constraints (needed for updating Lagrange Multiplier)
    for (int i = 0; i < fem.NonlinearConstraints(); ++i)
    {
        FENLConstraint* nlc = fem.NonlinearConstraint(i);
        if (nlc->IsActive()) nlc->Update(m_Ui, ui);
    }
    for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
    {
        FESurfacePairConstraint* spc = fem.SurfacePairConstraint(i);
        if (spc->IsActive()) spc->Update(ui);
    }
}

//-----------------------------------------------------------------------------
void FEThermoFluidSolver::Update2(const vector<double>& ui)
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
    scatter(U, mesh, m_dofEF[0]);
    scatter(U, mesh, m_dofT[0]);

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
void FEThermoFluidSolver::Update(vector<double>& ui)
{
    FEModel& fem = *GetFEModel();
    FETimeInfo& tp = fem.GetTime();
    tp.currentIteration = m_niter;
    
    // update kinematics
    UpdateKinematics(ui);
    
    // update model state
//    GetFEModel()->Update();
    UpdateModel();
}

//-----------------------------------------------------------------------------
//! Update DOF increments
void FEThermoFluidSolver::UpdateIncrements(vector<double>& Ui, vector<double>& ui, bool emap)
{
    FEModel& fem = *GetFEModel();
    
    // get the mesh
    FEMesh& mesh = fem.GetMesh();
    
    // extract the velocity and dilatation increments
    GetVelocityData(m_vi, ui);
    GetDilatationData(m_di, ui);
    GetTemperatureData(m_ti, ui);

    // update all degrees of freedom
    for (int i=0; i<m_neq; ++i) Ui[i] += ui[i];
        
    // update velocities
    for (int i = 0; i<m_nveq; ++i) m_Vi[i] += m_vi[i];

    // update dilatations
    for (int i = 0; i<m_ndeq; ++i) m_Di[i] += m_di[i];
        
    // update temperatures
    for (int i = 0; i<m_nteq; ++i) m_Ti[i] += m_ti[i];
        
    for (int i = 0; i < fem.NonlinearConstraints(); ++i)
    {
        FENLConstraint* plc = fem.NonlinearConstraint(i);
        if (plc && plc->IsActive()) plc->UpdateIncrements(Ui, ui);
    }
    
    // TODO: This is a hack!
    // The problem is that I only want to call the domain's IncrementalUpdate during
    // the quasi-Newtoon loop. However, this function is also called after the loop
    // converges. The emap parameter is used here to detect wether we are inside the
    // loop (emap == false), or not (emap == true).
    if (emap == false)
    {
        for (int i = 0; i < mesh.Domains(); ++i)
        {
            FEDomain& dom = mesh.Domain(i);
            dom.IncrementalUpdate(ui, true);
        }
    }
}

//-----------------------------------------------------------------------------
//! Update nonlinear constraints
void FEThermoFluidSolver::UpdateConstraints()
{
    FEModel& fem = *GetFEModel();
    FETimeInfo& tp = fem.GetTime();
    tp.currentIteration = m_niter;
    
    // Update all nonlinear constraints
    for (int i = 0; i<fem.NonlinearConstraints(); ++i)
    {
        FENLConstraint* pci = fem.NonlinearConstraint(i);
        if (pci->IsActive()) pci->Update();
    }
}

//-----------------------------------------------------------------------------
bool FEThermoFluidSolver::InitStep(double time)
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
void FEThermoFluidSolver::PrepStep()
{
    FEModel& fem = *GetFEModel();

    FETimeInfo& tp = fem.GetTime();
    double dt = tp.timeIncrement;
    tp.currentIteration = m_niter;

    // zero total DOFs
    zero(m_Ui);
    zero(m_Vi);
    zero(m_Di);
    zero(m_Ti);

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
                ni.set(m_dofAT, ni.get_prev(m_dofAT)*(m_gammaf-1)/m_gammaf);
            }
                break;
                
            case 1:
            {
                // initial guess at start of new time step (Zero Ydot)
                ni.set_vec3d(m_dofAW[0], m_dofAW[1], m_dofAW[2], vec3d(0,0,0));
                ni.set(m_dofAEF, 0);
                ni.set(m_dofAT, 0);

                vec3d vfp = ni.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2]);
                vec3d afp = ni.get_vec3d_prev(m_dofAW[0], m_dofAW[1], m_dofAW[2]);
                ni.set_vec3d(m_dofW[0], m_dofW[1], m_dofW[2], vfp + afp*dt*(1-m_gammaf)*m_alphaf);
                ni.set(m_dofEF[0], ni.get_prev(m_dofEF[0]) + ni.get_prev(m_dofAEF)*dt*(1-m_gammaf)*m_alphaf);
                ni.set(m_dofT[0] , ni.get_prev(m_dofT[0])  + ni.get_prev(m_dofAT )*dt*(1-m_gammaf)*m_alphaf);
            }
                break;
                
            case 2:
            {
                // initial guess at start of new time step (Same Ydot)
                vec3d afp = ni.get_vec3d_prev(m_dofAW[0], m_dofAW[1], m_dofAW[2]);
                ni.set_vec3d(m_dofAW[0], m_dofAW[1], m_dofAW[2], afp);
                ni.set(m_dofAEF, ni.get_prev(m_dofAEF));
                ni.set(m_dofAT, ni.get_prev(m_dofAT));

                vec3d vfp = ni.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2]);
                ni.set_vec3d(m_dofW[0], m_dofW[1], m_dofW[2], vfp + afp*dt);
                ni.set(m_dofEF[0], ni.get_prev(m_dofEF[0]) + ni.get_prev(m_dofAEF)*dt);
                ni.set(m_dofT[0] , ni.get_prev(m_dofT[0])  + ni.get_prev(m_dofAT)*dt);
            }
                break;
                
            default:
                break;
        }
    }
    
    // apply prescribed velocities
    // we save the prescribed velocity increments in the ui vector
    vector<double>& ui = m_ui;
    zero(ui);
    int nbc = fem.BoundaryConditions();
    for (int i=0; i<nbc; ++i)
    {
        FEBoundaryCondition& bc = *fem.BoundaryCondition(i);
        if (bc.IsActive() && HasActiveDofs(bc.GetDofList())) bc.PrepStep(ui);
    }
    
    // apply prescribed DOFs for specialized surface loads
    int nsl = fem.ModelLoads();
    for (int i = 0; i < nsl; ++i)
    {
        FEModelLoad& pml = *fem.ModelLoad(i);
        if (pml.IsActive()) pml.Update();
    }

    // initialize material point data
    // NOTE: do this before the stresses are updated
    // TODO: does it matter if the stresses are updated before
    //       the material point data is initialized
    // update domain data
    for (int i=0; i<mesh.Domains(); ++i) mesh.Domain(i).PreSolveUpdate(tp);

    // update stresses
    fem.Update();
    
    for (int i = 0; i < fem.NonlinearConstraints(); ++i)
    {
        FENLConstraint* plc = fem.NonlinearConstraint(i);
        if (plc && plc->IsActive()) plc->PrepStep();
    }
    
    // see if we need to do contact augmentations
    m_baugment = false;
    for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
    {
        FEContactInterface& ci = dynamic_cast<FEContactInterface&>(*fem.SurfacePairConstraint(i));
        if (ci.IsActive() && (ci.m_laugon == 1)) m_baugment = true;
    }
    
    // see if we have to do nonlinear constraint augmentations
    if (fem.NonlinearConstraints() != 0) m_baugment = true;
}

//-----------------------------------------------------------------------------
bool FEThermoFluidSolver::Quasin()
{
    FEModel& fem = *GetFEModel();

    // convergence norms
    double    normR1;       // residual norm
    double    normE1;       // energy norm
    double    normV;        // velocity norm
    double    normv;        // velocity increment norm
    double    normRi = 0;   // initial residual norm
    double    normVi = 0;   // initial velocity norm
    double    normEi = 0;   // initial energy norm
    double    normEm = 0;   // max energy norm
    double    normDi = 0;   // initial dilatation norm
    double    normD;        // current dilatation norm
    double    normd;        // incremement dilatation norm
    double    normTi = 0;   // initial temperature norm
    double    normT;        // current temperature norm
    double    normt;        // incremement temperature norm

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
            int veq = m_neq - m_nteq;
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
        
        // solve the equations
        SolveEquations(m_ui, m_R0);

        // do the line search
        double s = DoLineSearch();

        // for sequential solve, we set one of the residual components to zero
        if (m_solve_strategy == SOLVE_SEQUENTIAL)
        {
            int veq = m_neq - m_nteq;
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
        
        // set initial convergence norms
        if (m_niter == 0)
        {
            normRi = fabs(m_R0*m_R0);
            normEi = fabs(m_ui*m_R0);
            normVi = fabs(m_vi*m_vi);
            normDi = fabs(m_di*m_di);
            normTi = fabs(m_ti*m_ti);
            normEm = normEi;
        }
        
        // calculate actual increment
        // NOTE: We don't apply the line search directly to m_ui since we need the unscaled search direction for the QN update below
        int neq = (int)m_Ui.size();
        vector<double> ui(m_ui);
        for (int i = 0; i<neq; ++i) ui[i] *= s;

        // update other increments (e.g., Lagrange multipliers)
        UpdateIncrements(m_Ui, ui, false);
        
        // calculate the norms
        normR1 = m_R1*m_R1;
        normv  = m_vi*m_vi;
        normV  = m_Vi*m_Vi;
        normd  = m_di*m_di;
        normD  = m_Di*m_Di;
        normt  = m_ti*m_ti;
        normT  = m_Ti*m_Ti;
        normE1 = fabs(m_ui*m_R1);
        
        // check for nans
        if (ISNAN(normR1)) throw NANInResidualDetected();
        
        // check residual norm
        if ((m_Rtol > 0) && (normR1 > m_Rtol*normRi)) bconv = false;
        
        // check velocity norm
        if ((m_Vtol > 0) && (normv  > (m_Vtol*m_Vtol)*normV )) bconv = false;
        
        // check dilatation norm
        if ((m_Ftol > 0) && (normd  > (m_Ftol*m_Ftol)*normD )) bconv = false;
        
        // check temperature norm
        if ((m_Ttol > 0) && (normt  > (m_Ttol*m_Ttol)*normT )) bconv = false;
        
        // check energy norm
        if ((m_Etol > 0) && (normE1 > m_Etol*normEi)) bconv = false;
        
        // check linestep size
        if ((m_lineSearch->m_LStol > 0) && (s < m_lineSearch->m_LSmin)) bconv = false;
        
        // check energy divergence
        if (normE1 > normEm) bconv = false;
        
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
        feLog("\t   temperature      %15le %15le %15le \n", normTi, normt ,(m_Ttol*m_Ttol)*normT );

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
                normTi = normt;
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
                feLog("\n*** Velocity converged. Now solving for temperature.\n");
            }
        }
        
        // check for sudden temperature change
        if (bconv && m_sudden_T_change) {
            m_sudden_T_change = false;
            throw ConcentrationChangeDetected();
        }
        else m_sudden_T_change = false;
        
        // increase iteration number
        m_niter++;
        
        // do minor iterations callbacks
        fem.DoCallback(CB_MINOR_ITERS);
    }
    while (bconv == false);
    
    // if converged we update the total velocities
    if (bconv)
    {
        UpdateIncrements(m_Ut, m_Ui, true);
        zero(m_Ui);
        zero(m_Di); zero(m_Vi); zero(m_Ti);
    }
    
    return bconv;
}

//-----------------------------------------------------------------------------
//! Calculates global stiffness matrix.

bool FEThermoFluidSolver::StiffnessMatrix(FELinearSystem& LS)
{
    FEModel& fem = *GetFEModel();

    const FETimeInfo& tp = fem.GetTime();

    // get the mesh
    FEMesh& mesh = fem.GetMesh();
    
    // calculate the stiffness matrix for each domain
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        dom.StiffnessMatrix(LS);
    }
    
    // calculate the body force stiffness matrix for each domain
    int NML = fem.ModelLoads();
    for (int j = 0; j<NML; ++j)
    {
        FEModelLoad* pml = fem.ModelLoad(j);
        FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(pml);
        FEFluidHeatSupply* phs = dynamic_cast<FEFluidHeatSupply*>(pml);
        if (pbf && pbf->IsActive())
        {
            for (int i = 0; i<pbf->Domains(); ++i)
            {
                FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(*pbf->Domain(i));
                dom.BodyForceStiffness(LS, *pbf);
            }
        }
        else if (phs && phs->IsActive())
        {
            for (int i = 0; i<phs->Domains(); ++i)
            {
                FEThermoFluidDomain3D* tdom = dynamic_cast<FEThermoFluidDomain3D*>(phs->Domain(i));
                if (tdom) tdom->HeatSupplyStiffness(LS, *phs);
            }
        }
    }
    
    // calculate contact stiffness
    ContactStiffness(LS);
    
    // calculate stiffness matrix due to model loads
    int nsl = fem.ModelLoads();
    for (int i=0; i<nsl; ++i)
    {
        FEModelLoad* pml = fem.ModelLoad(i);
        if (pml->IsActive()) pml->StiffnessMatrix(LS);
    }
    // Add mass matrix
    // loop over all domains
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        dom.MassMatrix(LS);
    }
    
    // calculate nonlinear constraint stiffness
    // note that this is the contribution of the
    // constraints enforced with augmented lagrangian
    NonLinearConstraintStiffness(LS, tp);
    
    return true;
}

//-----------------------------------------------------------------------------
//! Calculate the stiffness contribution due to nonlinear constraints
void FEThermoFluidSolver::NonLinearConstraintStiffness(FELinearSystem& LS, const FETimeInfo& tp)
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

void FEThermoFluidSolver::ContactStiffness(FELinearSystem& LS)
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
void FEThermoFluidSolver::ContactForces(FEGlobalVector& R)
{
    FEModel& fem = *GetFEModel();

    const FETimeInfo& tp = fem.GetTime();
    for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
    {
        FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));
        if (pci->IsActive()) pci->LoadVector(R, tp);
    }
}

//-----------------------------------------------------------------------------
//! calculates the residual vector
//! Note that the concentrated nodal forces are not calculated here.
//! This is because they do not depend on the geometry
//! so we only calculate them once (in Quasin) and then add them here.

bool FEThermoFluidSolver::Residual(vector<double>& R)
{
    FEModel& fem = *GetFEModel();

    // get the time information
    const FETimeInfo& tp = fem.GetTime();

    // initialize residual with concentrated nodal loads
    zero(R);
    
    // zero nodal reaction forces
    zero(m_Fr);
    
    // setup the global vector
//    FEFluidResidualVector RHS(fem, R, m_Fr);
    FEResidualVector RHS(fem, R, m_Fr);

    // get the mesh
    FEMesh& mesh = fem.GetMesh();
    
    // calculate the internal (stress) forces
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        dom.InternalForces(RHS);
    }
    

    // calculate the model loads
    for (int j = 0; j<fem.ModelLoads(); ++j)
    {
        FEModelLoad* pml = fem.ModelLoad(j);
        FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(pml);
        FEFluidHeatSupply* phs = dynamic_cast<FEFluidHeatSupply*>(pml);
        if (pbf && pbf->IsActive())
        {
            for (int i = 0; i<pbf->Domains(); ++i)
            {
                FEFluidDomain* fdom = dynamic_cast<FEFluidDomain*>(pbf->Domain(i));
                fdom->BodyForce(RHS, *pbf);
            }
        }
        else if (phs && phs->IsActive())
        {
            for (int i = 0; i<phs->Domains(); ++i)
            {
                FEThermoFluidDomain3D* tdom = dynamic_cast<FEThermoFluidDomain3D*>(phs->Domain(i));
                if (tdom) tdom->HeatSupply(RHS, *phs);
            }
        }
    }
    
    // calculate inertial forces
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& fdom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        fdom.InertialForces(RHS);
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
            mli.LoadVector(RHS);
        }
    }
    
    // set the nodal reaction forces
    // TODO: Is this a good place to do this?
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        node.set_load(m_dofW[0], 0);
        node.set_load(m_dofW[1], 0);
        node.set_load(m_dofW[2], 0);

        int n;
        if ((n = -node.m_ID[m_dofW[0]] - 2) >= 0) node.set_load(m_dofW[0], -m_Fr[n]);
        if ((n = -node.m_ID[m_dofW[1]] - 2) >= 0) node.set_load(m_dofW[1], -m_Fr[n]);
        if ((n = -node.m_ID[m_dofW[2]] - 2) >= 0) node.set_load(m_dofW[2], -m_Fr[n]);
    }
    
    // increase RHS counter
    m_nrhs++;
    
    return true;
}

//-----------------------------------------------------------------------------
//! calculate the nonlinear constraint forces
void FEThermoFluidSolver::NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp)
{
    FEModel& fem = *GetFEModel();

    int N = fem.NonlinearConstraints();
    for (int i=0; i<N; ++i)
    {
        FENLConstraint* plc = fem.NonlinearConstraint(i);
        if (plc->IsActive()) plc->LoadVector(R, tp);
    }
}

//-----------------------------------------------------------------------------
//! Serialization
void FEThermoFluidSolver::Serialize(DumpStream& ar)
{
    FENewtonSolver::Serialize(ar);

    ar & m_nrhs;
    ar & m_niter;
    ar & m_nref & m_ntotref;

    ar & m_nveq & m_ndeq & m_nteq;

    ar & m_Fr & m_Ui &m_Ut;
    ar & m_Vi & m_Di & m_Ti;
    
    if (ar.IsLoading())
    {
        m_Fr.assign(m_neq, 0);
        m_Vi.assign(m_nveq,0);
        m_Di.assign(m_ndeq,0);
        m_Ti.assign(m_nteq,0);
    }
    
    if (ar.IsShallow()) return;
    
    ar & m_rhoi & m_alphaf & m_alpham;
    ar & m_gammaf;
    ar & m_pred;

    ar & m_dofW & m_dofEF & m_dofAEF & m_dofT & m_dofAT;    
}
