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
#include "FEFluidFSISolver.h"
#include "FEBioMech/FEElasticDomain.h"
#include "FEBioMech/FEPressureLoad.h"
#include "FEBioMech/FERigidConnector.h"
#include "FEBioMech/FESlidingInterfaceBW.h"
#include "FEBioMech/FESSIShellDomain.h"
#include "FEBioMech/FEResidualVector.h"
#include "FEFluidFSIDomain.h"
#include "FEFluidDomain.h"
#include "FECore/FEModel.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include "NumCore/NumCore.h"
#include <assert.h>
#include "FECore/FEGlobalMatrix.h"
#include "FECore/sys.h"
#include "FEBioMech/FE3FieldElasticSolidDomain.h"
#include "FEBioMech/FEUncoupledMaterial.h"
#include <FEBioMech/FEBodyForce.h>
#include <FECore/FEBoundaryCondition.h>
#include <FECore/FENodalLoad.h>
#include <FECore/FESurfaceLoad.h>
#include "FEFluidResistanceBC.h"
#include "FEBackFlowStabilization.h"
#include "FEFluidNormalVelocity.h"
#include "FEFluidVelocity.h"
#include "FEFluidRotationalVelocity.h"
#include <FECore/FEModelLoad.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FELinearConstraintManager.h>
#include <FECore/DumpStream.h>
#include <FEBioMech/FESolidLinearSystem.h>
#include "FEBioFSI.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_FECORE_CLASS(FEFluidFSISolver, FENewtonSolver)
	ADD_PARAMETER(m_Dtol, "dtol"        );
	ADD_PARAMETER(m_Vtol, "vtol"        );
	ADD_PARAMETER(m_Ftol, "ftol"        );
	ADD_PARAMETER(m_rhoi, "rhoi"        );
	ADD_PARAMETER(m_pred, "predictor"   );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! FEFluidFSISolver Construction
//
FEFluidFSISolver::FEFluidFSISolver(FEModel* pfem) : FENewtonSolver(pfem), m_rigidSolver(pfem), \
m_dofU(pfem), m_dofV(pfem), m_dofSU(pfem), m_dofSV(pfem), m_dofSA(pfem),m_dofR(pfem), m_dofVF(pfem),m_dofAF(pfem),m_dofW(pfem), m_dofAW(pfem)
{
    // default values
    m_Rtol = 0.001;
    m_Etol = 0.01;
    m_Dtol = 0.001;
    m_Vtol = 0.001;
    m_Ftol = 0.001;
    m_Rmin = 1.0e-20;
	m_Rmax = 0;	// not used if zero
    
    m_ndeq = 0;
    m_nveq = 0;
    m_nfeq = 0;
    m_niter = 0;
    m_nreq = 0;

	// assume non-symmetric
    m_msymm = REAL_UNSYMMETRIC;
    
    // default Newmark parameters for rhoi = 0
    m_rhoi = 0;
    m_alpha = m_alphaf = 1;
    m_alpham = 1.5;
    m_beta = 0.5625;
    m_gamma = 1;
    m_pred = 0;
    
	// Preferred strategy is Broyden's method
	SetDefaultStrategy(QN_BROYDEN);

    // turn off checking for a zero diagonal
    CheckZeroDiagonal(false);
    
    // Allocate degrees of freedom
    DOFS& dofs = pfem->GetDOFS();
    
    int varD = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::DISPLACEMENT), VAR_VEC3);
    dofs.SetDOFName(varD, 0, "x");
    dofs.SetDOFName(varD, 1, "y");
    dofs.SetDOFName(varD, 2, "z");

    int varV = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::VELOCITY), VAR_VEC3);
    dofs.SetDOFName(varV, 0, "vx");
    dofs.SetDOFName(varV, 1, "vy");
    dofs.SetDOFName(varV, 2, "vz");
    
    int varQ = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::SHELL_ROTATION), VAR_VEC3);
    dofs.SetDOFName(varQ, 0, "u");
    dofs.SetDOFName(varQ, 1, "v");
    dofs.SetDOFName(varQ, 2, "w");

    int varSD = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::SHELL_DISPLACEMENT), VAR_VEC3);
    dofs.SetDOFName(varSD, 0, "sx");
    dofs.SetDOFName(varSD, 1, "sy");
    dofs.SetDOFName(varSD, 2, "sz");

    int varQV = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::SHELL_VELOCITY), VAR_VEC3);
    dofs.SetDOFName(varQV, 0, "svx");
    dofs.SetDOFName(varQV, 1, "svy");
    dofs.SetDOFName(varQV, 2, "svz");

    int varQA = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::SHELL_ACCELERATION), VAR_VEC3);
    dofs.SetDOFName(varQA, 0, "sax");
    dofs.SetDOFName(varQA, 1, "say");
    dofs.SetDOFName(varQA, 2, "saz");

    int varQR = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::RIGID_ROTATION), VAR_VEC3);
    dofs.SetDOFName(varQR, 0, "Ru");
    dofs.SetDOFName(varQR, 1, "Rv");
    dofs.SetDOFName(varQR, 2, "Rw");
    
    int nW = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::RELATIVE_FLUID_VELOCITY), VAR_VEC3);
    dofs.SetDOFName(nW, 0, "wx");
    dofs.SetDOFName(nW, 1, "wy");
    dofs.SetDOFName(nW, 2, "wz");
    
    int nAW = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::RELATIVE_FLUID_ACCELERATION), VAR_VEC3);
    dofs.SetDOFName(nAW, 0, "awx");
    dofs.SetDOFName(nAW, 1, "awy");
    dofs.SetDOFName(nAW, 2, "awz");
    
    int nVF = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::FLUID_VELOCITY), VAR_VEC3);
    dofs.SetDOFName(nVF, 0, "vfx");
    dofs.SetDOFName(nVF, 1, "vfy");
    dofs.SetDOFName(nVF, 2, "vfz");
    
    int nAF = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::FLUID_ACCELERATION), VAR_VEC3);
    dofs.SetDOFName(nAF, 0, "afx");
    dofs.SetDOFName(nAF, 1, "afy");
    dofs.SetDOFName(nAF, 2, "afz");
    
    int nE = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::FLUID_DILATATION), VAR_SCALAR);
    dofs.SetDOFName(nE, 0, "ef");

	int nAE = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::FLUID_DILATATION_TDERIV), VAR_SCALAR);
    dofs.SetDOFName(nAE, 0, "aef");
     
    // get the dof indices
	m_dofU.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::DISPLACEMENT));
	m_dofV.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::VELOCITY));
	m_dofSU.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::SHELL_DISPLACEMENT));
	m_dofSV.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::SHELL_VELOCITY));
	m_dofSA.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::SHELL_ACCELERATION));
	m_dofR.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::RIGID_ROTATION));
	m_dofW.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::RELATIVE_FLUID_VELOCITY));
	m_dofAW.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::RELATIVE_FLUID_ACCELERATION));
	m_dofVF.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::FLUID_VELOCITY));
	m_dofAF.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::FLUID_ACCELERATION));
    m_dofEF  = pfem->GetDOFIndex(FEBioFSI::GetVariableName(FEBioFSI::FLUID_DILATATION), 0);
    m_dofAEF = pfem->GetDOFIndex(FEBioFSI::GetVariableName(FEBioFSI::FLUID_DILATATION_TDERIV), 0);
}

//-----------------------------------------------------------------------------
FEFluidFSISolver::~FEFluidFSISolver()
{
    
}

//-----------------------------------------------------------------------------
//! Generate warnings if needed
void FEFluidFSISolver:: SolverWarnings()
{
	FEModel& fem = *GetFEModel();

    // Generate warning if rigid connectors are used with symmetric stiffness
    if (m_msymm == REAL_SYMMETRIC) {
        for (int i=0; i<fem.NonlinearConstraints(); ++i)
        {
            FENLConstraint* plc = fem.NonlinearConstraint(i);
            FERigidConnector* prc = dynamic_cast<FERigidConnector*>(plc);
            if (prc) {
                feLogWarning("Rigid connectors require non-symmetric stiffness matrix.\nSet symmetric_stiffness flag to 0 in Control section.");
                break;
            }
        }
        
        // Generate warning if sliding-elastic contact is used with symmetric stiffness
        if (fem.SurfacePairConstraints() > 0)
        {
            // loop over all contact interfaces
            for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
            {
                FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));
                FESlidingInterfaceBW* pbw = dynamic_cast<FESlidingInterfaceBW*>(pci);
                if (pbw) {
					feLogWarning("The sliding-elastic contact algorithm runs better with a non-symmetric stiffness matrix.\nYou may set symmetric_stiffness 0 to false in Control section.");
                    break;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
//! Allocates and initializes the data structures used by the FEFluidFSISolver
//
bool FEFluidFSISolver::Init()
{
    // initialize base class
    if (FENewtonSolver::Init() == false) return false;
    
    // check parameters
    if (m_Dtol <  0.0) { feLogError("dtol must be nonnegative."); return false; }
    if (m_Vtol <  0.0) { feLogError("vtol must be nonnegative."); return false; }
    if (m_Ftol <  0.0) { feLogError("ftol must be nonnegative."); return false; }
    if (m_Etol <  0.0) { feLogError("etol must be nonnegative."); return false; }
    if (m_Rtol <  0.0) { feLogError("rtol must be nonnegative."); return false; }
    
    if (m_rhoi == -1) {
        m_alphaf = m_alpham = 1.0;
    }
    else if ((m_rhoi >= 0) && (m_rhoi <= 1)) {
        m_alphaf = 1.0/(1+m_rhoi);
//        m_alpham = (3-m_rhoi)/(1+m_rhoi)/2; // 1st-order system
        m_alpham = (2-m_rhoi)/(1+m_rhoi); // 2nd-order system
    }
    else { feLogError("rhoi must be -1 or between 0 and 1.\n"); return false; }
    m_alpha = m_alphaf;
    m_beta = pow(1 + m_alpham - m_alphaf,2)/4;
    m_gamma = 0.5 + m_alpham - m_alphaf;
    
    // allocate vectors
    int neq = m_neq;
    m_Fn.assign(neq, 0);
    m_Fr.assign(neq, 0);
    m_Ui.assign(neq, 0);
    m_Ut.assign(neq, 0);
    m_di.assign(m_ndeq,0);
    m_Di.assign(m_ndeq,0);
    m_vi.assign(m_nveq,0);
    m_Vi.assign(m_nveq,0);
    m_fi.assign(m_nfeq,0);
    m_Fi.assign(m_nfeq,0);
    
    // we need to fill the total DOF vector m_Ut
    // TODO: I need to find an easier way to do this
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
    gather(m_Ut, mesh, m_dofU[0]);
    gather(m_Ut, mesh, m_dofU[1]);
    gather(m_Ut, mesh, m_dofU[2]);
    gather(m_Ut, mesh, m_dofSU[0]);
    gather(m_Ut, mesh, m_dofSU[1]);
    gather(m_Ut, mesh, m_dofSU[2]);
    gather(m_Ut, mesh, m_dofW[0]);
    gather(m_Ut, mesh, m_dofW[1]);
    gather(m_Ut, mesh, m_dofW[2]);
    gather(m_Ut, mesh, m_dofEF);
    
    SolverWarnings();
    
    return true;
}

//-----------------------------------------------------------------------------
//! Initialize equations
bool FEFluidFSISolver::InitEquations()
{
    // base class initialization
    if (FENewtonSolver::InitEquations() == false) return false;

    // store the number of equations we currently have
    m_nreq = m_neq;
    
    // Next, we assign equation numbers to the rigid body degrees of freedom
    int neq = m_rigidSolver.InitEquations(m_neq);
    if (neq == -1) return false;
    else m_neq = neq;
    
    // determine the number of velocity and dilatation equations
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
    m_ndeq = m_nveq = m_nfeq = 0;
    
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& n = mesh.Node(i);
        if (n.m_ID[m_dofU[0] ] != -1) m_ndeq++;
        if (n.m_ID[m_dofU[1] ] != -1) m_ndeq++;
        if (n.m_ID[m_dofU[2] ] != -1) m_ndeq++;
        if (n.m_ID[m_dofSU[0]] != -1) m_ndeq++;
        if (n.m_ID[m_dofSU[1]] != -1) m_ndeq++;
        if (n.m_ID[m_dofSU[2]] != -1) m_ndeq++;
        if (n.m_ID[m_dofW[0] ] != -1) m_nveq++;
        if (n.m_ID[m_dofW[1] ] != -1) m_nveq++;
        if (n.m_ID[m_dofW[2] ] != -1) m_nveq++;
        if (n.m_ID[m_dofEF   ] != -1) m_nfeq++;
    }
    
    // All initialization is done
    return true;
}

//-----------------------------------------------------------------------------
void FEFluidFSISolver::GetDisplacementData(vector<double> &xi, vector<double> &ui)
{
	FEModel& fem = *GetFEModel();

    int N = fem.GetMesh().Nodes(), nid, m = 0;
    zero(xi);
    for (int i=0; i<N; ++i)
    {
        FENode& n = fem.GetMesh().Node(i);
        nid = n.m_ID[m_dofU[0]];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            xi[m++] = ui[nid];
            assert(m <= (int) xi.size());
        }
        nid = n.m_ID[m_dofU[1]];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            xi[m++] = ui[nid];
            assert(m <= (int) xi.size());
        }
        nid = n.m_ID[m_dofU[2]];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            xi[m++] = ui[nid];
            assert(m <= (int) xi.size());
        }
        nid = n.m_ID[m_dofSU[0]];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            xi[m++] = ui[nid];
            assert(m <= (int) xi.size());
        }
        nid = n.m_ID[m_dofSU[1]];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            xi[m++] = ui[nid];
            assert(m <= (int) xi.size());
        }
        nid = n.m_ID[m_dofSU[2]];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            xi[m++] = ui[nid];
            assert(m <= (int) xi.size());
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidFSISolver::GetVelocityData(vector<double> &vi, vector<double> &ui)
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
void FEFluidFSISolver::GetDilatationData(vector<double> &ei, vector<double> &ui)
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
//! Save data to dump file

void FEFluidFSISolver::Serialize(DumpStream& ar)
{
    // Serialize parameters
    FENewtonSolver::Serialize(ar);
    
    ar & m_nreq;
    ar & m_ndeq;
    ar & m_nfeq;
    ar & m_nveq;

	ar & m_rhoi & m_alphaf & m_alpham;
	ar & m_alpha & m_beta & m_gamma;
	ar & m_pred;

	ar & m_Fn & m_Ui & m_Ut & m_Fr;
	ar & m_di & m_Di;
	ar & m_vi & m_Vi;
	ar & m_fi & m_Fi;
    
    // serialize rigid solver
    m_rigidSolver.Serialize(ar);
}

//-----------------------------------------------------------------------------
//! Update the kinematics of the model, such as nodal positions, velocities,
//! accelerations, etc.
void FEFluidFSISolver::UpdateKinematics(vector<double>& ui)
{
	FEModel& fem = *GetFEModel();

    // get the mesh
    FEMesh& mesh = fem.GetMesh();
    
    // update rigid bodies
    m_rigidSolver.UpdateRigidBodies(m_Ui, ui);
    
    // update nodes
    vector<double> U(m_Ut.size());
    for (size_t i=0; i<m_Ut.size(); ++i) U[i] = ui[i] + m_Ui[i] + m_Ut[i];
    
    scatter(U, mesh, m_dofU[0]);
    scatter(U, mesh, m_dofU[1]);
    scatter(U, mesh, m_dofU[2]);
    scatter(U, mesh, m_dofSU[0]);
    scatter(U, mesh, m_dofSU[1]);
    scatter(U, mesh, m_dofSU[2]);
    scatter(U, mesh, m_dofW[0]);
    scatter(U, mesh, m_dofW[1]);
    scatter(U, mesh, m_dofW[2]);
    scatter(U, mesh, m_dofEF);
    
    // make sure the prescribed BCs are fullfilled
    int nvel = fem.BoundaryConditions();
    for (int i=0; i<nvel; ++i)
    {
        FEBoundaryCondition& bc = *fem.BoundaryCondition(i);
        if (bc.IsActive()) bc.Update();
    }
    
    // apply prescribed DOFs for specialized surface loads
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
    
    // Update the spatial nodal positions
    // Don't update rigid nodes since they are already updated
    for (int i = 0; i<mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        if (node.m_rid == -1)
            node.m_rt = node.m_r0 + node.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2]);
    }
    
    // update time derivatives of velocity and dilatation
    // for dynamic simulations
    FEAnalysis* pstep = fem.GetCurrentStep();
    if (pstep->m_nanalysis == FE_DYNAMIC)
    {
        int N = mesh.Nodes();
		double dt = fem.GetTime().timeIncrement;
        double a = 1.0 / (m_beta*dt);
        double b = a / dt;
        double c = 1.0 - 0.5/m_beta;
        double cgi = 1 - 1.0/m_gamma;
        for (int i=0; i<N; ++i)
        {
            FENode& n = mesh.Node(i);
            
            // solid acceleration
            n.m_at = (n.m_rt - n.m_rp)*b - n.m_vp*a + n.m_ap*c;
            // solid velocity
            vec3d vt = n.m_vp + (n.m_ap*(1.0 - m_gamma) + n.m_at*m_gamma)*dt;
            n.set_vec3d(m_dofV[0], m_dofV[1], m_dofV[2], vt);
            
            // shell kinematics
            vec3d qt = n.get_vec3d(m_dofSU[0], m_dofSU[1], m_dofSU[2]);
            vec3d qp = n.get_vec3d_prev(m_dofSU[0], m_dofSU[1], m_dofSU[2]);
            vec3d vqp = n.get_vec3d_prev(m_dofSV[0], m_dofSV[1], m_dofSV[2]);
            vec3d aqp = n.get_vec3d_prev(m_dofSA[0], m_dofSA[1], m_dofSA[2]);
            vec3d aqt = (qt - qp)*b - vqp*a + aqp*c;
            vec3d vqt = vqp + (aqp*(1.0 - m_gamma) + aqt*m_gamma)*dt;
            n.set_vec3d(m_dofSA[0], m_dofSA[1], m_dofSA[2], aqt);
            n.set_vec3d(m_dofSV[0], m_dofSV[1], m_dofSV[2], vqt);

            // relative fluid velocity material time derivative (in solid frame)
            vec3d wt = n.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2]);
            vec3d wp = n.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2]);
            vec3d awt = n.get_vec3d(m_dofAW[0], m_dofAW[1], m_dofAW[2]);
            vec3d awp = n.get_vec3d_prev(m_dofAW[0], m_dofAW[1], m_dofAW[2]);
            awt = awp*cgi + (wt - wp)/(m_gamma*dt);
            n.set_vec3d(m_dofAW[0], m_dofAW[1], m_dofAW[2], awt);
            
            // fluid velocity
            vec3d vft = vt + wt;
            n.set_vec3d(m_dofVF[0], m_dofVF[1], m_dofVF[2], vft);
            // material time derivative of fluid velocity (in solid frame)
            vec3d aft = n.m_at + awt;
            n.set_vec3d(m_dofAF[0], m_dofAF[1], m_dofAF[2], aft);
            
            // dilatation time derivative
            double eft = n.get(m_dofEF);
            double efp = n.get_prev(m_dofEF);
            double aefp = n.get_prev(m_dofAEF);
            double aeft = aefp*cgi + (eft - efp)/(m_gamma*dt);
            n.set(m_dofAEF, aeft);
        }
    }
}

//-----------------------------------------------------------------------------
//! Update DOF increments
void FEFluidFSISolver::UpdateIncrements(vector<double>& Ui, vector<double>& ui, bool emap)
{
	FEModel& fem = *GetFEModel();

    // get the mesh
    FEMesh& mesh = fem.GetMesh();
    
    // update rigid bodies
    m_rigidSolver.UpdateIncrements(Ui, ui, emap);
    
    // update flexible nodes
    int n;
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        
        // displacement dofs
        // current position = initial + total at prev conv step + total increment so far + current increment
        if ((n = node.m_ID[m_dofU[0]]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofU[1]]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofU[2]]) >= 0) Ui[n] += ui[n];
        
        // rotational dofs
        if ((n = node.m_ID[m_dofSU[0]]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofSU[1]]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofSU[2]]) >= 0) Ui[n] += ui[n];

        // fluid relative velocity
        if ((n = node.m_ID[m_dofW[0]]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofW[1]]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofW[2]]) >= 0) Ui[n] += ui[n];
        
        // fluid dilatation
        if ((n = node.m_ID[m_dofEF]) >= 0) Ui[n] += ui[n];
    }
}

//-----------------------------------------------------------------------------
//! Updates the current state of the model
void FEFluidFSISolver::Update(vector<double>& ui)
{
	FEModel& fem = *GetFEModel();

    // update EAS
    UpdateEAS(ui);
    UpdateIncrementsEAS(ui, true);
    
    // update kinematics
    UpdateKinematics(ui);
    
	// update contact
	if (fem.SurfacePairConstraints() > 0) UpdateContact();

	// update constraints
	if (fem.NonlinearConstraints() > 0) UpdateConstraints();

	// update element stresses
	UpdateModel();

	// update other stuff that may depend on the deformation
	int NBL = fem.BodyLoads();
	for (int i = 0; i<NBL; ++i)
	{
		FEBodyLoad* pbl = fem.GetBodyLoad(i);
		if (pbl->IsActive()) pbl->Update();
	}
}

//-----------------------------------------------------------------------------
//! Update EAS
void FEFluidFSISolver::UpdateEAS(vector<double>& ui)
{
	FEModel& fem = *GetFEModel();

    FEMesh& mesh = fem.GetMesh();
    
    // update EAS on shell domains
    for (int i=0; i<mesh.Domains(); ++i) {
        FESSIShellDomain* sdom = dynamic_cast<FESSIShellDomain*>(&mesh.Domain(i));
        if (sdom && sdom->IsActive()) sdom->UpdateEAS(ui);
    }
}

//-----------------------------------------------------------------------------
//! Update EAS
void FEFluidFSISolver::UpdateIncrementsEAS(vector<double>& ui, const bool binc)
{
	FEModel& fem = *GetFEModel();

    FEMesh& mesh = fem.GetMesh();
    
    // update EAS on shell domains
    for (int i=0; i<mesh.Domains(); ++i) {
        FESSIShellDomain* sdom = dynamic_cast<FESSIShellDomain*>(&mesh.Domain(i));
        if (sdom && sdom->IsActive()) sdom->UpdateIncrementsEAS(ui, binc);
    }
}


//-----------------------------------------------------------------------------
//!  Updates the element stresses
void FEFluidFSISolver::UpdateModel()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	const FETimeInfo& tp = fem.GetTime();

	// update the stresses on all domains
	for (int i = 0; i<mesh.Domains(); ++i)
	{
		if (mesh.Domain(i).IsActive()) mesh.Domain(i).Update(tp);
	}
}

//-----------------------------------------------------------------------------
//! Update contact interfaces.
void FEFluidFSISolver::UpdateContact()
{
	FEModel& fem = *GetFEModel();
	// Update all contact interfaces
	const FETimeInfo& tp = fem.GetTime();
	for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));
		if (pci->IsActive()) pci->Update();
	}
}

//-----------------------------------------------------------------------------
//! Update nonlinear constraints
void FEFluidFSISolver::UpdateConstraints()
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
//!  This functions performs the Lagrange augmentations
//!  It returns true if all the augmentation have converged,
//!	otherwise it returns false
//
//! \todo There is an inherent problem with this approach. Since
//!	      Lagrangian multipliers are inherited from previous timesteps
//!       they might not be zero in case a node-surface contact breaks.
//!       The node's gap value needs to become negative to a certain value
//!       before the Lagr. multipliers dissapears.
//
bool FEFluidFSISolver::Augment()
{
	FEModel& fem = *GetFEModel();

	const FETimeInfo& tp = fem.GetTime();
    
    // Assume we will pass (can't hurt to be optimistic)
    bool bconv = true;
    
    // Do contact augmentations
    // loop over all contact interfaces
    for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
    {
        FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));
        if (pci->IsActive()) bconv = (pci->Augment(m_naug, tp) && bconv);
    }
    
    // do nonlinear constraint augmentations
    int n = fem.NonlinearConstraints();
    for (int i=0; i<n; ++i)
    {
        FENLConstraint* plc = fem.NonlinearConstraint(i);
        if (plc->IsActive()) bconv = plc->Augment(m_naug, tp) && bconv;
    }
    
    // do incompressibility multipliers for 3Field domains
    FEMesh& mesh = fem.GetMesh();
    int ND = mesh.Domains();
    for (int i=0; i<ND; ++i)
    {
        FE3FieldElasticSolidDomain* pd = dynamic_cast<FE3FieldElasticSolidDomain*>(&mesh.Domain(i));
        if (pd) bconv = (pd->Augment(m_naug) && bconv);
    }
    
    return bconv;
}

//-----------------------------------------------------------------------------
bool FEFluidFSISolver::InitStep(double time)
{
	FEModel& fem = *GetFEModel();

    // get time integration parameters
    FETimeInfo& tp = fem.GetTime();
    tp.alpha = m_alpha;
    tp.beta  = m_beta;
    tp.gamma = m_gamma;
    tp.alphaf = m_alphaf;
    tp.alpham = m_alpham;
    
    // evaluate load curve values at current (or intermediate) time
    double t = tp.currentTime;
//    double dt = tp.timeIncrement;
//    double ta = (t > 0) ? t - (1-m_alpha)*dt : m_alpha*dt;
//    return FESolver::InitStep(ta);
    return FESolver::InitStep(t);
}

//-----------------------------------------------------------------------------
//! Prepares the data for the first BFGS-iteration.
void FEFluidFSISolver::PrepStep()
{
	FEModel& fem = *GetFEModel();

	const FETimeInfo& tp = fem.GetTime();
	double dt = tp.timeIncrement;
    
    // zero total DOFs
    zero(m_Ui);
    zero(m_Vi);
    zero(m_Di);
    zero(m_Fi);
    
    // store previous mesh state
    // we need them for strain and acceleration calculations
    FEMesh& mesh = fem.GetMesh();
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& ni = mesh.Node(i);
        ni.m_rp = ni.m_rt;
        ni.m_vp = ni.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);
        ni.m_ap = ni.m_at;
		ni.UpdateValues();
        
        switch (m_pred) {
            case 0:
            {
                // initial guess at start of new time step (default)
                // solid
                ni.m_at = ni.m_ap*(1-0.5/m_beta) - ni.m_vp/(m_beta*dt);
                vec3d vs = ni.m_vp + (ni.m_at*m_gamma + ni.m_ap*(1-m_gamma))*dt;
                ni.set_vec3d(m_dofV[0], m_dofV[1], m_dofV[2], vs);
                
                // solid shell
                vec3d aqp = ni.get_vec3d_prev(m_dofSA[0], m_dofSA[1], m_dofSA[2]);
                vec3d vqp = ni.get_vec3d_prev(m_dofSV[0], m_dofSV[1], m_dofSV[2]);
                vec3d aqt = aqp*(1-0.5/m_beta) - vqp/(m_beta*dt);
                ni.set_vec3d(m_dofSA[0], m_dofSA[1], m_dofSA[2], aqt);
                vec3d vqt = vqp + (aqt*m_gamma + aqp*(1-m_gamma))*dt;
                ni.set_vec3d(m_dofSV[0], m_dofSV[1], m_dofSV[2], vqt);
                
                // fluid
                vec3d awp = ni.get_vec3d_prev(m_dofAW[0], m_dofAW[1], m_dofAW[2]);
                ni.set_vec3d(m_dofAW[0], m_dofAW[1], m_dofAW[2], awp*(m_gamma-1)/m_gamma);
                
                ni.set(m_dofAEF, ni.get_prev(m_dofAEF)*(m_gamma-1)/m_gamma);
            }
                break;
                
            case 1:
            {
                // initial guess at start of new time step (Zero Ydot)
                ni.m_at = vec3d(0,0,0);
                ni.set_vec3d(m_dofAW[0], m_dofAW[1], m_dofAW[2],vec3d(0,0,0));
                ni.set(m_dofAEF, 0);

                ni.set_vec3d(m_dofV[0], m_dofV[1], m_dofV[2], ni.m_vp + ni.m_ap*dt*(1-m_gamma)*m_alphaf);
                vec3d wp = ni.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2]);
                vec3d awp = ni.get_vec3d_prev(m_dofAW[0], m_dofAW[1], m_dofAW[2]);
                ni.set_vec3d(m_dofW[0], m_dofW[1], m_dofW[2], wp + awp*dt*(1-m_gamma)*m_alphaf);
                ni.set(m_dofEF, ni.get_prev(m_dofEF) + ni.get_prev(m_dofAEF)*dt*(1-m_gamma)*m_alphaf);
            }
                break;
                
            case 2:
            {
                // initial guess at start of new time step (Same Ydot)
                vec3d awp = ni.get_vec3d_prev(m_dofAW[0], m_dofAW[1], m_dofAW[2]);
                ni.set_vec3d(m_dofAW[0], m_dofAW[1], m_dofAW[2], awp);
                ni.set(m_dofAEF, ni.get_prev(m_dofAEF));
                
                ni.set_vec3d(m_dofV[0], m_dofV[1], m_dofV[2], ni.m_vp + ni.m_ap*dt);
                vec3d wp = ni.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2]);
                ni.set_vec3d(m_dofW[0], m_dofW[1], m_dofW[2], wp + awp*dt);
                ni.set(m_dofEF, ni.get_prev(m_dofEF) + ni.get_prev(m_dofAEF)*dt);
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
	FEResidualVector Fn(*GetFEModel(), m_Fn, dummy);
	NodalForces(Fn, tp);
	    
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
    
    // do the linear constraints
    fem.GetLinearConstraintManager().PrepStep();
    
    // initialize rigid bodies
    m_rigidSolver.PrepStep(tp, ui);
    
   // initialize material point data
    // NOTE: do this before the stresses are updated
    // TODO: does it matter if the stresses are updated before
    //       the material point data is initialized
    // update domain data
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        if (dom.IsActive()) dom.PreSolveUpdate(tp);
    }

    // update model state
	UpdateModel();
    
    // see if we need to do contact augmentations
    m_baugment = false;
    for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
    {
        FEContactInterface& ci = dynamic_cast<FEContactInterface&>(*fem.SurfacePairConstraint(i));
        if (ci.IsActive() && ci.m_blaugon) m_baugment = true;
    }
    
    // see if we need to do incompressible augmentations
    int nmat = fem.Materials();
    for (int i = 0; i<nmat; ++i)
    {
        FEUncoupledMaterial* pmi = dynamic_cast<FEUncoupledMaterial*>(fem.GetMaterial(i));
        if (pmi && pmi->m_blaugon) m_baugment = true;
    }
    
    // see if we have to do nonlinear constraint augmentations
    if (fem.NonlinearConstraints() != 0) m_baugment = true;
}

//-----------------------------------------------------------------------------
bool FEFluidFSISolver::Quasin()
{
	FEModel& fem = *GetFEModel();

    vector<double> u0(m_neq);
    vector<double> Rold(m_neq);
    
    // convergence norms
    double	normR1;		// residual norm
    double	normE1;		// energy norm
    double	normD;		// displacement norm
    double	normd;		// displacement increment norm
    double	normV;		// velocity norm
    double	normv;		// velocity increment norm
    double	normRi = 0;	// initial residual norm
    double	normDi = 0;	// initial displacement norm
    double	normVi = 0;	// initial velocity norm
    double	normEi = 0; // initial energy norm
    double	normEm = 0;	// max energy norm
    double	normFi = 0;	// initial dilatation norm
    double	normF;		// current dilatation norm
    double	normf;		// incremement dilatation norm
    
    // Get the current step
    FEAnalysis* pstep = fem.GetCurrentStep();
    
    // prepare for the first iteration
    const FETimeInfo& tp = fem.GetTime();
    PrepStep();
    
    // init QN method
	if (QNInit() == false) return false;
    
    // loop until converged or when max nr of reformations reached
	bool bconv = false;		// convergence flag
	do
    {
        feLog(" %d\n", m_niter+1);
        
        // assume we'll converge.
        bconv = true;
        
		// solve the equations (returns line search; solution stored in m_ui)
		double s = QNSolve();
            
        // extract the velocity and dilatation increments
        GetDisplacementData(m_di, m_ui);
        GetVelocityData(m_vi, m_ui);
        GetDilatationData(m_fi, m_ui);
            
        // set initial convergence norms
        if (m_niter == 0)
        {
            normRi = fabs(m_R0*m_R0);
            normEi = fabs(m_ui*m_R0);
            normDi = fabs(m_di*m_di);
            normVi = fabs(m_vi*m_vi);
            normFi = fabs(m_fi*m_fi);
            normEm = normEi;
        }
        
        // calculate norms
        // update all degrees of freedom
        for (int i=0; i<m_neq; ++i) m_Ui[i] += s*m_ui[i];
            
        // update displacements
		for (int i = 0; i<m_ndeq; ++i) m_Di[i] += s*m_di[i];
            
        // update velocities
		for (int i = 0; i<m_nveq; ++i) m_Vi[i] += s*m_vi[i];
            
        // update dilatations
		for (int i = 0; i<m_nfeq; ++i) m_Fi[i] += s*m_fi[i];
            
        // calculate the norms
        normR1 = m_R1*m_R1;
        normd  = (m_di*m_di)*(s*s);
        normD  = m_Di*m_Di;
        normv  = (m_vi*m_vi)*(s*s);
        normV  = m_Vi*m_Vi;
        normf  = (m_fi*m_fi)*(s*s);
        normF  = m_Fi*m_Fi;
        normE1 = s*fabs(m_ui*m_R1);
            
        // check for nans
        if (ISNAN(normR1)) throw NANDetected();
        
        // check residual norm
        if ((m_Rtol > 0) && (normR1 > m_Rtol*normRi)) bconv = false;
        
        // check displacement norm
        if ((m_Dtol > 0) && (normd  > (m_Dtol*m_Dtol)*normD )) bconv = false;
        
        // check velocity norm
        if ((m_Vtol > 0) && (normv  > (m_Vtol*m_Vtol)*normV )) bconv = false;
        
        // check dilatation norm
        if ((m_Ftol > 0) && (normf  > (m_Ftol*m_Ftol)*normF )) bconv = false;
        
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
        feLog("\t   displacement     %15le %15le %15le \n", normDi, normd ,(m_Dtol*m_Dtol)*normD );
        feLog("\t   velocity         %15le %15le %15le \n", normVi, normv ,(m_Vtol*m_Vtol)*normV );
        feLog("\t   dilatation       %15le %15le %15le \n", normFi, normf ,(m_Ftol*m_Ftol)*normF );
        
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
                normDi = normd;
                normVi = normv;
                normFi = normf;
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
        
        // increase iteration number
        m_niter++;
        
        // do minor iterations callbacks
        fem.DoCallback(CB_MINOR_ITERS);
    }
    while (bconv == false);
    
    // if converged we update the total velocities
    if (bconv)
    {
        UpdateIncrementsEAS(m_Ui, false);
        UpdateIncrements(m_Ut, m_Ui, true);
    }
    
    return bconv;
}

//-----------------------------------------------------------------------------
//! Calculates global stiffness matrix.

bool FEFluidFSISolver::StiffnessMatrix()
{
	FEModel& fem = *GetFEModel();

	const FETimeInfo& tp = fem.GetTime();

    // get the mesh
    FEMesh& mesh = fem.GetMesh();

	FESolidLinearSystem LS(this, &m_rigidSolver, *m_pK, m_Fd, m_ui, (m_msymm == REAL_SYMMETRIC), m_alpha, m_nreq);
    
    // calculate the stiffness matrix for each domain
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        if (dom.IsActive()) {
            FEFluidDomain* fdom = dynamic_cast<FEFluidDomain*>(&dom);
            FEFluidFSIDomain* fsidom = dynamic_cast<FEFluidFSIDomain*>(&dom);
            FEElasticDomain* edom = dynamic_cast<FEElasticDomain*>(&dom);
            if (fdom) fdom->StiffnessMatrix(LS, tp);
            else if (fsidom) fsidom->StiffnessMatrix(LS, tp);
            else if (edom) edom->StiffnessMatrix(LS);
        }
    }
    
    // calculate the body force stiffness matrix for each domain
    // but not for solid domains (since they have no mass in FSI)
	int NBL = fem.BodyLoads();
	for (int j = 0; j<NBL; ++j)
	{
		FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.GetBodyLoad(j));
		if (pbf && pbf->IsActive())
		{
			for (int i = 0; i<pbf->Domains(); ++i)
			{
				FEDomain* dom = pbf->Domain(i);
				if (dom->IsActive())
				{
					FEFluidDomain* fdom = dynamic_cast<FEFluidDomain*>(dom);
					FEFluidFSIDomain* fsidom = dynamic_cast<FEFluidFSIDomain*>(dom);
					FEElasticDomain* edom = dynamic_cast<FEElasticDomain*>(dom);
					if (fdom) fdom->BodyForceStiffness(LS, tp, *pbf);
					else if (fsidom) fsidom->BodyForceStiffness(LS, tp, *pbf);
					else if (edom)
					{
						FESolidMaterial* mat = dynamic_cast<FESolidMaterial*>(dom->GetMaterial());
						if (mat && (mat->IsRigid()==false)) edom->BodyForceStiffness(LS, *pbf);
					}
                }
            }
        }
    }
    
    // TODO: add body force stiffness for rigid bodies
    
    // Add mass matrix
    FEAnalysis* pstep = fem.GetCurrentStep();
    if (pstep->m_nanalysis == FE_DYNAMIC)
    {
        // scale factor
        double dt = tp.timeIncrement;
        double a = tp.alpham / (m_beta*dt*dt);
        // loop over all domains (except rigid)
		for (int i = 0; i<mesh.Domains(); ++i)
        {
            FEDomain& dom = mesh.Domain(i);
            if (dom.IsActive())
			{
				FEFluidDomain* fdom = dynamic_cast<FEFluidDomain*>(&dom);
				FEFluidFSIDomain* fsidom = dynamic_cast<FEFluidFSIDomain*>(&dom);
				FEElasticDomain* edom = dynamic_cast<FEElasticDomain*>(&dom);
				if (fdom) fdom->MassMatrix(LS, tp);
				else if (fsidom) fsidom->MassMatrix(LS, tp);
				else if (edom)
				{
					FESolidMaterial* mat = dynamic_cast<FESolidMaterial*>(dom.GetMaterial());
					if (mat && (mat->IsRigid() == false)) edom->MassMatrix(LS, a);
				}
            }
        }
        
        m_rigidSolver.RigidMassMatrix(LS, tp);
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
    
    // calculate nonlinear constraint stiffness
    // note that this is the contribution of the
    // constrainst enforced with augmented lagrangian
    NonLinearConstraintStiffness(LS, tp);
    
    // calculate the stiffness contributions for the rigid forces
    for (int i = 0; i<fem.ModelLoads(); ++i) fem.ModelLoad(i)->StiffnessMatrix(LS, tp);
    
    // add contributions from rigid bodies
    m_rigidSolver.StiffnessMatrix(*m_pK, tp);
    
    return true;
}

//-----------------------------------------------------------------------------
//! Calculate the stiffness contribution due to nonlinear constraints
void FEFluidFSISolver::NonLinearConstraintStiffness(FELinearSystem& LS, const FETimeInfo& tp)
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

void FEFluidFSISolver::ContactStiffness(FELinearSystem& LS)
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
void FEFluidFSISolver::ContactForces(FEGlobalVector& R)
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

bool FEFluidFSISolver::Residual(vector<double>& R)
{
	FEModel& fem = *GetFEModel();

    // get the time information
	const FETimeInfo& tp = fem.GetTime();
    
    // initialize residual with concentrated nodal loads
    R = m_Fn;
    
    // zero nodal reaction forces
    zero(m_Fr);
    
    // setup the global vector
    FEResidualVector RHS(fem, R, m_Fr);

    // zero rigid body reaction forces
    m_rigidSolver.Residual();
    
    // get the mesh
    FEMesh& mesh = fem.GetMesh();
    
    // set flag for transient or steady-state analyses
    FEAnalysis* pstep = fem.GetCurrentStep();
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        if (dom.IsActive())
		{
			FEFluidDomain* fdom = dynamic_cast<FEFluidDomain*>(&dom);
			FEFluidFSIDomain* fsidom = dynamic_cast<FEFluidFSIDomain*>(&dom);
			if (fdom) {
				if (pstep->m_nanalysis == FE_STEADY_STATE)
					fdom->SetSteadyStateAnalysis();
				else
					fdom->SetTransientAnalysis();
			}
			else if (fsidom) {
				if (pstep->m_nanalysis == FE_STEADY_STATE)
					fsidom->SetSteadyStateAnalysis();
				else
					fsidom->SetTransientAnalysis();
			}
        }
    }
    
    // calculate the internal (stress) forces
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        if (dom.IsActive())
		{
			FEFluidDomain* fdom = dynamic_cast<FEFluidDomain*>(&dom);
			FEFluidFSIDomain* fsidom = dynamic_cast<FEFluidFSIDomain*>(&dom);
			FEElasticDomain* edom = dynamic_cast<FEElasticDomain*>(&dom);
			if (fdom) fdom->InternalForces(RHS, tp);
			else if (fsidom) fsidom->InternalForces(RHS, tp);
			else if (edom)
			{
				FESolidMaterial* mat = dynamic_cast<FESolidMaterial*>(dom.GetMaterial());
				if (mat && (mat->IsRigid() == false)) edom->InternalForces(RHS);
			}
        }
    }
    
    // calculate the body forces
	for (int j = 0; j<fem.BodyLoads(); ++j)
	{
		FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.GetBodyLoad(j));
		if (pbf && pbf->IsActive())
		{
			for (int i = 0; i<pbf->Domains(); ++i)
			{
				FEDomain* dom = pbf->Domain(i);
				if (dom->IsActive())
				{
					FEFluidDomain* fdom = dynamic_cast<FEFluidDomain*>(dom);
					FEFluidFSIDomain* fsidom = dynamic_cast<FEFluidFSIDomain*>(dom);
					FEElasticDomain* edom = dynamic_cast<FEElasticDomain*>(dom);
					if (fdom) fdom->BodyForce(RHS, tp, *pbf);
					else if (fsidom) fsidom->BodyForce(RHS, tp, *pbf);
					else if (edom)
					{
						FESolidMaterial* mat = dynamic_cast<FESolidMaterial*>(dom->GetMaterial());
						if (mat && (mat->IsRigid()==false)) edom->BodyForce(RHS, *pbf);
					}
				}
			}
		}
    }
    
    // calculate body forces for rigid bodies
    for (int j=0; j<fem.BodyLoads(); ++j)
    {
        FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.GetBodyLoad(j));
		if (pbf && pbf->IsActive())
			m_rigidSolver.BodyForces(RHS, tp, *pbf);
    }

    // allocate F
    vector<double> F;
    
    // calculate inertial forces
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        if (dom.IsActive())
		{
			FEFluidDomain* fdom = dynamic_cast<FEFluidDomain*>(&dom);
			FEFluidFSIDomain* fsidom = dynamic_cast<FEFluidFSIDomain*>(&dom);
			FEElasticDomain* edom = dynamic_cast<FEElasticDomain*>(&dom);
			if (fdom) fdom->InertialForces(RHS, tp);
			else if (fsidom) fsidom->InertialForces(RHS, tp);
			else if (edom && (pstep->m_nanalysis == FE_DYNAMIC))
			{
				FESolidMaterial* mat = dynamic_cast<FESolidMaterial*>(dom.GetMaterial());
				if (mat && (mat->IsRigid()==false)) edom->InertialForces(RHS, F);
			}
        }
    }
    
    // update rigid bodies
    if (pstep->m_nanalysis == FE_DYNAMIC) m_rigidSolver.InertialForces(RHS, tp);

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
        if ((n = -node.m_ID[m_dofU[0]]-2) >= 0) node.m_Fr.x = -m_Fr[n];
        if ((n = -node.m_ID[m_dofU[1]]-2) >= 0) node.m_Fr.y = -m_Fr[n];
        if ((n = -node.m_ID[m_dofU[2]]-2) >= 0) node.m_Fr.z = -m_Fr[n];
    }
    
    // increase RHS counter
    m_nrhs++;
    
    return true;
}

//-----------------------------------------------------------------------------
//! calculate the nonlinear constraint forces
void FEFluidFSISolver::NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp)
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
//! calculates the concentrated nodal forces

void FEFluidFSISolver::NodalForces(FEGlobalVector& R, const FETimeInfo& tp)
{
    // loop over nodal loads
	FEModel& fem = *GetFEModel();
	int NNL = fem.NodalLoads();
    for (int i=0; i<NNL; ++i)
    {
        FENodalLoad& fc = *fem.NodalLoad(i);
		if (fc.IsActive()) fc.Residual(R, tp);
    }
}
