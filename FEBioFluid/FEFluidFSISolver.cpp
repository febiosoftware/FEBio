//
//  FEFluidFSISolver.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 8/13/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

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
#include "FECore/BFGSSolver2.h"
#include "FECore/sys.h"
#include "FEBioMech/FE3FieldElasticSolidDomain.h"
#include "FEBioMech/FEUncoupledMaterial.h"
#include <FEBioMech/FEBodyForce.h>
#include <FECore/BC.h>
#include <FECore/FESurfaceLoad.h>
#include "FEFluidResistanceBC.h"
#include "FEBackFlowStabilization.h"
#include "FEFluidNormalVelocity.h"
#include "FEFluidVelocity.h"
#include "FEFluidRotationalVelocity.h"
#include <FECore/FEModelLoad.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FELinearConstraintManager.h>

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_PARAMETER_LIST(FEFluidFSISolver, FENewtonSolver)
ADD_PARAMETER(m_Dtol         , FE_PARAM_DOUBLE, "dtol"        );
ADD_PARAMETER(m_Vtol         , FE_PARAM_DOUBLE, "vtol"        );
ADD_PARAMETER(m_Ftol         , FE_PARAM_DOUBLE, "ftol"        );
ADD_PARAMETER(m_Etol         , FE_PARAM_DOUBLE, "etol"        );
ADD_PARAMETER(m_Rtol         , FE_PARAM_DOUBLE, "rtol"        );
ADD_PARAMETER(m_Rmin         , FE_PARAM_DOUBLE, "min_residual");
ADD_PARAMETER(m_bdivreform   , FE_PARAM_BOOL  , "diverge_reform");
ADD_PARAMETER(m_bdoreforms   , FE_PARAM_BOOL  , "do_reforms"  );
ADD_PARAMETER(m_bsymm        , FE_PARAM_BOOL  , "symmetric_stiffness");
ADD_PARAMETER(m_breformtimestep, FE_PARAM_BOOL, "reform_each_time_step");
ADD_PARAMETER(m_rhoi         , FE_PARAM_DOUBLE, "rhoi"        );
ADD_PARAMETER(m_pred         , FE_PARAM_INT   , "predictor"   );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! FEFluidFSISolver Construction
//
FEFluidFSISolver::FEFluidFSISolver(FEModel* pfem) : FENewtonSolver(pfem), m_rigidSolver(pfem)
{
    // default values
    m_Rtol = 0.001;
    m_Etol = 0.01;
    m_Dtol = 0.001;
    m_Vtol = 0.001;
    m_Ftol = 0.001;
    m_Rmin = 1.0e-20;
    
    m_ndeq = 0;
    m_nveq = 0;
    m_nfeq = 0;
    m_niter = 0;
    m_nreq = 0;

    m_bsymm = false;
    m_bdivreform = true;
    m_bdoreforms = true;
    m_breformtimestep = true;
    
    // default Newmark parameters for rhoi = 0
    m_rhoi = 0;
    m_alpha = m_alphaf = 1;
    m_alpham = 1.5;
    m_beta = 0.5625;
    m_gamma = 1;
    m_pred = 0;
    
    m_baugment = false;
    
    // a different solution strategy is used here
    m_nqnsolver = QN_BROYDEN;
    
    // turn off checking for a zero diagonal
    CheckZeroDiagonal(false);
    
    // Allocate degrees of freedom
    DOFS& dofs = pfem->GetDOFS();
    
    int varD = dofs.AddVariable("displacement", VAR_VEC3);
    dofs.SetDOFName(varD, 0, "x");
    dofs.SetDOFName(varD, 1, "y");
    dofs.SetDOFName(varD, 2, "z");
    int varV = dofs.AddVariable("velocity", VAR_VEC3);
    dofs.SetDOFName(varV, 0, "vx");
    dofs.SetDOFName(varV, 1, "vy");
    dofs.SetDOFName(varV, 2, "vz");
    
    int varQ = dofs.AddVariable("rotation", VAR_VEC3);
    dofs.SetDOFName(varQ, 0, "u");
    dofs.SetDOFName(varQ, 1, "v");
    dofs.SetDOFName(varQ, 2, "w");
    int varQP = dofs.AddVariable("previous rotation", VAR_VEC3);
    dofs.SetDOFName(varQP, 0, "up");
    dofs.SetDOFName(varQP, 1, "vp");
    dofs.SetDOFName(varQP, 2, "wp");
    
    int varSD = dofs.AddVariable("shell displacement", VAR_VEC3);
    dofs.SetDOFName(varSD, 0, "sx");
    dofs.SetDOFName(varSD, 1, "sy");
    dofs.SetDOFName(varSD, 2, "sz");
    int varSDP = dofs.AddVariable("previous shell displacement", VAR_VEC3);
    dofs.SetDOFName(varSDP, 0, "sxp");
    dofs.SetDOFName(varSDP, 1, "syp");
    dofs.SetDOFName(varSDP, 2, "szp");
    int varQV = dofs.AddVariable("shell velocity", VAR_VEC3);
    dofs.SetDOFName(varQV, 0, "svx");
    dofs.SetDOFName(varQV, 1, "svy");
    dofs.SetDOFName(varQV, 2, "svz");
    int varQVP = dofs.AddVariable("previous shell velocity", VAR_VEC3);
    dofs.SetDOFName(varQVP, 0, "svxp");
    dofs.SetDOFName(varQVP, 1, "svyp");
    dofs.SetDOFName(varQVP, 2, "svzp");
    
    int varQA = dofs.AddVariable("shell acceleration", VAR_VEC3);
    dofs.SetDOFName(varQA, 0, "sax");
    dofs.SetDOFName(varQA, 1, "say");
    dofs.SetDOFName(varQA, 2, "saz");
    int varQAP = dofs.AddVariable("previous shell acceleration", VAR_VEC3);
    dofs.SetDOFName(varQAP, 0, "saxp");
    dofs.SetDOFName(varQAP, 1, "sayp");
    dofs.SetDOFName(varQAP, 2, "sazp");
    
    int varQR = dofs.AddVariable("rigid rotation", VAR_VEC3);
    dofs.SetDOFName(varQR, 0, "Ru");
    dofs.SetDOFName(varQR, 1, "Rv");
    dofs.SetDOFName(varQR, 2, "Rw");
    
    int nW = dofs.AddVariable("relative fluid velocity", VAR_VEC3);
    dofs.SetDOFName(nW, 0, "wx");
    dofs.SetDOFName(nW, 1, "wy");
    dofs.SetDOFName(nW, 2, "wz");
    
    int nWP = dofs.AddVariable("previous relative fluid velocity", VAR_VEC3);
    dofs.SetDOFName(nWP, 0, "wxp");
    dofs.SetDOFName(nWP, 1, "wyp");
    dofs.SetDOFName(nWP, 2, "wzp");
    
    int nAW = dofs.AddVariable("relative fluid acceleration", VAR_VEC3);
    dofs.SetDOFName(nAW, 0, "awx");
    dofs.SetDOFName(nAW, 1, "awy");
    dofs.SetDOFName(nAW, 2, "awz");
    
    int nAWP = dofs.AddVariable("previous relative fluid acceleration", VAR_VEC3);
    dofs.SetDOFName(nAWP, 0, "awxp");
    dofs.SetDOFName(nAWP, 1, "awyp");
    dofs.SetDOFName(nAWP, 2, "awzp");

    int nVF = dofs.AddVariable("fluid velocity", VAR_VEC3);
    dofs.SetDOFName(nVF, 0, "vfx");
    dofs.SetDOFName(nVF, 1, "vfy");
    dofs.SetDOFName(nVF, 2, "vfz");
    
    int nAF = dofs.AddVariable("fluid acceleration", VAR_VEC3);
    dofs.SetDOFName(nAF, 0, "afx");
    dofs.SetDOFName(nAF, 1, "afy");
    dofs.SetDOFName(nAF, 2, "afz");
    
    int nE = dofs.AddVariable("fluid dilation", VAR_SCALAR);
    dofs.SetDOFName(nE, 0, "ef");
    int nEP = dofs.AddVariable("previous fluid dilation", VAR_SCALAR);
    dofs.SetDOFName(nEP, 0, "efp");
    int nAE = dofs.AddVariable("fluid dilation tderiv", VAR_SCALAR);
    dofs.SetDOFName(nAE, 0, "aef");
    int nAEP = dofs.AddVariable("previous fluid dilation tderiv", VAR_SCALAR);
    dofs.SetDOFName(nAEP, 0, "aefp");
    
    // get the dof indices
    m_dofX  = pfem->GetDOFIndex("x");
    m_dofY  = pfem->GetDOFIndex("y");
    m_dofZ  = pfem->GetDOFIndex("z");
    m_dofVX = pfem->GetDOFIndex("vx");
    m_dofVY = pfem->GetDOFIndex("vy");
    m_dofVZ = pfem->GetDOFIndex("vz");
    
    m_dofSX  = pfem->GetDOFIndex("sx");
    m_dofSY  = pfem->GetDOFIndex("sy");
    m_dofSZ  = pfem->GetDOFIndex("sz");
    m_dofSXP  = pfem->GetDOFIndex("sxp");
    m_dofSYP  = pfem->GetDOFIndex("syp");
    m_dofSZP  = pfem->GetDOFIndex("szp");
    m_dofSVX  = pfem->GetDOFIndex("svx");
    m_dofSVY  = pfem->GetDOFIndex("svy");
    m_dofSVZ  = pfem->GetDOFIndex("svz");
    m_dofSVXP  = pfem->GetDOFIndex("svxp");
    m_dofSVYP  = pfem->GetDOFIndex("svyp");
    m_dofSVZP  = pfem->GetDOFIndex("svzp");
    m_dofSAX  = pfem->GetDOFIndex("sax");
    m_dofSAY  = pfem->GetDOFIndex("say");
    m_dofSAZ  = pfem->GetDOFIndex("saz");
    m_dofSAXP  = pfem->GetDOFIndex("saxp");
    m_dofSAYP  = pfem->GetDOFIndex("sayp");
    m_dofSAZP  = pfem->GetDOFIndex("sazp");

    m_dofRU = pfem->GetDOFIndex("Ru");
    m_dofRV = pfem->GetDOFIndex("Rv");
    m_dofRW = pfem->GetDOFIndex("Rw");
    
    m_dofWX = pfem->GetDOFIndex("wx");
    m_dofWY = pfem->GetDOFIndex("wy");
    m_dofWZ = pfem->GetDOFIndex("wz");
    
    m_dofWXP = pfem->GetDOFIndex("wxp");
    m_dofWYP = pfem->GetDOFIndex("wyp");
    m_dofWZP = pfem->GetDOFIndex("wzp");
    
    m_dofAWX = pfem->GetDOFIndex("awx");
    m_dofAWY = pfem->GetDOFIndex("awy");
    m_dofAWZ = pfem->GetDOFIndex("awz");
    
    m_dofAWXP = pfem->GetDOFIndex("awxp");
    m_dofAWYP = pfem->GetDOFIndex("awyp");
    m_dofAWZP = pfem->GetDOFIndex("awzp");

    m_dofVFX = pfem->GetDOFIndex("vfx");
    m_dofVFY = pfem->GetDOFIndex("vfy");
    m_dofVFZ = pfem->GetDOFIndex("vfz");
    
    m_dofAFX = pfem->GetDOFIndex("afx");
    m_dofAFY = pfem->GetDOFIndex("afy");
    m_dofAFZ = pfem->GetDOFIndex("afz");
    
    m_dofEF  = pfem->GetDOFIndex("ef");
    m_dofEFP  = pfem->GetDOFIndex("efp");
    m_dofAEF = pfem->GetDOFIndex("aef");
    m_dofAEFP = pfem->GetDOFIndex("aefp");
}

//-----------------------------------------------------------------------------
FEFluidFSISolver::~FEFluidFSISolver()
{
    
}

//-----------------------------------------------------------------------------
//! Generate warnings if needed
void FEFluidFSISolver:: SolverWarnings()
{
    // Generate warning if rigid connectors are used with symmetric stiffness
    if (m_bsymm) {
        for (int i=0; i<m_fem.NonlinearConstraints(); ++i)
        {
            FENLConstraint* plc = m_fem.NonlinearConstraint(i);
            FERigidConnector* prc = dynamic_cast<FERigidConnector*>(plc);
            if (prc) {
                felog.printbox("WARNING", "Rigid connectors require non-symmetric stiffness matrix.\nSet symmetric_stiffness flag to false in Control section.");
                break;
            }
        }
        
        // Generate warning if sliding-elastic contact is used with symmetric stiffness
        if (m_fem.SurfacePairConstraints() > 0)
        {
            // loop over all contact interfaces
            for (int i = 0; i<m_fem.SurfacePairConstraints(); ++i)
            {
                FEContactInterface* pci = dynamic_cast<FEContactInterface*>(m_fem.SurfacePairConstraint(i));
                FESlidingInterfaceBW* pbw = dynamic_cast<FESlidingInterfaceBW*>(pci);
                if (pbw) {
                    felog.printbox("WARNING", "The sliding-elastic contact algorithm \nruns better with a non-symmetric stiffness matrix.\nYou may set symmetric_stiffness flag to false in Control section.");
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
    if (m_Dtol <  0.0) { felog.printf("Error: dtol must be nonnegative.\n"); return false; }
    if (m_Vtol <  0.0) { felog.printf("Error: vtol must be nonnegative.\n"); return false; }
    if (m_Ftol <  0.0) { felog.printf("Error: ftol must be nonnegative.\n"); return false; }
    if (m_Etol <  0.0) { felog.printf("Error: etol must be nonnegative.\n"); return false; }
    if (m_Rtol <  0.0) { felog.printf("Error: rtol must be nonnegative.\n"); return false; }
    if (m_Rmin <  0.0) { felog.printf("Error: min_residual must be nonnegative.\n"  ); return false; }
    
    if (m_rhoi == -1) {
        m_alphaf = m_alpham = 1.0;
    }
    else if ((m_rhoi >= 0) && (m_rhoi <= 1)) {
        m_alphaf = 1.0/(1+m_rhoi);
//        m_alpham = (3-m_rhoi)/(1+m_rhoi)/2; // 1st-order system
        m_alpham = (2-m_rhoi)/(1+m_rhoi); // 2nd-order system
    }
    else { felog.printf("Error: rhoi must be -1 or between 0 and 1.\n"); return false; }
    m_alpha = m_alphaf;
    m_beta = pow(1 + m_alpham - m_alphaf,2)/4;
    m_gamma = 0.5 + m_alpham - m_alphaf;
    
    // allocate vectors
    int neq = m_neq;
    m_Fn.assign(neq, 0);
    m_Fd.assign(neq, 0);
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
    FEMesh& mesh = m_fem.GetMesh();
    gather(m_Ut, mesh, m_dofX);
    gather(m_Ut, mesh, m_dofY);
    gather(m_Ut, mesh, m_dofZ);
    gather(m_Ut, mesh, m_dofSX);
    gather(m_Ut, mesh, m_dofSY);
    gather(m_Ut, mesh, m_dofSZ);
    gather(m_Ut, mesh, m_dofWX);
    gather(m_Ut, mesh, m_dofWY);
    gather(m_Ut, mesh, m_dofWZ);
    gather(m_Ut, mesh, m_dofEF);
    
    // TODO: move this somewhere else
    int nsl = m_fem.SurfaceLoads();
    for (int i=0; i<nsl; ++i)
    {
        FESurfaceLoad* psl = m_fem.SurfaceLoad(i);
        if (psl->IsActive()) {
            FEFluidResistanceBC* pfr = dynamic_cast<FEFluidResistanceBC*>(psl);
            FEFluidNormalVelocity* pnv = dynamic_cast<FEFluidNormalVelocity*>(psl);
            FEFluidRotationalVelocity* prv = dynamic_cast<FEFluidRotationalVelocity*>(psl);
            FEFluidVelocity* pv = dynamic_cast<FEFluidVelocity*>(psl);
            if (pfr) pfr->MarkDilatation();
            else if (pnv && pnv->m_bpv) pnv->MarkVelocity();
            else if (prv) prv->MarkVelocity();
            else if (pv) pv->MarkVelocity();
        }
    }
    
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
    FEMesh& mesh = m_fem.GetMesh();
    m_ndeq = m_nveq = m_nfeq = 0;
    
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& n = mesh.Node(i);
        if (n.m_ID[m_dofX ] != -1) m_ndeq++;
        if (n.m_ID[m_dofY ] != -1) m_ndeq++;
        if (n.m_ID[m_dofZ ] != -1) m_ndeq++;
        if (n.m_ID[m_dofSX] != -1) m_ndeq++;
        if (n.m_ID[m_dofSY] != -1) m_ndeq++;
        if (n.m_ID[m_dofSZ] != -1) m_ndeq++;
        if (n.m_ID[m_dofWX] != -1) m_nveq++;
        if (n.m_ID[m_dofWY] != -1) m_nveq++;
        if (n.m_ID[m_dofWZ] != -1) m_nveq++;
        if (n.m_ID[m_dofEF] != -1) m_nfeq++;
    }
    
    // All initialization is done
    return true;
}

//-----------------------------------------------------------------------------
void FEFluidFSISolver::GetDisplacementData(vector<double> &xi, vector<double> &ui)
{
    int N = m_fem.GetMesh().Nodes(), nid, m = 0;
    zero(xi);
    for (int i=0; i<N; ++i)
    {
        FENode& n = m_fem.GetMesh().Node(i);
        nid = n.m_ID[m_dofX];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            xi[m++] = ui[nid];
            assert(m <= (int) xi.size());
        }
        nid = n.m_ID[m_dofY];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            xi[m++] = ui[nid];
            assert(m <= (int) xi.size());
        }
        nid = n.m_ID[m_dofZ];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            xi[m++] = ui[nid];
            assert(m <= (int) xi.size());
        }
        nid = n.m_ID[m_dofSX];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            xi[m++] = ui[nid];
            assert(m <= (int) xi.size());
        }
        nid = n.m_ID[m_dofSY];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            xi[m++] = ui[nid];
            assert(m <= (int) xi.size());
        }
        nid = n.m_ID[m_dofSZ];
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
    int N = m_fem.GetMesh().Nodes(), nid, m = 0;
    zero(vi);
    for (int i=0; i<N; ++i)
    {
        FENode& n = m_fem.GetMesh().Node(i);
        nid = n.m_ID[m_dofWX];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            vi[m++] = ui[nid];
            assert(m <= (int) vi.size());
        }
        nid = n.m_ID[m_dofWY];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            vi[m++] = ui[nid];
            assert(m <= (int) vi.size());
        }
        nid = n.m_ID[m_dofWZ];
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
    int N = m_fem.GetMesh().Nodes(), nid, m = 0;
    zero(ei);
    for (int i=0; i<N; ++i)
    {
        FENode& n = m_fem.GetMesh().Node(i);
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
    
    if (ar.IsSaving())
    {
        ar << m_nrhs;
        ar << m_niter;
        ar << m_nref << m_ntotref;
        ar << m_naug;
        ar << m_nreq;
        ar << m_ndeq;
        ar << m_nfeq;
        ar << m_nveq;
    }
    else
    {
        ar >> m_nrhs;
        ar >> m_niter;
        ar >> m_nref >> m_ntotref;
        ar >> m_naug;
        ar >> m_nreq;
        ar >> m_ndeq;
        ar >> m_nfeq;
        ar >> m_nveq;
        
        // initialize data structures
        // (only when number of equations is non-zero.
        // This can be zero in a multi-step analysis for steps that have not yet been initialized.)
        if (m_neq > 0) Init();
    }
    
    // serialize rigid solver
    m_rigidSolver.Serialize(ar);
}

//-----------------------------------------------------------------------------
//! Update the kinematics of the model, such as nodal positions, velocities,
//! accelerations, etc.
void FEFluidFSISolver::UpdateKinematics(vector<double>& ui)
{
    // get the mesh
    FEMesh& mesh = m_fem.GetMesh();
    
    // update rigid bodies
    m_rigidSolver.UpdateRigidBodies(m_Ui, ui);
    
    // update nodes
    vector<double> U(m_Ut.size());
    for (size_t i=0; i<m_Ut.size(); ++i) U[i] = ui[i] + m_Ui[i] + m_Ut[i];
    
    scatter(U, mesh, m_dofX);
    scatter(U, mesh, m_dofY);
    scatter(U, mesh, m_dofZ);
    scatter(U, mesh, m_dofSX);
    scatter(U, mesh, m_dofSY);
    scatter(U, mesh, m_dofSZ);
    scatter(U, mesh, m_dofWX);
    scatter(U, mesh, m_dofWY);
    scatter(U, mesh, m_dofWZ);
    scatter(U, mesh, m_dofEF);
    
    // make sure the prescribed BCs are fullfilled
    int nvel = m_fem.PrescribedBCs();
    for (int i=0; i<nvel; ++i)
    {
        FEPrescribedBC& dc = *m_fem.PrescribedBC(i);
        if (dc.IsActive()) dc.Update();
    }
    
    // enforce the linear constraints
    // TODO: do we really have to do this? Shouldn't the algorithm
    // already guarantee that the linear constraints are satisfied?
    FELinearConstraintManager& LCM = m_fem.GetLinearConstraintManager();
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
            node.m_rt = node.m_r0 + node.get_vec3d(m_dofX, m_dofY, m_dofZ);
    }
    
    // update time derivatives of velocity and dilatation
    // for dynamic simulations
    FEAnalysis* pstep = m_fem.GetCurrentStep();
    if (pstep->m_nanalysis == FE_DYNAMIC)
    {
        int N = mesh.Nodes();
        double dt = pstep->m_dt;
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
            n.set_vec3d(m_dofVX, m_dofVY, m_dofVZ, vt);
            
            // shell kinematics
            vec3d qt = n.get_vec3d(m_dofSX, m_dofSY, m_dofSZ);
            vec3d qp = n.get_vec3d(m_dofSXP, m_dofSYP, m_dofSZP);
            vec3d vqp = n.get_vec3d(m_dofSVXP, m_dofSVYP, m_dofSVZP);
            vec3d aqp = n.get_vec3d(m_dofSAXP, m_dofSAYP, m_dofSAZP);
            vec3d aqt = (qt - qp)*b - vqp*a + aqp*c;
            vec3d vqt = vqp + (aqp*(1.0 - m_gamma) + aqt*m_gamma)*dt;
            n.set_vec3d(m_dofSAX, m_dofSAY, m_dofSAZ, aqt);
            n.set_vec3d(m_dofSVX, m_dofSVY, m_dofSVZ, vqt);

            // relative fluid velocity material time derivative (in solid frame)
            vec3d wt = n.get_vec3d(m_dofWX, m_dofWY, m_dofWZ);
            vec3d wp = n.get_vec3d(m_dofWXP, m_dofWYP, m_dofWZP);
            vec3d awt = n.get_vec3d(m_dofAWX, m_dofAWY, m_dofAWZ);
            vec3d awp = n.get_vec3d(m_dofAWXP, m_dofAWYP, m_dofAWZP);
            awt = awp*cgi + (wt - wp)/(m_gamma*dt);
            n.set_vec3d(m_dofAWX, m_dofAWY, m_dofAWZ, awt);
            
            // fluid velocity
            vec3d vft = vt + wt;
            n.set_vec3d(m_dofVFX, m_dofVFY, m_dofVFZ, vft);
            // material time derivative of fluid velocity (in solid frame)
            vec3d aft = n.m_at + awt;
            n.set_vec3d(m_dofAFX, m_dofAFY, m_dofAFZ, aft);
            
            // dilatation time derivative
            double eft = n.get(m_dofEF);
            double efp = n.get(m_dofEFP);
            double aefp = n.get(m_dofAEFP);
            double aeft = aefp*cgi + (eft - efp)/(m_gamma*dt);
            n.set(m_dofAEF, aeft);
        }
    }
}

//-----------------------------------------------------------------------------
//! Update DOF increments
void FEFluidFSISolver::UpdateIncrements(vector<double>& Ui, vector<double>& ui, bool emap)
{
    // get the mesh
    FEMesh& mesh = m_fem.GetMesh();
    
    // update rigid bodies
    m_rigidSolver.UpdateIncrements(Ui, ui, emap);
    
    // update flexible nodes
    int n;
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        
        // displacement dofs
        // current position = initial + total at prev conv step + total increment so far + current increment
        if ((n = node.m_ID[m_dofX]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofY]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofZ]) >= 0) Ui[n] += ui[n];
        
        // rotational dofs
        if ((n = node.m_ID[m_dofSX]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofSY]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofSZ]) >= 0) Ui[n] += ui[n];

        // fluid relative velocity
        if ((n = node.m_ID[m_dofWX]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofWY]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofWZ]) >= 0) Ui[n] += ui[n];
        
        // fluid dilatation
        if ((n = node.m_ID[m_dofEF]) >= 0) Ui[n] += ui[n];
    }
}

//-----------------------------------------------------------------------------
//! Updates the current state of the model
void FEFluidFSISolver::Update(vector<double>& ui)
{
    TimerTracker t(m_UpdateTime);
    
    // update EAS
    UpdateEAS(ui);
    UpdateIncrementsEAS(ui, true);
    
    // update kinematics
    UpdateKinematics(ui);
    
    // TODO: move this somewhere else
    int nsl = m_fem.SurfaceLoads();
    for (int i=0; i<nsl; ++i)
    {
        FESurfaceLoad* psl = m_fem.SurfaceLoad(i);
        if (psl->IsActive()) {
            FEFluidResistanceBC* pfr = dynamic_cast<FEFluidResistanceBC*>(psl);
            FEFluidNormalVelocity* pnv = dynamic_cast<FEFluidNormalVelocity*>(psl);
            FEFluidRotationalVelocity* prv = dynamic_cast<FEFluidRotationalVelocity*>(psl);
            FEFluidVelocity* pv = dynamic_cast<FEFluidVelocity*>(psl);
            if (pfr) pfr->SetDilatation();
            else if (pnv && pnv->m_bpv) pnv->SetVelocity();
            else if (prv) prv->SetVelocity();
            else if (pv) pv->SetVelocity();
        }
    }
    
    // update contact
    if (m_fem.SurfacePairConstraints() > 0) UpdateContact();
    
    // update constraints
    if (m_fem.NonlinearConstraints() > 0) UpdateConstraints();
    
    // update element stresses
    UpdateStresses();
    
    // update other stuff that may depend on the deformation
    int NBL = m_fem.BodyLoads();
    for (int i=0; i<NBL; ++i)
    {
        FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(i));
        if (pbf) pbf->Update();
    }
}

//-----------------------------------------------------------------------------
//! Update EAS
void FEFluidFSISolver::UpdateEAS(vector<double>& ui)
{
    FEMesh& mesh = m_fem.GetMesh();
    
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
    FEMesh& mesh = m_fem.GetMesh();
    
    // update EAS on shell domains
    for (int i=0; i<mesh.Domains(); ++i) {
        FESSIShellDomain* sdom = dynamic_cast<FESSIShellDomain*>(&mesh.Domain(i));
        if (sdom && sdom->IsActive()) sdom->UpdateIncrementsEAS(ui, binc);
    }
}

//-----------------------------------------------------------------------------
//!  Updates the element stresses
void FEFluidFSISolver::UpdateStresses()
{
    FEMesh& mesh = m_fem.GetMesh();
    FETimeInfo tp = m_fem.GetTime();
    tp.alpha = m_alpha;
    tp.beta  = m_beta;
    tp.gamma = m_gamma;
    tp.alphaf = m_alphaf;
    tp.alpham = m_alpham;
    
    // update the stresses on all domains
    for (int i=0; i<mesh.Domains(); ++i)
    {
        if (mesh.Domain(i).IsActive()) mesh.Domain(i).Update(tp);
    }
}

//-----------------------------------------------------------------------------
//! Update contact interfaces.
void FEFluidFSISolver::UpdateContact()
{
    // Update all contact interfaces
    FETimeInfo tp = GetFEModel().GetTime();
    tp.alpha = m_alpha;
    tp.beta  = m_beta;
    tp.gamma = m_gamma;
    tp.alphaf = m_alphaf;
    tp.alpham = m_alpham;
    for (int i = 0; i<m_fem.SurfacePairConstraints(); ++i)
    {
        FEContactInterface* pci = dynamic_cast<FEContactInterface*>(m_fem.SurfacePairConstraint(i));
        if (pci->IsActive()) pci->Update(m_niter, tp);
    }
}

//-----------------------------------------------------------------------------
//! Update nonlinear constraints
void FEFluidFSISolver::UpdateConstraints()
{
    FETimeInfo tp = m_fem.GetTime();
    tp.alpha = m_alpha;
    tp.beta  = m_beta;
    tp.gamma = m_gamma;
    tp.alphaf = m_alphaf;
    tp.alpham = m_alpham;
    tp.currentIteration = m_niter;
    
    // Update all nonlinear constraints
    for (int i=0; i<m_fem.NonlinearConstraints(); ++i)
    {
        FENLConstraint* pci = m_fem.NonlinearConstraint(i);
        if (pci->IsActive()) pci->Update(m_niter, tp);
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
    FETimeInfo tp = m_fem.GetTime();
    tp.alpha = m_alpha;
    tp.beta  = m_beta;
    tp.gamma = m_gamma;
    tp.alphaf = m_alphaf;
    tp.alpham = m_alpham;
    
    // Assume we will pass (can't hurt to be optimistic)
    bool bconv = true;
    
    // Do contact augmentations
    // loop over all contact interfaces
    for (int i = 0; i<m_fem.SurfacePairConstraints(); ++i)
    {
        FEContactInterface* pci = dynamic_cast<FEContactInterface*>(m_fem.SurfacePairConstraint(i));
        if (pci->IsActive()) bconv = (pci->Augment(m_naug, tp) && bconv);
    }
    
    // do nonlinear constraint augmentations
    int n = m_fem.NonlinearConstraints();
    for (int i=0; i<n; ++i)
    {
        FENLConstraint* plc = m_fem.NonlinearConstraint(i);
        if (plc->IsActive()) bconv = plc->Augment(m_naug, tp) && bconv;
    }
    
    // do incompressibility multipliers for 3Field domains
    FEMesh& mesh = m_fem.GetMesh();
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
    FEModel& fem = GetFEModel();
    
    // get the time information
    FETimeInfo tp = fem.GetTime();
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
void FEFluidFSISolver::PrepStep(const FETimeInfo& timeInfo)
{
    TimerTracker t(m_UpdateTime);
    double dt = timeInfo.timeIncrement;
    
    // initialize counters
    m_niter = 0;	// nr of iterations
    m_nrhs  = 0;	// nr of RHS evaluations
    m_nref  = 0;	// nr of stiffness reformations
    m_ntotref = 0;
    m_naug  = 0;	// nr of augmentations
    
    // zero total DOFs
    zero(m_Ui);
    zero(m_Vi);
    zero(m_Di);
    zero(m_Fi);
    
    // store previous mesh state
    // we need them for strain and acceleration calculations
    FEMesh& mesh = m_fem.GetMesh();
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& ni = mesh.Node(i);
        ni.m_rp = ni.m_rt;
        ni.m_vp = ni.get_vec3d(m_dofVX, m_dofVY, m_dofVZ);
        ni.m_ap = ni.m_at;
        ni.set_vec3d(m_dofSXP, m_dofSYP, m_dofSZP, ni.get_vec3d(m_dofSX, m_dofSY, m_dofSZ));
        ni.set_vec3d(m_dofSVXP, m_dofSVYP, m_dofSVZP, ni.get_vec3d(m_dofSVX, m_dofSVY, m_dofSVZ));
        ni.set_vec3d(m_dofSAXP, m_dofSAYP, m_dofSAZP, ni.get_vec3d(m_dofSAX, m_dofSAY, m_dofSAZ));
        ni.set_vec3d(m_dofWXP, m_dofWYP, m_dofWZP, ni.get_vec3d(m_dofWX, m_dofWY, m_dofWZ));
        ni.set_vec3d(m_dofAWXP, m_dofAWYP, m_dofAWZP, ni.get_vec3d(m_dofAWX, m_dofAWY, m_dofAWZ));
        ni.set(m_dofEFP, ni.get(m_dofEF));
        ni.set(m_dofAEFP, ni.get(m_dofAEF));
        
        switch (m_pred) {
            case 0:
            {
                // initial guess at start of new time step (default)
                // solid
                ni.m_at = ni.m_ap*(1-0.5/m_beta) - ni.m_vp/(m_beta*dt);
                vec3d vs = ni.m_vp + (ni.m_at*m_gamma + ni.m_ap*(1-m_gamma))*dt;
                ni.set_vec3d(m_dofVX, m_dofVY, m_dofVZ, vs);
                
                // solid shell
                vec3d aqp = ni.get_vec3d(m_dofSAXP, m_dofSAYP, m_dofSAZP);
                vec3d vqp = ni.get_vec3d(m_dofSVXP, m_dofSVYP, m_dofSVZP);
                vec3d aqt = aqp*(1-0.5/m_beta) - vqp/(m_beta*dt);
                ni.set_vec3d(m_dofSAX, m_dofSAY, m_dofSAZ, aqt);
                vec3d vqt = vqp + (aqt*m_gamma + aqp*(1-m_gamma))*dt;
                ni.set_vec3d(m_dofSVX, m_dofSVY, m_dofSVZ, vqt);
                
                // fluid
                vec3d awp = ni.get_vec3d(m_dofAWXP, m_dofAWYP, m_dofAWZP);
                ni.set_vec3d(m_dofAWX, m_dofAWY, m_dofAWZ, awp*(m_gamma-1)/m_gamma);
                
                ni.set(m_dofAEF, ni.get(m_dofAEFP)*(m_gamma-1)/m_gamma);
            }
                break;
                
            case 1:
            {
                // initial guess at start of new time step (Zero Ydot)
                ni.m_at = vec3d(0,0,0);
                ni.set_vec3d(m_dofAWX, m_dofAWY, m_dofAWZ,vec3d(0,0,0));
                ni.set(m_dofAEF, 0);

                ni.set_vec3d(m_dofVX, m_dofVY, m_dofVZ, ni.m_vp + ni.m_ap*dt*(1-m_gamma)*m_alphaf);
                vec3d wp = ni.get_vec3d(m_dofWXP, m_dofWYP, m_dofWZP);
                vec3d awp = ni.get_vec3d(m_dofAWXP, m_dofAWYP, m_dofAWZP);
                ni.set_vec3d(m_dofWX, m_dofWY, m_dofWZ, wp + awp*dt*(1-m_gamma)*m_alphaf);
                ni.set(m_dofEF, ni.get(m_dofEFP) + ni.get(m_dofAEFP)*dt*(1-m_gamma)*m_alphaf);
            }
                break;
                
            case 2:
            {
                // initial guess at start of new time step (Same Ydot)
                vec3d awp = ni.get_vec3d(m_dofAWXP, m_dofAWYP, m_dofAWZP);
                ni.set_vec3d(m_dofAWX, m_dofAWY, m_dofAWZ, awp);
                ni.set(m_dofAEF, ni.get(m_dofAEFP));
                
                ni.set_vec3d(m_dofVX, m_dofVY, m_dofVZ, ni.m_vp + ni.m_ap*dt);
                vec3d wp = ni.get_vec3d(m_dofWXP, m_dofWYP, m_dofWZP);
                ni.set_vec3d(m_dofWX, m_dofWY, m_dofWZ, wp + awp*dt);
                ni.set(m_dofEF, ni.get(m_dofEFP) + ni.get(m_dofAEFP)*dt);
            }
                break;
                
            default:
                break;
        }
    }
    
    FETimeInfo tp = m_fem.GetTime();
    tp.alpha = m_alpha;
    tp.beta  = m_beta;
    tp.gamma = m_gamma;
    tp.alphaf = m_alphaf;
    tp.alpham = m_alpham;
    
    // apply concentrated nodal forces
    // since these forces do not depend on the geometry
    // we can do this once outside the NR loop.
    NodalForces(m_Fn, tp);
    
    // apply prescribed velocities
    // we save the prescribed velocity increments in the ui vector
    vector<double>& ui = m_ui;
    zero(ui);
    int nbc = m_fem.PrescribedBCs();
    for (int i=0; i<nbc; ++i)
    {
        FEPrescribedBC& dc = *m_fem.PrescribedBC(i);
        if (dc.IsActive()) dc.PrepStep(ui);
    }
    
    // do the linear constraints
    m_fem.GetLinearConstraintManager().PrepStep();
    
    // initialize rigid bodies
    m_rigidSolver.PrepStep(timeInfo, ui);
    
    // initialize contact
    if (m_fem.SurfacePairConstraints() > 0) UpdateContact();
    
    // initialize nonlinear constraints
    if (m_fem.NonlinearConstraints() > 0) UpdateConstraints();
    
   // initialize material point data
    // NOTE: do this before the stresses are updated
    // TODO: does it matter if the stresses are updated before
    //       the material point data is initialized
    // update domain data
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        if (dom.IsActive()) dom.PreSolveUpdate(timeInfo);
    }

    // update stresses
    UpdateStresses();
    
    // see if we need to do contact augmentations
    m_baugment = false;
    for (int i = 0; i<m_fem.SurfacePairConstraints(); ++i)
    {
        FEContactInterface& ci = dynamic_cast<FEContactInterface&>(*m_fem.SurfacePairConstraint(i));
        if (ci.IsActive() && ci.m_blaugon) m_baugment = true;
    }
    
    // see if we need to do incompressible augmentations
    int nmat = m_fem.Materials();
    for (int i = 0; i<nmat; ++i)
    {
        FEUncoupledMaterial* pmi = dynamic_cast<FEUncoupledMaterial*>(m_fem.GetMaterial(i));
        if (pmi && pmi->m_blaugon) m_baugment = true;
    }
    
    // see if we have to do nonlinear constraint augmentations
    if (m_fem.NonlinearConstraints() != 0) m_baugment = true;
}

//-----------------------------------------------------------------------------
//! Implements the BFGS2 algorithm to solve the nonlinear FE equations.
bool FEFluidFSISolver::Quasin(double time)
{
    int i;
    
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
    
    // initialize flags
    bool bconv = false;		// convergence flag
    
    static int nretries;
    
    // Get the current step
    FEAnalysis* pstep = m_fem.GetCurrentStep();
    
    // prepare for the first iteration
    FETimeInfo tp = m_fem.GetTime();
    tp.alpha = m_alpha;
    tp.beta  = m_beta;
    tp.gamma = m_gamma;
    tp.alphaf = m_alphaf;
    tp.alpham = m_alpham;
    PrepStep(tp);
    
    // calculate initial stiffness matrix
    bool breform = m_breformtimestep;
    if (pstep->m_ntotiter == 0) breform = true;
    // force reformation on a retry, if m_breformtimestep is set to false
    if ((m_breformtimestep == false) && (pstep->m_nretries > nretries)) breform = true;
    nretries = pstep->m_nretries;
    if (breform)
    {
        // reset the bfgs updates
        if (ReformStiffness(tp) == false) return false;
    }
    
    // reset reformation flag to false so that we won't reform until necessary
    breform = false;
    
    // calculate initial residual
    if (Residual(m_R0) == false) return false;
    
    // Add the "reaction forces" from prescribed dofs.
    // This vector is created by bringing the stiffness contributions
    // from the prescribed dofs to the RHS.
    m_R0 += m_Fd;
    
    // TODO: I can check here if the residual is zero.
    // If it is than there is probably no force acting on the system
    // if (m_R0*m_R0 < eps) bconv = true;
    
    //	double r0 = m_R0*m_R0;
    
    // set the initial step length estimates to 1.0
    double s = 1.0;
    
    // loop until converged or when max nr of reformations reached
    do
    {
        Logfile::MODE oldmode = felog.GetMode();
        if ((pstep->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
            (pstep->GetPrintLevel() != FE_PRINT_NEVER)) felog.SetMode(Logfile::LOG_FILE);
        
        felog.printf(" %d\n", m_niter+1);
        felog.SetMode(oldmode);
        
        // assume we'll converge.
        bconv = true;
        
        // solve the equations
        m_SolverTime.start();
        {
            m_pbfgs->SolveEquations(m_ui, m_R0);
        }
        m_SolverTime.stop();
        
        // check for nans
        m_UpdateTime.start();
        {
            double du = m_ui*m_ui;
            if (ISNAN(du)) throw NANDetected();
            
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
        }
        m_UpdateTime.stop();
        
        // perform a linesearch
        // the geometry is also updated in the line search
        if (m_LStol > 0) s = LineSearch(1.0);
        else
        {
            s = 1;
            
            // Update geometry
            Update(m_ui);
            
            // calculate residual at this point
            Residual(m_R1);
        }
        
        // calculate norms
        m_UpdateTime.start();
        {
            // update all degrees of freedom
            for (i=0; i<m_neq; ++i) m_Ui[i] += s*m_ui[i];
            
            // update displacements
            for (i=0; i<m_ndeq; ++i) m_Di[i] += s*m_di[i];
            
            // update velocities
            for (i=0; i<m_nveq; ++i) m_Vi[i] += s*m_vi[i];
            
            // update dilatations
            for (i=0; i<m_nfeq; ++i) m_Fi[i] += s*m_fi[i];
            
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
        }
        m_UpdateTime.stop();
        
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
        if ((m_LStol > 0) && (s < m_LSmin)) bconv = false;
        
        // check energy divergence
        if (normE1 > normEm) bconv = false;
        
        // print convergence summary
        oldmode = felog.GetMode();
        if ((pstep->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
            (pstep->GetPrintLevel() != FE_PRINT_NEVER)) felog.SetMode(Logfile::LOG_FILE);
        
        felog.printf(" Nonlinear solution status: time= %lg\n", time);
        felog.printf("\tstiffness updates             = %d\n", m_pbfgs->m_nups);
        felog.printf("\tright hand side evaluations   = %d\n", m_nrhs);
        felog.printf("\tstiffness matrix reformations = %d\n", m_nref);
        if (m_LStol > 0) felog.printf("\tstep from line search         = %lf\n", s);
        felog.printf("\tconvergence norms :     INITIAL         CURRENT         REQUIRED\n");
        felog.printf("\t   residual         %15le %15le %15le \n", normRi, normR1, m_Rtol*normRi);
        felog.printf("\t   energy           %15le %15le %15le \n", normEi, normE1, m_Etol*normEi);
        felog.printf("\t   displacement     %15le %15le %15le \n", normDi, normd ,(m_Dtol*m_Dtol)*normD );
        felog.printf("\t   velocity         %15le %15le %15le \n", normVi, normv ,(m_Vtol*m_Vtol)*normV );
        felog.printf("\t   dilatation       %15le %15le %15le \n", normFi, normf ,(m_Ftol*m_Ftol)*normF );
        
        felog.SetMode(oldmode);
        
        // see if we may have a small residual
        if ((bconv == false) && (normR1 < m_Rmin))
        {
            // check for almost zero-residual on the first iteration
            // this might be an indication that there is no force on the system
            felog.printbox("WARNING", "No force acting on the system.");
            bconv = true;
        }
        
        // check if we have converged.
        // If not, calculate the BFGS update vectors
        if (bconv == false)
        {
            if (s < m_LSmin)
            {
                // check for zero linestep size
                felog.printbox("WARNING", "Zero linestep size. Stiffness matrix will now be reformed");
                breform = true;
            }
            else if ((normE1 > normEm) && m_bdivreform)
            {
                // check for diverging
                felog.printbox("WARNING", "Problem is diverging. Stiffness matrix will now be reformed");
                normEm = normE1;
                normEi = normE1;
                normRi = normR1;
                normDi = normd;
                normVi = normv;
                normFi = normf;
                breform = true;
            }
            else
            {
                // If we havn't reached max nr of BFGS updates
                // do an update
                if (!breform)
                {
                    if (m_pbfgs->m_nups < m_pbfgs->m_maxups-1)
                    {
                        m_QNTime.start();
                        if (m_pbfgs->Update(s, m_ui, m_R0, m_R1) == false)
                        {
                            // Stiffness update has failed.
                            // this might be due a too large condition number
                            // or the update was no longer positive definite.
                            felog.printbox("WARNING", "The BFGS update has failed.\nStiffness matrix will now be reformed.");
                            breform = true;
                        }
                        m_QNTime.stop();
                    }
                    else
                    {
                        // we've reached the max nr of BFGS updates, so
                        // we need to do a stiffness reformation
                        breform = true;
                        
                        // print a warning only if the user did not intent full-Newton
                        if (m_pbfgs->m_maxups > 0)
                            felog.printbox("WARNING", "Max nr of iterations reached.\nStiffness matrix will now be reformed.");
                        
                    }
                }
            }
            
            // zero velocity increments
            // we must set this to zero before the reformation
            // because we assume that the prescribed velocities are stored
            // in the m_ui vector.
            m_UpdateTime.start();
            zero(m_ui);
            m_UpdateTime.stop();
            
            // reform stiffness matrices if necessary
            if (breform && m_bdoreforms)
            {
                felog.printf("Reforming stiffness matrix: reformation #%d\n\n", m_nref);
                
                // reform the matrix
                if (ReformStiffness(tp) == false) break;
                
                // reset reformation flag
                breform = false;
            }
            
            // copy last calculated residual
            m_RHSTime.start();
            m_R0 = m_R1;
            m_RHSTime.stop();
        }
        else if (m_baugment)
        {
            // we have converged, so let's see if the augmentations have converged as well
            felog.printf("\n........................ augmentation # %d\n", m_naug+1);
            
            // plot states before augmentations.
            // The reason we store the state prior to the augmentations
            // is because the augmentations are going to change things such that
            // the system no longer in equilibrium. Since the model has to be converged
            // before we do augmentations, storing the model now will store an actual converged state.
            pstep->GetFEModel().DoCallback(CB_AUGMENT);
            
            // do the augmentations
            bconv = Augment();
            
            // update counter
            ++m_naug;
            
            // we reset the reformations counter
            m_nref = 0;
            
            // If we havn't converged we prepare for the next iteration
            if (!bconv)
            {
                // Since the Lagrange multipliers have changed, we can't just copy
                // the last residual but have to recalculate the residual
                // we also recalculate the stresses in case we are doing augmentations
                // for incompressible materials
                UpdateStresses();
                Residual(m_R0);
                
                // reform the matrix if we are using full-Newton
                if (m_pbfgs->m_maxups == 0)
                {
                    felog.printf("Reforming stiffness matrix: reformation #%d\n\n", m_nref);
                    if (ReformStiffness(tp) == false) break;
                }
            }
        }
        
        // increase iteration number
        m_niter++;
        
        // let's flush the logfile to make sure the last output will not get lost
        felog.flush();
        
        // do minor iterations callbacks
        m_fem.DoCallback(CB_MINOR_ITERS);
    }
    while (bconv == false);
    
    // when converged,
    // print a convergence summary to the felog file
    if (bconv)
    {
        Logfile::MODE mode = felog.GetMode();
        if (mode != Logfile::LOG_NEVER)
        {
            felog.SetMode(Logfile::LOG_FILE);
            felog.printf("\nconvergence summary\n");
            felog.printf("    number of iterations   : %d\n", m_niter);
            felog.printf("    number of reformations : %d\n", m_nref);
            felog.SetMode(mode);
        }
    }
    
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

bool FEFluidFSISolver::StiffnessMatrix(const FETimeInfo& tp)
{
    // get the stiffness matrix
    SparseMatrix& K = *m_pK;
    
    // zero stiffness matrix
    K.zero();
    
    // zero the residual adjustment vector
    zero(m_Fd);
    
    // nodal degrees of freedom
    int i;
    
    // get the mesh
    FEMesh& mesh = m_fem.GetMesh();
    
    // calculate the stiffness matrix for each domain
    for (i=0; i<mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        if (dom.IsActive()) {
            FEFluidDomain* fdom = dynamic_cast<FEFluidDomain*>(&dom);
            FEFluidFSIDomain* fsidom = dynamic_cast<FEFluidFSIDomain*>(&dom);
            FEElasticDomain* edom = dynamic_cast<FEElasticDomain*>(&dom);
            if (fdom) fdom->StiffnessMatrix(this, tp);
            else if (fsidom) fsidom->StiffnessMatrix(this, tp);
            else if (edom) edom->StiffnessMatrix(this);
        }
    }
    
    // calculate the body force stiffness matrix for each domain
    // but not for solid domains (since they have no mass in FSI)
    for (i=0; i<mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        if (dom.IsActive() && (dom.GetMaterial()->IsRigid() == false)) {
            FEFluidDomain* fdom = dynamic_cast<FEFluidDomain*>(&dom);
            FEFluidFSIDomain* fsidom = dynamic_cast<FEFluidFSIDomain*>(&dom);
            FEElasticDomain* edom = dynamic_cast<FEElasticDomain*>(&dom);
            int NBL = m_fem.BodyLoads();
            for (int j=0; j<NBL; ++j)
            {
                FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(j));
                if (pbf) {
                    if (fdom) fdom->BodyForceStiffness(this, tp, *pbf);
                    else if (fsidom) fsidom->BodyForceStiffness(this, tp, *pbf);
                    else if (edom) edom->BodyForceStiffness(this, *pbf);
                }
            }
        }
    }
    
    // TODO: add body force stiffness for rigid bodies
    
    // Add mass matrix
    FEAnalysis* pstep = m_fem.GetCurrentStep();
    if (pstep->m_nanalysis == FE_DYNAMIC)
    {
        // scale factor
        double dt = tp.timeIncrement;
        double a = tp.alpham / (m_beta*dt*dt);
        // loop over all domains (except rigid)
        for (i=0; i<mesh.Domains(); ++i)
        {
            FEDomain& dom = mesh.Domain(i);
            if (dom.IsActive() && dom.GetMaterial()->IsRigid() == false) {
                FEFluidDomain* fdom = dynamic_cast<FEFluidDomain*>(&dom);
                FEFluidFSIDomain* fsidom = dynamic_cast<FEFluidFSIDomain*>(&dom);
                FEElasticDomain* edom = dynamic_cast<FEElasticDomain*>(&dom);
                if (fdom) fdom->MassMatrix(this, tp);
                else if (fsidom) fsidom->MassMatrix(this, tp);
                else if (edom) edom->MassMatrix(this, a);
            }
        }
        
        m_rigidSolver.RigidMassMatrix(this, tp);
    }
    
    // calculate contact stiffness
    ContactStiffness();

    // calculate stiffness matrix due to surface loads
    int nsl = m_fem.SurfaceLoads();
    for (int i=0; i<nsl; ++i)
    {
        FESurfaceLoad* psl = m_fem.SurfaceLoad(i);
        if (psl->IsActive()) psl->StiffnessMatrix(tp, this);
    }
    
    // calculate nonlinear constraint stiffness
    // note that this is the contribution of the
    // constrainst enforced with augmented lagrangian
    NonLinearConstraintStiffness(tp);
    
    // calculate the stiffness contributions for the rigid forces
    for (int i = 0; i<m_fem.ModelLoads(); ++i) m_fem.ModelLoad(i)->StiffnessMatrix(this, tp);
    
    // add contributions from rigid bodies
    m_rigidSolver.StiffnessMatrix(*m_pK, tp);
    
    return true;
}

//-----------------------------------------------------------------------------
//! Calculate the stiffness contribution due to nonlinear constraints
void FEFluidFSISolver::NonLinearConstraintStiffness(const FETimeInfo& tp)
{
    int N = m_fem.NonlinearConstraints();
    for (int i=0; i<N; ++i)
    {
        FENLConstraint* plc = m_fem.NonlinearConstraint(i);
        if (plc->IsActive()) plc->StiffnessMatrix(this, tp);
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the contact stiffness matrix

void FEFluidFSISolver::ContactStiffness()
{
    FETimeInfo tp = GetFEModel().GetTime();
    tp.alpha = m_alpha;
    tp.beta  = m_beta;
    tp.gamma = m_gamma;
    tp.alphaf = m_alphaf;
    tp.alpham = m_alpham;
    for (int i = 0; i<m_fem.SurfacePairConstraints(); ++i)
    {
        FEContactInterface* pci = dynamic_cast<FEContactInterface*>(m_fem.SurfacePairConstraint(i));
        if (pci->IsActive()) pci->StiffnessMatrix(this, tp);
    }
}

//-----------------------------------------------------------------------------
void FEFluidFSISolver::AssembleResidual(int node_id, int dof, double f, vector<double>& R)
{
    // get the mesh
    FEMesh& mesh = m_fem.GetMesh();
    
    // get the equation number
    FENode& node = mesh.Node(node_id);
    int n = node.m_ID[dof];
    
    // assemble into global vector
    if (n >= 0) R[n] += f;
    else m_rigidSolver.AssembleResidual(node_id, dof, f, R);
}

//-----------------------------------------------------------------------------
//! \todo This function is only used for rigid joints. I need to figure out if
//!       I can use the other assembly function.
void FEFluidFSISolver::AssembleStiffness(std::vector<int>& lm, matrix& ke)
{
    m_pK->Assemble(ke, lm);
}

//-----------------------------------------------------------------------------
//!  Assembles the element stiffness matrix into the global stiffness matrix.
//!  Also adjusts the global stiffness matrix and residual to take the
//!  prescribed velocities into account.

//! \todo In stead of changing the global stiffness matrix to accomodate for
//!       the rigid bodies and linear constraints, can I modify the element stiffness
//!       matrix prior to assembly? I might have to change the elm vector as well as
//!       the element matrix size.

void FEFluidFSISolver::AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke)
{
    // assemble into global stiffness matrix
    m_pK->Assemble(ke, elm);
    
    vector<double>& ui = m_ui;
    
    // adjust for linear constraints
    FELinearConstraintManager& LCM = m_fem.GetLinearConstraintManager();
    if (LCM.LinearConstraints() > 0)
    {
        LCM.AssembleStiffness(*m_pK, m_Fd, m_ui, en, elm, ke);
    }
    
    // adjust stiffness matrix for prescribed degrees of freedom
    // NOTE: I had to comment this if statement out since otherwise
    //       poroelastic DOF's that are set as free-draining in the
    //       sliding2 contact code are skipt and zeroes will appear
    //       on the diagonal of the stiffness matrix.
    //	if (m_fem.m_DC.size() > 0)
    {
        int i, j;
        int I, J;
        
        SparseMatrix& K = *m_pK;
        
        int N = ke.rows();
        
        // loop over columns
        for (j=0; j<N; ++j)
        {
            J = -elm[j]-2;
            if ((J >= 0) && (J<m_nreq))
            {
                // dof j is a prescribed degree of freedom
                
                // loop over rows
                for (i=0; i<N; ++i)
                {
                    I = elm[i];
                    if (I >= 0)
                    {
                        // dof i is not a prescribed degree of freedom
                        m_Fd[I] -= ke[i][j]*ui[J];
                    }
                }
                
                // set the diagonal element of K to 1
                K.set(J,J, 1);
            }
        }
    }
    
    // see if there are any rigid body dofs here
    m_rigidSolver.RigidStiffness(*m_pK, m_ui, m_Fd, en, elm, ke, m_alpha);
}

//-----------------------------------------------------------------------------
// \todo adjust for rigid bodies
void FEFluidFSISolver::AssembleStiffness2(vector<int>& lmi, vector<int>& lmj, matrix& ke)
{
    m_pK->Assemble(ke, lmi, lmj);
    
    // adjust for prescribed dofs
    vector<double>& ui = m_ui;
    
    SparseMatrix& K = *m_pK;
    // loop over columns
    int cols = ke.columns();
    int rows = ke.rows();
    for (int j=0; j<cols; ++j)
    {
        int J = -lmj[j] - 2;
        if ((J >= 0) && (J<m_nreq))
        {
            // dof j is a prescribed degree of freedom
            
            // loop over rows
            for (int i=0; i<rows; ++i)
            {
                int I = lmi[i];
                if (I >= 0)
                {
                    // dof i is not a prescribed degree of freedom
                    m_Fd[I] -= ke[i][j]*ui[J];
                }
            }
            
            // set the diagonal element of K to 1
            K.set(J,J, 1);
        }
    }
    
    // TODO: implement this in the FERigidSolver
}

//-----------------------------------------------------------------------------
//! Calculates the contact forces
void FEFluidFSISolver::ContactForces(FEGlobalVector& R)
{
    FETimeInfo tp = GetFEModel().GetTime();
    tp.alpha = m_alpha;
    tp.beta  = m_beta;
    tp.gamma = m_gamma;
    tp.alphaf = m_alphaf;
    tp.alpham = m_alpham;
    for (int i = 0; i<m_fem.SurfacePairConstraints(); ++i)
    {
        FEContactInterface* pci = dynamic_cast<FEContactInterface*>(m_fem.SurfacePairConstraint(i));
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
    TimerTracker t(m_RHSTime);
    
    // get the time information
    FETimeInfo tp = m_fem.GetTime();
    tp.alpha = m_alpha;
    tp.beta  = m_beta;
    tp.gamma = m_gamma;
    tp.alphaf = m_alphaf;
    tp.alpham = m_alpham;
    
    // initialize residual with concentrated nodal loads
    R = m_Fn;
    
    // zero nodal reaction forces
    zero(m_Fr);
    
    // setup the global vector
    FEResidualVector RHS(GetFEModel(), R, m_Fr);

    // zero rigid body reaction forces
    m_rigidSolver.Residual();
    
    // get the mesh
    FEMesh& mesh = m_fem.GetMesh();
    
    // set flag for transient or steady-state analyses
    FEAnalysis* pstep = m_fem.GetCurrentStep();
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        if (dom.IsActive() && dom.GetMaterial()->IsRigid() == false) {
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
        if (dom.IsActive() && dom.GetMaterial()->IsRigid() == false) {
            FEFluidDomain* fdom = dynamic_cast<FEFluidDomain*>(&dom);
            FEFluidFSIDomain* fsidom = dynamic_cast<FEFluidFSIDomain*>(&dom);
            FEElasticDomain* edom = dynamic_cast<FEElasticDomain*>(&dom);
            if (fdom) fdom->InternalForces(RHS, tp);
            else if (fsidom) fsidom->InternalForces(RHS, tp);
            else if (edom) edom->InternalForces(RHS);
        }
    }
    
    // calculate the body forces
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        if (dom.IsActive() && dom.GetMaterial()->IsRigid() == false) {
            FEFluidDomain* fdom = dynamic_cast<FEFluidDomain*>(&dom);
            FEFluidFSIDomain* fsidom = dynamic_cast<FEFluidFSIDomain*>(&dom);
            FEElasticDomain* edom = dynamic_cast<FEElasticDomain*>(&dom);
            for (int j=0; j<m_fem.BodyLoads(); ++j)
            {
                FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(j));
                if (fdom) fdom->BodyForce(RHS, tp, *pbf);
                else if (fsidom) fsidom->BodyForce(RHS, tp, *pbf);
                else if (edom) edom->BodyForce(RHS, *pbf);
            }
        }
    }
    
    // calculate body forces for rigid bodies
    for (int j=0; j<m_fem.BodyLoads(); ++j)
    {
        FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(j));
        m_rigidSolver.BodyForces(RHS, tp, *pbf);
    }

    // allocate F
    vector<double> F;
    
    // calculate inertial forces
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        if (dom.IsActive() && dom.GetMaterial()->IsRigid() == false) {
            FEFluidDomain* fdom = dynamic_cast<FEFluidDomain*>(&dom);
            FEFluidFSIDomain* fsidom = dynamic_cast<FEFluidFSIDomain*>(&dom);
            FEElasticDomain* edom = dynamic_cast<FEElasticDomain*>(&dom);
            if (fdom) fdom->InertialForces(RHS, tp);
            else if (fsidom) fsidom->InertialForces(RHS, tp);
            else if (edom && (pstep->m_nanalysis == FE_DYNAMIC)) edom->InertialForces(RHS, F);
        }
    }
    
    // update rigid bodies
    if (pstep->m_nanalysis == FE_DYNAMIC) m_rigidSolver.InertialForces(RHS, tp);

    // calculate forces due to surface loads
    int nsl = m_fem.SurfaceLoads();
    for (int i=0; i<nsl; ++i)
    {
        FESurfaceLoad* psl = m_fem.SurfaceLoad(i);
        if (psl->IsActive()) psl->Residual(tp, RHS);
    }
    
    // calculate contact forces
    ContactForces(RHS);
    
    // calculate nonlinear constraint forces
    // note that these are the linear constraints
    // enforced using the augmented lagrangian
    NonLinearConstraintForces(RHS, tp);
    
    // add model loads
    int NML = m_fem.ModelLoads();
    for (int i=0; i<NML; ++i)
    {
        FEModelLoad& mli = *m_fem.ModelLoad(i);
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
        if ((n = -node.m_ID[m_dofX]-2) >= 0) node.m_Fr.x = -m_Fr[n];
        if ((n = -node.m_ID[m_dofY]-2) >= 0) node.m_Fr.y = -m_Fr[n];
        if ((n = -node.m_ID[m_dofZ]-2) >= 0) node.m_Fr.z = -m_Fr[n];
    }
    
    // increase RHS counter
    m_nrhs++;
    
    return true;
}

//-----------------------------------------------------------------------------
//! calculate the nonlinear constraint forces
void FEFluidFSISolver::NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp)
{
    int N = m_fem.NonlinearConstraints();
    for (int i=0; i<N; ++i)
    {
        FENLConstraint* plc = m_fem.NonlinearConstraint(i);
        if (plc->IsActive()) plc->Residual(R, tp);
    }
}

//-----------------------------------------------------------------------------
//! calculates the concentrated nodal forces

void FEFluidFSISolver::NodalForces(vector<double>& F, const FETimeInfo& tp)
{
    // zero nodal force vector
    zero(F);
    
    // loop over nodal loads
    int NNL = m_fem.NodalLoads();
    for (int i=0; i<NNL; ++i)
    {
        const FENodalLoad& fc = *m_fem.NodalLoad(i);
        if (fc.IsActive())
        {
            int dof = fc.GetDOF();
            int N = fc.Nodes();
            for (int j=0; j<N; ++j)
            {
                int nid = fc.NodeID(j);
                
                // get the nodal load value
                double f = fc.NodeValue(j);
                
                // assemble into residual
                AssembleResidual(nid, dof, f, F);
            }
        }
    }
}
