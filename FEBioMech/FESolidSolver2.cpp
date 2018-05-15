#include "stdafx.h"
#include "FESolidSolver2.h"
#include "FERigidMaterial.h"
#include "FERigidConnector.h"
#include "FESlidingInterfaceBW.h"
#include "FE3FieldElasticSolidDomain.h"
#include "FE3FieldElasticShellDomain.h"
#include "FEBodyForce.h"
#include "FEPressureLoad.h"
#include "FEResidualVector.h"
#include "FECore/FENodeReorder.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include "FEUncoupledMaterial.h"
#include "NumCore/NumCore.h"
#include "FECore/FEGlobalMatrix.h"
#include "FEContactInterface.h"
#include <FECore/sys.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/BC.h>
#include <FECore/RigidBC.h>
#include <FECore/FEModelLoad.h>
#include <FECore/FELinearConstraintManager.h>
#include "FESSIShellDomain.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_PARAMETER_LIST(FESolidSolver2, FENewtonSolver)
	ADD_PARAMETER2(m_Dtol        , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "dtol"        );
	ADD_PARAMETER2(m_Etol        , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "etol"        );
	ADD_PARAMETER2(m_Rtol        , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "rtol"        );
	ADD_PARAMETER2(m_Rmin        , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "min_residual");
	ADD_PARAMETER2(m_Rmax        , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "max_residual");
    ADD_PARAMETER(m_rhoi         , FE_PARAM_DOUBLE, "rhoi"        );
    ADD_PARAMETER(m_alpha        , FE_PARAM_DOUBLE, "alpha"       );
	ADD_PARAMETER(m_beta         , FE_PARAM_DOUBLE, "beta"        );
	ADD_PARAMETER(m_gamma        , FE_PARAM_DOUBLE, "gamma"       );
    ADD_PARAMETER(m_bsymm        , FE_PARAM_BOOL  , "symmetric_stiffness");
	ADD_PARAMETER(m_logSolve     , FE_PARAM_BOOL  ,"logSolve");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! FESolidSolver2 Construction
//
FESolidSolver2::FESolidSolver2(FEModel* pfem) : FENewtonSolver(pfem), m_rigidSolver(pfem)
{
	// default values
	m_Rtol = 0;	// deactivate residual convergence 
	m_Dtol = 0.001;
	m_Etol = 0.01;
	m_Rmin = 1.0e-20;
	m_Rmax = 0;	// not used if zero

	m_niter = 0;
	m_nreq = 0;

	m_logSolve = false;

	// default Newmark parameters (trapezoidal rule)
    m_rhoi = -2;
    m_alpha = m_alphaf = 1.0;
    m_alpham = 1.0;
	m_beta  = 0.25;
	m_gamma = 0.5;

	// Allocate degrees of freedom
	DOFS& dofs = pfem->GetDOFS();
	int varD = dofs.AddVariable("displacement", VAR_VEC3);
	dofs.SetDOFName(varD, 0, "x");
	dofs.SetDOFName(varD, 1, "y");
	dofs.SetDOFName(varD, 2, "z");
	int varQ = dofs.AddVariable("rotation", VAR_VEC3);
	dofs.SetDOFName(varQ, 0, "u");
	dofs.SetDOFName(varQ, 1, "v");
	dofs.SetDOFName(varQ, 2, "w");
	int varQR = dofs.AddVariable("rigid rotation", VAR_VEC3);
	dofs.SetDOFName(varQR, 0, "Ru");
	dofs.SetDOFName(varQR, 1, "Rv");
	dofs.SetDOFName(varQR, 2, "Rw");
    int varSD = dofs.AddVariable("shell displacement", VAR_VEC3);
    dofs.SetDOFName(varSD, 0, "sx");
    dofs.SetDOFName(varSD, 1, "sy");
    dofs.SetDOFName(varSD, 2, "sz");
	int varV = dofs.AddVariable("velocity", VAR_VEC3);
	dofs.SetDOFName(varV, 0, "vx");
	dofs.SetDOFName(varV, 1, "vy");
	dofs.SetDOFName(varV, 2, "vz");
    int varQP = dofs.AddVariable("previous rotation", VAR_VEC3);
    dofs.SetDOFName(varQP, 0, "up");
    dofs.SetDOFName(varQP, 1, "vp");
    dofs.SetDOFName(varQP, 2, "wp");
    int varSDP = dofs.AddVariable("previous shell displacement", VAR_VEC3);
    dofs.SetDOFName(varSDP, 0, "sxp");
    dofs.SetDOFName(varSDP, 1, "syp");
    dofs.SetDOFName(varSDP, 2, "szp");
    int varQV = dofs.AddVariable("shell velocity", VAR_VEC3);
    dofs.SetDOFName(varQV, 0, "svx");
    dofs.SetDOFName(varQV, 1, "svy");
    dofs.SetDOFName(varQV, 2, "svz");
    int varQA = dofs.AddVariable("shell acceleration", VAR_VEC3);
    dofs.SetDOFName(varQA, 0, "sax");
    dofs.SetDOFName(varQA, 1, "say");
    dofs.SetDOFName(varQA, 2, "saz");
    int varQVP = dofs.AddVariable("previous shell velocity", VAR_VEC3);
    dofs.SetDOFName(varQVP, 0, "svxp");
    dofs.SetDOFName(varQVP, 1, "svyp");
    dofs.SetDOFName(varQVP, 2, "svzp");
    int varQAP = dofs.AddVariable("previous shell acceleration", VAR_VEC3);
    dofs.SetDOFName(varQAP, 0, "saxp");
    dofs.SetDOFName(varQAP, 1, "sayp");
    dofs.SetDOFName(varQAP, 2, "sazp");
    
    // get the DOF indices
	m_dofX  = pfem->GetDOFIndex("x");
	m_dofY  = pfem->GetDOFIndex("y");
	m_dofZ  = pfem->GetDOFIndex("z");
	m_dofVX = pfem->GetDOFIndex("vx");
	m_dofVY = pfem->GetDOFIndex("vy");
	m_dofVZ = pfem->GetDOFIndex("vz");
	m_dofU  = pfem->GetDOFIndex("u");
	m_dofV  = pfem->GetDOFIndex("v");
	m_dofW  = pfem->GetDOFIndex("w");
	m_dofRU = pfem->GetDOFIndex("Ru");
	m_dofRV = pfem->GetDOFIndex("Rv");
	m_dofRW = pfem->GetDOFIndex("Rw");
    
    m_dofSX  = pfem->GetDOFIndex("sx");
    m_dofSY  = pfem->GetDOFIndex("sy");
    m_dofSZ  = pfem->GetDOFIndex("sz");
    m_dofSVX  = pfem->GetDOFIndex("svx");
    m_dofSVY  = pfem->GetDOFIndex("svy");
    m_dofSVZ  = pfem->GetDOFIndex("svz");
    m_dofSAX  = pfem->GetDOFIndex("sax");
    m_dofSAY  = pfem->GetDOFIndex("say");
    m_dofSAZ  = pfem->GetDOFIndex("saz");
    
    m_dofSXP  = pfem->GetDOFIndex("sxp");
    m_dofSYP  = pfem->GetDOFIndex("syp");
    m_dofSZP  = pfem->GetDOFIndex("szp");
    m_dofSVXP  = pfem->GetDOFIndex("svxp");
    m_dofSVYP  = pfem->GetDOFIndex("svyp");
    m_dofSVZP  = pfem->GetDOFIndex("svzp");
    m_dofSAXP  = pfem->GetDOFIndex("saxp");
    m_dofSAYP  = pfem->GetDOFIndex("sayp");
    m_dofSAZP  = pfem->GetDOFIndex("sazp");
}

//-----------------------------------------------------------------------------
FESolidSolver2::~FESolidSolver2()
{

}

//-----------------------------------------------------------------------------
//! Generate warnings if needed
void FESolidSolver2:: SolverWarnings()
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
//! Allocates and initializes the data structures used by the FESolidSolver2
//
bool FESolidSolver2::Init()
{
	// initialize base class
	if (FENewtonSolver::Init() == false) return false;


    if (m_rhoi == -1) {
        // Euler integration
        m_alpha = m_alphaf = m_alpham = 1.0;
        m_beta = pow(1 + m_alpham - m_alphaf,2)/4;
        m_gamma = 0.5 + m_alpham - m_alphaf;
    }
    else if ((m_rhoi >= 0) && (m_rhoi <= 1)) {
        // Generalized-alpha integration (2nd order system)
        m_alpha = m_alphaf = 1.0/(1+m_rhoi);
        m_alpham = (2-m_rhoi)/(1+m_rhoi);
        m_beta = pow(1 + m_alpham - m_alphaf,2)/4;
        m_gamma = 0.5 + m_alpham - m_alphaf;
    }
    else {
        // for any other value of rhoi, use the user-defined alpha, beta, gamma parameters
        m_alphaf = m_alpham = m_alpha;
    }
    
	// allocate vectors
	int neq = m_neq;
	m_Fn.assign(neq, 0);
	m_Fr.assign(neq, 0);
	m_Ui.assign(neq, 0);
	m_Ut.assign(neq, 0);

	// we need to fill the total displacement vector m_Ut
	FEMesh& mesh = m_fem.GetMesh();
	gather(m_Ut, mesh, m_dofX);
	gather(m_Ut, mesh, m_dofY);
	gather(m_Ut, mesh, m_dofZ);
	gather(m_Ut, mesh, m_dofU);
	gather(m_Ut, mesh, m_dofV);
	gather(m_Ut, mesh, m_dofW);
    gather(m_Ut, mesh, m_dofSX);
    gather(m_Ut, mesh, m_dofSY);
    gather(m_Ut, mesh, m_dofSZ);

    SolverWarnings();
    
	return true;
}

//-----------------------------------------------------------------------------
//! Save data to dump file

void FESolidSolver2::Serialize(DumpStream& ar)
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
	}
	else
	{
		ar >> m_nrhs;
		ar >> m_niter;
		ar >> m_nref >> m_ntotref;
		ar >> m_naug;
		ar >> m_nreq;

		// initialize data structures
		// (only when number of equations is non-zero.
		// This can be zero in a multi-step analysis for steps that have not yet been initialized.)
		if (m_neq > 0) Init();
	}

	// serialize rigid solver
	m_rigidSolver.Serialize(ar);
}

//-----------------------------------------------------------------------------
bool FESolidSolver2::InitEquations()
{
	// First call the base class.
	// This will initialize all equation numbers, except the rigid body equation numbers
	if (FENewtonSolver::InitEquations() == false) return false;

	// store the number of equations we currently have
	m_nreq = m_neq;

	// Next, we assign equation numbers to the rigid body degrees of freedom
	int neq = m_rigidSolver.InitEquations(m_neq);
	if (neq == -1) return false; 
	else m_neq = neq;

	// All initialization is done
	return true;
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
bool FESolidSolver2::Augment()
{
	const FETimeInfo& tp = m_fem.GetTime();

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
        FE3FieldElasticShellDomain* ps = dynamic_cast<FE3FieldElasticShellDomain*>(&mesh.Domain(i));
		if (pd && pd->IsActive()) bconv = (pd->Augment(m_naug) && bconv);
        else if (ps && ps->IsActive()) bconv = (ps->Augment(m_naug) && bconv);
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! Update the kinematics of the model, such as nodal positions, velocities,
//! accelerations, etc.
void FESolidSolver2::UpdateKinematics(vector<double>& ui)
{
	// get the mesh
	FEMesh& mesh = m_fem.GetMesh();

	// update rigid bodies
	m_rigidSolver.UpdateRigidBodies(m_Ui, ui);

	// total displacements
	vector<double> U(m_Ut.size());
	for (size_t i=0; i<m_Ut.size(); ++i) U[i] = ui[i] + m_Ui[i] + m_Ut[i];

	// update flexible nodes
	// translational dofs
	scatter(U, mesh, m_dofX);
	scatter(U, mesh, m_dofY);
	scatter(U, mesh, m_dofZ);
	// rotational dofs
	scatter(U, mesh, m_dofU);
	scatter(U, mesh, m_dofV);
	scatter(U, mesh, m_dofW);
    // shell dofs
    scatter(U, mesh, m_dofSX);
    scatter(U, mesh, m_dofSY);
    scatter(U, mesh, m_dofSZ);

	// make sure the prescribed displacements are fullfilled
	int ndis = m_fem.PrescribedBCs();
	for (int i = 0; i<ndis; ++i)
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

	// update velocity and accelerations
	// for dynamic simulations
	FEAnalysis* pstep = m_fem.GetCurrentStep();
	if (pstep->m_nanalysis == FE_DYNAMIC)
	{
		int N = mesh.Nodes();
		double dt = m_fem.GetTime().timeIncrement;
		double a = 1.0 / (m_beta*dt);
		double b = a / dt;
		double c = 1.0 - 0.5/m_beta;
		for (int i=0; i<N; ++i)
		{
			FENode& n = mesh.Node(i);
			n.m_at = (n.m_rt - n.m_rp)*b - n.m_vp*a + n.m_ap*c;
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
        }
    }
}

//-----------------------------------------------------------------------------
//! Update DOF increments
void FESolidSolver2::UpdateIncrements(vector<double>& Ui, vector<double>& ui, bool emap)
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
        if ((n = node.m_ID[m_dofU]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofV]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofW]) >= 0) Ui[n] += ui[n];
        
        // shell dofs
        if ((n = node.m_ID[m_dofSX]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofSY]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofSZ]) >= 0) Ui[n] += ui[n];
	}
}

//-----------------------------------------------------------------------------
//! Updates the current state of the model
void FESolidSolver2::Update(vector<double>& ui)
{
	TRACK_TIME("update");

    // update EAS
    UpdateEAS(ui);
    UpdateIncrementsEAS(ui, true);

	// update kinematics
	UpdateKinematics(ui);

	// update contact
	if (m_fem.SurfacePairConstraints() > 0) UpdateContact();

	// update constraints
	if (m_fem.NonlinearConstraints() > 0) UpdateConstraints();

	// update element stresses
	UpdateModel();

	// update other stuff that may depend on the deformation
	int NBL = m_fem.BodyLoads();
	for (int i=0; i<NBL; ++i)
	{
		FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(i));
		if (pbf && pbf->IsActive()) pbf->Update();
	}
}

//-----------------------------------------------------------------------------
//! Update EAS
void FESolidSolver2::UpdateEAS(vector<double>& ui)
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
void FESolidSolver2::UpdateIncrementsEAS(vector<double>& ui, const bool binc)
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
void FESolidSolver2::UpdateModel()
{
	FEMesh& mesh = m_fem.GetMesh();
	const FETimeInfo& tp = m_fem.GetTime();

    // update the stresses on all domains
	for (int i=0; i<mesh.Domains(); ++i) 
	{
		FEDomain& dom = mesh.Domain(i);
		if (dom.IsActive()) dom.Update(tp);
	}
}

//-----------------------------------------------------------------------------
//! Update contact interfaces.
void FESolidSolver2::UpdateContact()
{
	// Update all contact interfaces
	const FETimeInfo& tp = m_fem.GetTime();
	for (int i = 0; i<m_fem.SurfacePairConstraints(); ++i)
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(m_fem.SurfacePairConstraint(i));
		if (pci->IsActive()) pci->Update(m_niter, tp);
	}
}

//-----------------------------------------------------------------------------
//! Update nonlinear constraints
void FESolidSolver2::UpdateConstraints()
{
	FETimeInfo& tp = m_fem.GetTime();
	tp.currentIteration = m_niter;

	// Update all nonlinear constraints
	for (int i=0; i<m_fem.NonlinearConstraints(); ++i) 
	{
		FENLConstraint* pci = m_fem.NonlinearConstraint(i);
		if (pci->IsActive()) pci->Update(m_niter, tp);
	}
}

//-----------------------------------------------------------------------------
bool FESolidSolver2::InitStep(double time)
{
	FEModel& fem = GetFEModel();

	// set time integration parameters
	FETimeInfo& tp = m_fem.GetTime();
	tp.alpha = m_alpha;
	tp.beta = m_beta;
	tp.gamma = m_gamma;
	tp.alphaf = m_alphaf;
	tp.alpham = m_alpham;

    // evaluate load curve values at current (or intermediate) time
	double t = tp.currentTime;
//	double dt = tp.timeIncrement;
//	double ta = (t > 0) ? t - (1-m_alpha)*dt : m_alpha*dt;
//	return FESolver::InitStep(ta);
    return FESolver::InitStep(t);
}

//-----------------------------------------------------------------------------
//! Prepares the data for the first BFGS-iteration. 
void FESolidSolver2::PrepStep()
{
	TRACK_TIME("update");

    const FETimeInfo& tp = m_fem.GetTime();
	double dt = tp.timeIncrement;
    
	// zero total displacements
	zero(m_Ui);

	// store previous mesh state
	// we need them for velocity and acceleration calculations
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

        // initial guess at start of new time step
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
    }

    // apply concentrated nodal forces
	// since these forces do not depend on the geometry
	// we can do this once outside the NR loop.
	NodalForces(m_Fn, tp);

	// apply prescribed displacements
	// we save the prescribed displacements increments in the ui vector
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
	m_rigidSolver.PrepStep(tp, ui);

	// initialize contact
	if (m_fem.SurfacePairConstraints() > 0) UpdateContact();

	// initialize nonlinear constraints
	if (m_fem.NonlinearConstraints() > 0) UpdateConstraints();

	// intialize material point data
	// NOTE: do this before the stresses are updated
	// TODO: does it matter if the stresses are updated before
	//       the material point data is initialized
	for (int i=0; i<mesh.Domains(); ++i) 
	{
		FEDomain& dom = mesh.Domain(i);
		if (dom.IsActive()) dom.PreSolveUpdate(tp);
	}

	// update stresses
	UpdateModel();

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
// Performs the quasi-newton iterations.
bool FESolidSolver2::Quasin()
{
	vector<double> u0(m_neq);
	vector<double> Rold(m_neq);

	// convergence norms
	double	normR1;		// residual norm
	double	normE1;		// energy norm
	double	normU;		// displacement norm
	double	normu;		// displacement increment norm
	double	normRi;		// initial residual norm
	double	normEi;		// initial energy norm
	double	normEm;		// max energy norm
	double	normUi;		// initial displacement norm

	// Get the current step
	FEAnalysis* pstep = m_fem.GetCurrentStep();

	// set the time information
	const FETimeInfo& tp = m_fem.GetTime();

	// prepare for the first iteration
	PrepStep();

	// Initialize the QN-method
	if (QNInit() == false) return false;

	// loop until converged or when max nr of reformations reached
	bool bconv = false;		// convergence flag
	do
	{
		Logfile::MODE oldmode = felog.GetMode();
		if ((pstep->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
			(pstep->GetPrintLevel() != FE_PRINT_NEVER)) felog.SetMode(Logfile::LOG_FILE);

		felog.printf(" %d\n", m_niter+1);
		felog.SetMode(oldmode);

		// assume we'll converge. 
		bconv = true;

		// solve the equations (returns line search; solution stored in m_ui)
		double s = QNSolve();

		// set initial convergence norms
		if (m_niter == 0)
		{
			normRi = fabs(m_R0*m_R0);
			normEi = fabs(m_ui*m_R0);
			normUi = fabs(m_ui*m_ui);
			normEm = normEi;
		}

		// calculate actual displacement increment
		// NOTE: We don't apply the line search directly to m_ui since we need the unscaled search direction for the QN update below
		int neq = (int)m_Ui.size();
		vector<double> ui(m_ui);
		for (int i = 0; i<neq; ++i) ui[i] *= s;

		// update total displacements
		UpdateIncrements(m_Ui, ui, false);

		// calculate norms
		normR1 = m_R1*m_R1;
		normu  = ui*ui;
		normU  = m_Ui*m_Ui;
		normE1 = fabs(ui*m_R1);

		// check for nans
		if (ISNAN(normR1) || ISNAN(normu)) throw NANDetected();

		// check residual norm
		if ((m_Rtol > 0) && (normR1 > m_Rtol*normRi)) bconv = false;	

		// check displacement norm
		if ((m_Dtol > 0) && (normu  > (m_Dtol*m_Dtol)*normU )) bconv = false;

		// check energy norm
		if ((m_Etol > 0) && (normE1 > m_Etol*normEi)) bconv = false;

		// check linestep size
		if ((m_lineSearch->m_LStol > 0) && (s < m_lineSearch->m_LSmin)) bconv = false;

		// check energy divergence
		if (m_bdivreform)
		{
			if (normE1 > normEm) bconv = false;
		}

		// print convergence summary
		oldmode = felog.GetMode();
		if ((pstep->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
			(pstep->GetPrintLevel() != FE_PRINT_NEVER)) felog.SetMode(Logfile::LOG_FILE);

		felog.printf(" Nonlinear solution status: time= %lg\n", tp.currentTime);
		felog.printf("\tstiffness updates             = %d\n", m_strategy->m_nups);
		felog.printf("\tright hand side evaluations   = %d\n", m_nrhs);
		felog.printf("\tstiffness matrix reformations = %d\n", m_nref);
		if (m_lineSearch->m_LStol > 0) felog.printf("\tstep from line search         = %lf\n", s);
		felog.printf("\tconvergence norms :     INITIAL         CURRENT         REQUIRED\n");
		felog.printf("\t   residual         %15le %15le %15le \n", normRi, normR1, m_Rtol*normRi);
		felog.printf("\t   energy           %15le %15le %15le \n", normEi, normE1, m_Etol*normEi);
		felog.printf("\t   displacement     %15le %15le %15le \n", normUi, normu ,(m_Dtol*m_Dtol)*normU );

		felog.SetMode(oldmode);

		// see if we may have a small residual
		if ((bconv == false) && (normR1 < m_Rmin))
		{
			// check for almost zero-residual on the first iteration
			// this might be an indication that there is no force on the system
			felog.printbox("WARNING", "No force acting on the system.");
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
			// do additional checks that may trigger a stiffness reformation
			if (s < m_lineSearch->m_LSmin)
			{
				// check for zero linestep size
				felog.printbox("WARNING", "Zero linestep size. Stiffness matrix will now be reformed");
				QNForceReform(true);
			}
			else if ((normE1 > normEm) && m_bdivreform)
			{
				// check for diverging
				felog.printbox("WARNING", "Problem is diverging. Stiffness matrix will now be reformed");
				normEm = normE1;
				normEi = normE1;
				normRi = normR1;
				QNForceReform(true);
			}

			// Do the QN update (This may also do a stiffness reformation if necessary)
			bool bret = QNUpdate();

			// something went wrong with the update, so we'll need to break
			if (bret == false) break;
		}
		else if (m_baugment)
		{
			// do the augmentations
			bconv = DoAugmentations();
		}
	
		// increase iteration number
		m_niter++;

		// let's flush the logfile to make sure the last output will not get lost
		felog.flush();

		// do minor iterations callbacks
		m_fem.DoCallback(CB_MINOR_ITERS);
	}
	while (bconv == false);

	// if converged we update the total displacements
	if (bconv)
	{
        UpdateIncrementsEAS(m_Ui, false);
        UpdateIncrements(m_Ut, m_Ui, true);
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! Calculates global stiffness matrix.

bool FESolidSolver2::StiffnessMatrix()
{
	const FETimeInfo& tp = GetFEModel().GetTime();

	// get the mesh
	FEMesh& mesh = m_fem.GetMesh();

	// calculate the stiffness matrix for each domain
	for (int i=0; i<mesh.Domains(); ++i) 
	{
		if (mesh.Domain(i).IsActive()) 
		{
			FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
			dom.StiffnessMatrix(this);
		}
	}

	// calculate the body force stiffness matrix for each non-rigid domain
	for (int j = 0; j<m_fem.BodyLoads(); ++j)
	{
		FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(j));
		if (pbf && pbf->IsActive())
		{
			for (int i = 0; i<pbf->Domains(); ++i)
			{
				FEDomain& dom = *pbf->Domain(i);
				if (dom.IsActive() && dom.GetMaterial()->IsRigid() == false)
				{
					FEElasticDomain& edom = dynamic_cast<FEElasticDomain&>(dom);
					if (pbf) edom.BodyForceStiffness(this, *pbf);
				}
			}
        }
	}
    
    // TODO: add body force stiffness for rigid bodies

	// Add mass matrix for dynamic problems
	FEAnalysis* pstep = m_fem.GetCurrentStep();
	if (pstep->m_nanalysis == FE_DYNAMIC)
	{
		// scale factor
		double dt = tp.timeIncrement;
		double a = tp.alpham / (m_beta*dt*dt);

		// loop over all domains (except rigid)
		for (int i = 0; i<mesh.Domains(); ++i)
		{
			FEDomain& dom = mesh.Domain(i);
			if (dom.IsActive() && dom.GetMaterial()->IsRigid() == false)
			{
				FEElasticDomain& edom = dynamic_cast<FEElasticDomain&>(dom);
				edom.MassMatrix(this, a);
			}
		}

		m_rigidSolver.RigidMassMatrix(this, tp);
	}

	// calculate contact stiffness
	ContactStiffness();

	// calculate stiffness matrices for surface loads
	int nsl = m_fem.SurfaceLoads();
	for (int i = 0; i<nsl; ++i)
	{
		FESurfaceLoad* psl = m_fem.SurfaceLoad(i);
		if (psl->IsActive())
		{
			psl->StiffnessMatrix(tp, this); 
		}
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
void FESolidSolver2::NonLinearConstraintStiffness(const FETimeInfo& tp)
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

void FESolidSolver2::ContactStiffness()
{
	const FETimeInfo& tp = m_fem.GetTime();
	for (int i = 0; i<m_fem.SurfacePairConstraints(); ++i)
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(m_fem.SurfacePairConstraint(i));
		if (pci->IsActive()) pci->StiffnessMatrix(this, tp);
	}
}

//-----------------------------------------------------------------------------
// \todo I'd like to do something different with this. Right now, if a nodal load
//       it applied to a rigid body, the load has to be translated to a force and 
//       torque applied to the rigid body. Perhaps we should really define two types
//       of nodal loads, one for the deformable body and for the rigid body. This can
//       be done in a pre-processor phase. That way, standard assembly routines can be
//       used to assemble to loads into the global vector.
void FESolidSolver2::AssembleResidual(int node_id, int dof, double f, vector<double>& R)
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
void FESolidSolver2::AssembleStiffness(std::vector<int>& lm, matrix& ke)
{
	m_pK->Assemble(ke, lm);
}

//-----------------------------------------------------------------------------
// \todo adjust for rigid bodies
void FESolidSolver2::AssembleStiffness2(vector<int>& lmi, vector<int>& lmj, matrix& ke)
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
//!  Assembles the element stiffness matrix into the global stiffness matrix.
//!  Also adjusts the global stiffness matrix and residual to take the 
//!  prescribed displacements into account.

//! \todo In stead of changing the global stiffness matrix to accomodate for 
//!       the rigid bodies and linear constraints, can I modify the element stiffness
//!       matrix prior to assembly? I might have to change the elm vector as well as 
//!       the element matrix size.

void FESolidSolver2::AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke)
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
//! Calculates the contact forces
void FESolidSolver2::ContactForces(FEGlobalVector& R)
{
	const FETimeInfo& tp = m_fem.GetTime();
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

bool FESolidSolver2::Residual(vector<double>& R)
{
	TRACK_TIME("residual");

	// get the time information
	const FETimeInfo& tp = m_fem.GetTime();

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

	// calculate the internal (stress) forces
	for (int i=0; i<mesh.Domains(); ++i)
	{
        FEDomain& dom = mesh.Domain(i);
        if (dom.IsActive() && dom.GetMaterial()->IsRigid() == false)
        {
            FEElasticDomain& edom = dynamic_cast<FEElasticDomain&>(dom);
            edom.InternalForces(RHS);
        }
	}

	// extract the internal forces
	// (only when we really need it, below)
	bool logSolve = m_logSolve;
	vector<double> Rint;
	if (m_logSolve && m_fem.GetCurrentStep()->m_ntimesteps > 0)
	{
		Rint.resize(R.size());

		// we need to subtract the point forces since they are part of the external forces.
		for (int i = 0; i<Rint.size(); ++i)
			Rint[i] = RHS[i] - m_Fn[i];
	}

	// calculate the body forces
	for (int j = 0; j<m_fem.BodyLoads(); ++j)
	{
		FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(j));
		if (pbf && pbf->IsActive())
		{
			for (int i = 0; i<pbf->Domains(); ++i)
			{
				FEDomain& dom = *pbf->Domain(i);
				if (dom.IsActive() && dom.GetMaterial()->IsRigid() == false)
				{
					FEElasticDomain& edom = dynamic_cast<FEElasticDomain&>(dom);
					edom.BodyForce(RHS, *pbf);
				}
			}
        }
	}
    
    // calculate body forces for rigid bodies
    for (int j=0; j<m_fem.BodyLoads(); ++j)
    {
        FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(j));
		if (pbf && pbf->IsActive())
			m_rigidSolver.BodyForces(RHS, tp, *pbf);
    } 
    
    // calculate inertial forces for dynamic problems
    if (m_fem.GetCurrentStep()->m_nanalysis == FE_DYNAMIC)
    {
        // allocate F
        vector<double> F;
        
        // calculate the inertial forces for all elastic domains (except rigid domains)
        for (int nd = 0; nd < mesh.Domains(); ++nd)
        {
            FEDomain& dom = mesh.Domain(nd);
            if (dom.IsActive() && dom.GetMaterial()->IsRigid() == false)
            {
                FEElasticDomain& edom = dynamic_cast<FEElasticDomain&>(dom);
                edom.InertialForces(RHS, F);
            }
        }
        
        // update rigid bodies
        m_rigidSolver.InertialForces(RHS, tp);
    }
    
    // calculate forces due to surface loads
	int nsl = m_fem.SurfaceLoads();
	for (int i = 0; i<nsl; ++i)
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

	// forces due to point constraints
//	for (i=0; i<(int) fem.m_PC.size(); ++i) fem.m_PC[i]->Residual(this, R);

	// add model loads
	int NML = m_fem.ModelLoads();
	for (int i = 0; i<NML; ++i)
	{
		FEModelLoad& mli = *m_fem.ModelLoad(i);
		if (mli.IsActive())
		{
			mli.Residual(RHS, tp);
		}
	}

	// set the nodal reaction forces
	// TODO: Is this a good place to do this?
	for (int i = 0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		node.m_Fr = vec3d(0,0,0);

		int n;
		if ((n = -node.m_ID[m_dofX]-2) >= 0) node.m_Fr.x = -m_Fr[n];
		if ((n = -node.m_ID[m_dofY]-2) >= 0) node.m_Fr.y = -m_Fr[n];
		if ((n = -node.m_ID[m_dofZ]-2) >= 0) node.m_Fr.z = -m_Fr[n];
	}

	// apply the residual transformation
	// NOTE: This is an implementation of Ankush Aggarwal method to accelerate the Newton convergence
	if (m_logSolve && m_fem.GetCurrentStep()->m_ntimesteps > 0)
	{
		double TOL = 1.e-8;
		bool logused = false;
		vector<double> RHSlog;
		RHSlog.resize(R.size());
		for (int i = 0; i<Rint.size(); ++i)
		{
			if (fabs(RHS[i] - Rint[i])>TOL && fabs(Rint[i])>TOL && (Rint[i] - RHS[i]) / Rint[i]>0)
			{
				RHSlog[i] = -Rint[i] * log((Rint[i] - RHS[i]) / Rint[i]);
				logused = true;
			}
			else
			{
				RHSlog[i] = RHS[i];
			}
		}
		for (int i = 0; i<Rint.size(); ++i) R[i] = RHSlog[i];
		if (logused)
			felog.printf("Log method used\n");
	}

	// increase RHS counter
	m_nrhs++;

	return true;
}

//-----------------------------------------------------------------------------
//! calculate the nonlinear constraint forces 
void FESolidSolver2::NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp)
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
void FESolidSolver2::NodalForces(vector<double>& F, const FETimeInfo& tp)
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

