#include "stdafx.h"
#include "FEFluidSolver.h"
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
#include <FECore/BC.h>
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

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_PARAMETER_LIST(FEFluidSolver, FENewtonSolver)
	ADD_PARAMETER(m_Vtol         , FE_PARAM_DOUBLE, "vtol"        );
    ADD_PARAMETER(m_Ftol         , FE_PARAM_DOUBLE, "ftol"        );
    ADD_PARAMETER(m_Etol         , FE_PARAM_DOUBLE, "etol"        );
	ADD_PARAMETER(m_Rtol         , FE_PARAM_DOUBLE, "rtol"        );
	ADD_PARAMETER2(m_Rmin        , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "min_residual");
	ADD_PARAMETER2(m_Rmax        , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "max_residual");
	ADD_PARAMETER(m_bsymm        , FE_PARAM_BOOL  , "symmetric_stiffness");
    ADD_PARAMETER(m_rhoi         , FE_PARAM_DOUBLE, "rhoi"        );
    ADD_PARAMETER(m_pred         , FE_PARAM_INT   , "predictor"   );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! FEFluidSolver Construction
//
FEFluidSolver::FEFluidSolver(FEModel* pfem) : FENewtonSolver(pfem)
{
    // default values
    m_Rtol = 0.001;
    m_Etol = 0.01;
    m_Vtol = 0.001;
    m_Ftol = 0.001;
    m_Rmin = 1.0e-20;
	m_Rmax = 0;	// not used if zero
    
    m_nveq = 0;
    m_ndeq = 0;
    m_niter = 0;
    
    m_bsymm = false;

    m_rhoi = 0;
    m_pred = 0;
    
	// Preferred strategy is Broyden's method
	SetDefaultStrategy(QN_BROYDEN);

	// turn off checking for a zero diagonal
	CheckZeroDiagonal(false);

	// Allocate degrees of freedom
	DOFS& dofs = pfem->GetDOFS();
    int varD = dofs.AddVariable("displacement", VAR_VEC3);
    dofs.SetDOFName(varD, 0, "x");
    dofs.SetDOFName(varD, 1, "y");
    dofs.SetDOFName(varD, 2, "z");

    int nW = dofs.AddVariable("relative fluid velocity", VAR_VEC3);
    dofs.SetDOFName(nW, 0, "wx");
    dofs.SetDOFName(nW, 1, "wy");
    dofs.SetDOFName(nW, 2, "wz");
    int nE = dofs.AddVariable("fluid dilation", VAR_SCALAR);
	dofs.SetDOFName(nE, 0, "ef");
    
    int nWP = dofs.AddVariable("previous relative fluid velocity", VAR_VEC3);
    dofs.SetDOFName(nWP, 0, "wxp");
    dofs.SetDOFName(nWP, 1, "wyp");
    dofs.SetDOFName(nWP, 2, "wzp");
    int nEP = dofs.AddVariable("previous fluid dilation", VAR_SCALAR);
    dofs.SetDOFName(nEP, 0, "efp");
    
    int nAW = dofs.AddVariable("relative fluid acceleration", VAR_VEC3);
    dofs.SetDOFName(nAW, 0, "awx");
    dofs.SetDOFName(nAW, 1, "awy");
    dofs.SetDOFName(nAW, 2, "awz");
    int nAE = dofs.AddVariable("fluid dilation tderiv", VAR_SCALAR);
    dofs.SetDOFName(nAE, 0, "aef");
    
    int nAWP = dofs.AddVariable("previous relative fluid acceleration", VAR_VEC3);
    dofs.SetDOFName(nAWP, 0, "awxp");
    dofs.SetDOFName(nAWP, 1, "awyp");
    dofs.SetDOFName(nAWP, 2, "awzp");
    int nAEP = dofs.AddVariable("previous fluid dilation tderiv", VAR_SCALAR);
    dofs.SetDOFName(nAEP, 0, "aefp");

	// get the dof indices
    m_dofWX = pfem->GetDOFIndex("wx");
    m_dofWY = pfem->GetDOFIndex("wy");
    m_dofWZ = pfem->GetDOFIndex("wz");
	m_dofEF  = pfem->GetDOFIndex("ef");
    
    m_dofWXP = pfem->GetDOFIndex("wxp");
    m_dofWYP = pfem->GetDOFIndex("wyp");
    m_dofWZP = pfem->GetDOFIndex("wzp");
    m_dofEFP  = pfem->GetDOFIndex("efp");
    
    m_dofAWX = pfem->GetDOFIndex("awx");
    m_dofAWY = pfem->GetDOFIndex("awy");
    m_dofAWZ = pfem->GetDOFIndex("awz");
    m_dofAEF = pfem->GetDOFIndex("aef");

    m_dofAWXP = pfem->GetDOFIndex("awxp");
    m_dofAWYP = pfem->GetDOFIndex("awyp");
    m_dofAWZP = pfem->GetDOFIndex("awzp");
    m_dofAEFP = pfem->GetDOFIndex("aefp");
}

//-----------------------------------------------------------------------------
FEFluidSolver::~FEFluidSolver()
{

}

//-----------------------------------------------------------------------------
//! Allocates and initializes the data structures used by the FEFluidSolver
//
bool FEFluidSolver::Init()
{
	// initialize base class
	if (FENewtonSolver::Init() == false) return false;

	// set the block size of the sparse matrix
	m_plinsolve->SetPartition(m_nveq);

    // check parameters
    if (m_Vtol <  0.0) { felog.printf("Error: vtol must be nonnegative.\n"); return false; }
    if (m_Ftol <  0.0) { felog.printf("Error: dtol must be nonnegative.\n"); return false; }
    if (m_Etol <  0.0) { felog.printf("Error: etol must be nonnegative.\n"); return false; }
    if (m_Rtol <  0.0) { felog.printf("Error: rtol must be nonnegative.\n"); return false; }
    
    if (m_rhoi == -1) {
        m_alphaf = m_alpham = m_gammaf = 1.0;
    }
    else if ((m_rhoi >= 0) && (m_rhoi <= 1)) {
        m_alphaf = 1.0/(1+m_rhoi);
        m_alpham = (3-m_rhoi)/(1+m_rhoi)/2;
        m_gammaf = 0.5 + m_alpham - m_alphaf;
    }
    else { felog.printf("Error: rhoi must be -1 or between 0 and 1.\n"); return false; }
    
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
    
    // we need to fill the total DOF vector m_Ut
    // TODO: I need to find an easier way to do this
    FEMesh& mesh = m_fem.GetMesh();
    gather(m_Ut, mesh, m_dofWX);
    gather(m_Ut, mesh, m_dofWY);
    gather(m_Ut, mesh, m_dofWZ);
    gather(m_Ut, mesh, m_dofEF);
    
    return true;
}

//-----------------------------------------------------------------------------
//! Initialize equations
bool FEFluidSolver::InitEquations()
{
    // base class initialization
    FENewtonSolver::InitEquations();
    
    // determined the nr of velocity and dilatation equations
    FEMesh& mesh = m_fem.GetMesh();
    m_nveq = m_ndeq = 0;
    
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& n = mesh.Node(i);
        if (n.m_ID[m_dofWX] != -1) m_nveq++;
        if (n.m_ID[m_dofWY] != -1) m_nveq++;
        if (n.m_ID[m_dofWZ] != -1) m_nveq++;
        if (n.m_ID[m_dofEF ] != -1) m_ndeq++;
    }

    return true;
}

//-----------------------------------------------------------------------------
void FEFluidSolver::GetVelocityData(vector<double> &vi, vector<double> &ui)
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
void FEFluidSolver::GetDilatationData(vector<double> &ei, vector<double> &ui)
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

void FEFluidSolver::Serialize(DumpStream& ar)
{
    // Serialize parameters
    FENewtonSolver::Serialize(ar);
    
    if (ar.IsSaving())
    {
        ar << m_Vtol << m_Ftol << m_Etol << m_Rtol << m_Rmin << m_Rmax;
        ar << m_bsymm;
        ar << m_nrhs;
        ar << m_niter;
        ar << m_nref << m_ntotref;
        ar << m_naug;
    }
    else
    {
		ar >> m_Vtol >> m_Ftol >> m_Etol >> m_Rtol >> m_Rmin >> m_Rmax;
        ar >> m_bsymm;
        ar >> m_nrhs;
        ar >> m_niter;
        ar >> m_nref >> m_ntotref;
        ar >> m_naug;
    }
}

//-----------------------------------------------------------------------------
//! Update the kinematics of the model, such as nodal positions, velocities,
//! accelerations, etc.
void FEFluidSolver::UpdateKinematics(vector<double>& ui)
{
    // get the mesh
    FEMesh& mesh = m_fem.GetMesh();
    
    // update nodes
    vector<double> U(m_Ut.size());
    for (size_t i=0; i<m_Ut.size(); ++i) U[i] = ui[i] + m_Ui[i] + m_Ut[i];
    
    scatter(U, mesh, m_dofWX);
    scatter(U, mesh, m_dofWY);
    scatter(U, mesh, m_dofWZ);
    scatter(U, mesh, m_dofEF);
    
    // make sure the prescribed velocities are fullfilled
    int nvel = m_fem.PrescribedBCs();
    for (int i=0; i<nvel; ++i)
    {
        FEPrescribedBC& dc = *m_fem.PrescribedBC(i);
        if (dc.IsActive()) dc.Update();
    }

    // prescribe DOFs for specialized surface loads
    int nsl = m_fem.SurfaceLoads();
    for (int i=0; i<nsl; ++i)
    {
        FESurfaceLoad& psl = *m_fem.SurfaceLoad(i);
        if (psl.IsActive()) psl.Update();
    }
    
	// enforce the linear constraints
	// TODO: do we really have to do this? Shouldn't the algorithm
	// already guarantee that the linear constraints are satisfied?
	FELinearConstraintManager& LCM = m_fem.GetLinearConstraintManager();
	if (LCM.LinearConstraints() > 0)
	{
		LCM.Update();
	}
    
    // update time derivatives of velocity and dilatation
    // for dynamic simulations
    FEAnalysis* pstep = m_fem.GetCurrentStep();
    if (pstep->m_nanalysis == FE_DYNAMIC)
    {
        int N = mesh.Nodes();
		double dt = m_fem.GetTime().timeIncrement;
        double cgi = 1 - 1.0/m_gammaf;
        for (int i=0; i<N; ++i)
        {
            FENode& n = mesh.Node(i);
            
            // velocity time derivative
            vec3d vft = n.get_vec3d(m_dofWX, m_dofWY, m_dofWZ);
            vec3d vfp = n.get_vec3d(m_dofWXP, m_dofWYP, m_dofWZP);
            vec3d aft = n.get_vec3d(m_dofAWX, m_dofAWY, m_dofAWZ);
            vec3d afp = n.get_vec3d(m_dofAWXP, m_dofAWYP, m_dofAWZP);
            aft = afp*cgi + (vft - vfp)/(m_gammaf*dt);
            n.set_vec3d(m_dofAWX, m_dofAWY, m_dofAWZ, aft);
            
            // dilatation time derivative
            double eft = n.get(m_dofEF);
            double efp = n.get(m_dofEFP);
            double aefp = n.get(m_dofAEFP);
            double aeft = aefp*cgi + (eft - efp)/(m_gammaf*dt);
            n.set(m_dofAEF, aeft);
        }
    }
}

//-----------------------------------------------------------------------------
//! Updates the current state of the model
void FEFluidSolver::Update(vector<double>& ui)
{
	TRACK_TIME("update");

    // update kinematics
    UpdateKinematics(ui);
    
    // update contact
    if (m_fem.SurfacePairConstraints() > 0) UpdateContact();
    
    // update element stresses
	UpdateModel();
}

//-----------------------------------------------------------------------------
//!  Updates the element stresses
void FEFluidSolver::UpdateModel()
{
    FEMesh& mesh = m_fem.GetMesh();
	const FETimeInfo& tp = GetFEModel().GetTime();
	
    // update the stresses on all domains
    for (int i=0; i<mesh.Domains(); ++i) mesh.Domain(i).Update(tp);
}

//-----------------------------------------------------------------------------
//! Update contact interfaces.
void FEFluidSolver::UpdateContact()
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
bool FEFluidSolver::Augment()
{
	const FETimeInfo& tp = GetFEModel().GetTime();
    
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
    
    return bconv;
}

//-----------------------------------------------------------------------------
bool FEFluidSolver::InitStep(double time)
{
    FEModel& fem = GetFEModel();
    
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
void FEFluidSolver::PrepStep()
{
	TRACK_TIME("update");

	const FETimeInfo& tp = m_fem.GetTime();
	double dt = tp.timeIncrement;

    // zero total DOFs
    zero(m_Ui);
    zero(m_Vi);
    zero(m_Di);
    
    // store previous mesh state
    // we need them for strain and acceleration calculations
    FEMesh& mesh = m_fem.GetMesh();
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& ni = mesh.Node(i);
        ni.m_rp = ni.m_rt = ni.m_r0;
        ni.set_vec3d(m_dofWXP, m_dofWYP, m_dofWZP, ni.get_vec3d(m_dofWX, m_dofWY, m_dofWZ));
        ni.set_vec3d(m_dofAWXP, m_dofAWYP, m_dofAWZP, ni.get_vec3d(m_dofAWX, m_dofAWY, m_dofAWZ));
        ni.set(m_dofEFP, ni.get(m_dofEF));
        ni.set(m_dofAEFP, ni.get(m_dofAEF));
        
        switch (m_pred) {
            case 0:
            {
                // initial guess at start of new time step (default)
                vec3d afp = ni.get_vec3d(m_dofAWXP, m_dofAWYP, m_dofAWZP);
                ni.set_vec3d(m_dofAWX, m_dofAWY, m_dofAWZ, afp*(m_gammaf-1)/m_gammaf);
                ni.set(m_dofAEF, ni.get(m_dofAEFP)*(m_gammaf-1)/m_gammaf);
            }
                break;
                
            case 1:
            {
                // initial guess at start of new time step (Zero Ydot)
                ni.set_vec3d(m_dofAWX, m_dofAWY, m_dofAWZ,vec3d(0,0,0));
                ni.set(m_dofAEF, 0);
                
                vec3d vfp = ni.get_vec3d(m_dofWXP, m_dofWYP, m_dofWZP);
                vec3d afp = ni.get_vec3d(m_dofAWXP, m_dofAWYP, m_dofAWZP);
                ni.set_vec3d(m_dofWX, m_dofWY, m_dofWZ, vfp + afp*dt*(1-m_gammaf)*m_alphaf);
                ni.set(m_dofEF, ni.get(m_dofEFP) + ni.get(m_dofAEFP)*dt*(1-m_gammaf)*m_alphaf);
            }
                break;
                
            case 2:
            {
                // initial guess at start of new time step (Same Ydot)
                vec3d afp = ni.get_vec3d(m_dofAWXP, m_dofAWYP, m_dofAWZP);
                ni.set_vec3d(m_dofAWX, m_dofAWY, m_dofAWZ, afp);
                ni.set(m_dofAEF, ni.get(m_dofAEFP));
                
                vec3d vfp = ni.get_vec3d(m_dofWXP, m_dofWYP, m_dofWZP);
                ni.set_vec3d(m_dofWX, m_dofWY, m_dofWZ, vfp + afp*dt);
                ni.set(m_dofEF, ni.get(m_dofEFP) + ni.get(m_dofAEFP)*dt);
            }
                break;
                
            default:
                break;
        }
    }
    
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
    
    // apply prescribed DOFs for specialized surface loads
    int nsl = m_fem.SurfaceLoads();
    for (int i=0; i<nsl; ++i)
    {
        FESurfaceLoad& psl = *m_fem.SurfaceLoad(i);
        if (psl.IsActive()) psl.Update();
    }
    
    // initialize contact
    if (m_fem.SurfacePairConstraints() > 0) UpdateContact();
    
    // intialize material point data
    // NOTE: do this before the stresses are updated
    // TODO: does it matter if the stresses are updated before
    //       the material point data is initialized
	// update domain data
    for (int i=0; i<mesh.Domains(); ++i) mesh.Domain(i).PreSolveUpdate(tp);
    
    // update stresses
	UpdateModel();
    
    // see if we need to do contact augmentations
    m_baugment = false;
    for (int i = 0; i<m_fem.SurfacePairConstraints(); ++i)
    {
        FEContactInterface& ci = dynamic_cast<FEContactInterface&>(*m_fem.SurfacePairConstraint(i));
        if (ci.IsActive() && ci.m_blaugon) m_baugment = true;
    }
    
    // see if we have to do nonlinear constraint augmentations
    if (m_fem.NonlinearConstraints() != 0) m_baugment = true;
}

//-----------------------------------------------------------------------------
bool FEFluidSolver::Quasin()
{
    // convergence norms
    double	normR1;		// residual norm
    double	normE1;		// energy norm
    double	normV;		// velocity norm
    double	normv;		// velocity increment norm
    double	normRi = 0;	// initial residual norm
    double	normVi = 0;	// initial velocity norm
    double	normEi = 0; // initial energy norm
    double	normEm = 0;	// max energy norm
    double	normDi = 0;	// initial dilatation norm
    double	normD;		// current dilatation norm
    double	normd;		// incremement dilatation norm
    
    // Get the current step
    FEAnalysis* pstep = m_fem.GetCurrentStep();
    
    // prepare for the first iteration
	const FETimeInfo& tp = m_fem.GetTime();
    PrepStep();
    
	// Init QN method
	if (QNInit() == false) return false;
    
    // loop until converged or when max nr of reformations reached
	bool bconv = false; // convergence flag
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
        felog.printf("\t   velocity         %15le %15le %15le \n", normVi, normv ,(m_Vtol*m_Vtol)*normV );
        felog.printf("\t   dilatation       %15le %15le %15le \n", normDi, normd ,(m_Ftol*m_Ftol)*normD );
        
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
                normVi = normv;
                normDi = normd;
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
        
        // let's flush the logfile to make sure the last output will not get lost
        felog.flush();
        
        // do minor iterations callbacks
        m_fem.DoCallback(CB_MINOR_ITERS);
    }
    while (bconv == false);
    
    // if converged we update the total velocities
    if (bconv)
    {
        m_Ut += m_Ui;
    }
    
    return bconv;
}

//-----------------------------------------------------------------------------
//! Calculates global stiffness matrix.

bool FEFluidSolver::StiffnessMatrix()
{
	const FETimeInfo& tp = GetFEModel().GetTime();

    // get the mesh
    FEMesh& mesh = m_fem.GetMesh();
    
    // calculate the stiffness matrix for each domain
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        dom.StiffnessMatrix(this, tp);
    }
    
    // calculate the body force stiffness matrix for each domain
	int NBL = m_fem.BodyLoads();
	for (int j = 0; j<NBL; ++j)
	{
		FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(j));
		if (pbf && pbf->IsActive())
		{
			for (int i = 0; i<pbf->Domains(); ++i)
			{
				FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(*pbf->Domain(i));
				dom.BodyForceStiffness(this, tp, *pbf);
			}
        }
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
    
    // Add mass matrix
    // loop over all domains
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        dom.MassMatrix(this, tp);
    }
    
    // calculate nonlinear constraint stiffness
    // note that this is the contribution of the
    // constrainst enforced with augmented lagrangian
    NonLinearConstraintStiffness(tp);
    
    return true;
}

//-----------------------------------------------------------------------------
//! Calculate the stiffness contribution due to nonlinear constraints
void FEFluidSolver::NonLinearConstraintStiffness(const FETimeInfo& tp)
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

void FEFluidSolver::ContactStiffness()
{
    const FETimeInfo& tp = m_fem.GetTime();
    for (int i = 0; i<m_fem.SurfacePairConstraints(); ++i)
    {
        FEContactInterface* pci = dynamic_cast<FEContactInterface*>(m_fem.SurfacePairConstraint(i));
        if (pci->IsActive()) pci->StiffnessMatrix(this, tp);
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolver::AssembleResidual(int node_id, int dof, double f, vector<double>& R)
{
    // get the mesh
    FEMesh& mesh = m_fem.GetMesh();
    
    // get the equation number
    FENode& node = mesh.Node(node_id);
    int n = node.m_ID[dof];
    
    // assemble into global vector
    if (n >= 0)
#pragma omp atomic
        R[n] += f;
}

//-----------------------------------------------------------------------------
//! \todo This function is only used for rigid joints. I need to figure out if
//!       I can use the other assembly function.
void FEFluidSolver::AssembleStiffness(std::vector<int>& lm, matrix& ke)
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

void FEFluidSolver::AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke)
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
            if (J >= 0)
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
}

//-----------------------------------------------------------------------------
//! Calculates the contact forces
void FEFluidSolver::ContactForces(FEGlobalVector& R)
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

bool FEFluidSolver::Residual(vector<double>& R)
{
	TRACK_TIME("residual");

    // get the time information
	const FETimeInfo& tp = GetFEModel().GetTime();

    // initialize residual with concentrated nodal loads
    R = m_Fn;
    
    // zero nodal reaction forces
    zero(m_Fr);
    
    // setup the global vector
    FEFluidResidualVector RHS(GetFEModel(), R, m_Fr);
    
    // get the mesh
    FEMesh& mesh = m_fem.GetMesh();
    
    // set flag for transient or steady-state analyses
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        if (m_fem.GetCurrentStep()->m_nanalysis == FE_STEADY_STATE)
            dom.SetSteadyStateAnalysis();
        else
            dom.SetTransientAnalysis();
    }
    
    // calculate the internal (stress) forces
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        dom.InternalForces(RHS, tp);
    }
    
    // calculate the body forces
	for (int j = 0; j<m_fem.BodyLoads(); ++j)
	{
		FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(j));
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
        if ((n = -node.m_ID[m_dofWX]-2) >= 0) node.m_Fr.x = -m_Fr[n];
        if ((n = -node.m_ID[m_dofWY]-2) >= 0) node.m_Fr.y = -m_Fr[n];
        if ((n = -node.m_ID[m_dofWZ]-2) >= 0) node.m_Fr.z = -m_Fr[n];
    }
    
    // increase RHS counter
    m_nrhs++;
    
    return true;
}

//-----------------------------------------------------------------------------
//! calculate the nonlinear constraint forces
void FEFluidSolver::NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp)
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

void FEFluidSolver::NodalForces(vector<double>& F, const FETimeInfo& tp)
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
