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
#include "FECore/BFGSSolver2.h"
#include "FECore/sys.h"
#include <FEBioMech/FEBodyForce.h>
#include <FECore/BC.h>
#include <FECore/FESurfaceLoad.h>
#include "FEFluidResistanceBC.h"
#include "FEBackFlowStabilization.h"
#include "FEFluidNormalVelocity.h"
#include <FECore/FEModelLoad.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FELinearConstraintManager.h>

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_PARAMETER_LIST(FEFluidSolver, FENewtonSolver)
	ADD_PARAMETER(m_Vtol         , FE_PARAM_DOUBLE, "vtol"        );
    ADD_PARAMETER(m_Dtol         , FE_PARAM_DOUBLE, "dtol"        );
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
//! FEFluidSolver Construction
//
FEFluidSolver::FEFluidSolver(FEModel* pfem) : FENewtonSolver(pfem)
{
    // default values
    m_Rtol = 0.001;
    m_Etol = 0.01;
    m_Vtol = 0.001;
    m_Dtol = 0.001;
    m_Rmin = 1.0e-20;
    
    m_nveq = 0;
    m_ndeq = 0;
    m_niter = 0;
    
    m_bsymm = false;
    m_bdivreform = true;
    m_bdoreforms = true;
	m_breformtimestep = true;

    m_rhoi = -1;
    m_pred = 0;
    
    m_baugment = false;
    
	// a different solution strategy is used here
	m_nqnsolver = QN_BFGS2;

	// turn off checking for a zero diagonal
	CheckZeroDiagonal(false);

	// Allocate degrees of freedom
	DOFS& dofs = pfem->GetDOFS();
	int nV = dofs.AddVariable("fluid velocity", VAR_VEC3);
	dofs.SetDOFName(nV, 0, "vx");
	dofs.SetDOFName(nV, 1, "vy");
	dofs.SetDOFName(nV, 2, "vz");
    int nE = dofs.AddVariable("fluid dilation", VAR_SCALAR);
	dofs.SetDOFName(nE, 0, "e");
    int nEP = dofs.AddVariable("previous fluid dilation", VAR_SCALAR);
    dofs.SetDOFName(nEP, 0, "ep");
    int nAE = dofs.AddVariable("fluid dilation tderiv", VAR_SCALAR);
    dofs.SetDOFName(nAE, 0, "ae");
    int nAEP = dofs.AddVariable("previous fluid dilation tderiv", VAR_SCALAR);
    dofs.SetDOFName(nAEP, 0, "aep");

	// get the dof indices
	m_dofVX = pfem->GetDOFIndex("vx");
	m_dofVY = pfem->GetDOFIndex("vy");
	m_dofVZ = pfem->GetDOFIndex("vz");
	m_dofE  = pfem->GetDOFIndex("e");
    
    m_dofEP  = pfem->GetDOFIndex("ep");
    m_dofAE  = pfem->GetDOFIndex("ae");
    m_dofAEP = pfem->GetDOFIndex("aep");
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

    // check parameters
    if (m_Vtol <  0.0) { felog.printf("Error: vtol must be nonnegative.\n"); return false; }
    if (m_Dtol <  0.0) { felog.printf("Error: dtol must be nonnegative.\n"); return false; }
    if (m_Etol <  0.0) { felog.printf("Error: etol must be nonnegative.\n"); return false; }
    if (m_Rtol <  0.0) { felog.printf("Error: rtol must be nonnegative.\n"); return false; }
    if (m_Rmin <  0.0) { felog.printf("Error: min_residual must be nonnegative.\n"  ); return false; }
    
    if (m_rhoi == -1) {
        m_alphaf = m_alpham = m_gamma = 1.0;
    }
    else if ((m_rhoi >= 0) && (m_rhoi <= 1)) {
        m_alphaf = 1.0/(1+m_rhoi);
        m_alpham = (3-m_rhoi)/(1+m_rhoi)/2;
        m_gamma = 0.5 + m_alpham - m_alphaf;
    }
    else { felog.printf("Error: rhoi must be -1 or between 0 and 1.\n"); return false; }
    
    // allocate vectors
    int neq = m_neq;
    m_Fn.assign(neq, 0);
    m_Fd.assign(neq, 0);
    m_Fr.assign(neq, 0);
    m_Ui.assign(neq, 0);
    m_Ut.assign(neq, 0);
    m_vi.assign(m_nveq,0);
    m_Vi.assign(m_nveq,0);
    m_di.assign(m_ndeq,0);
    m_Di.assign(m_ndeq,0);
    
    int i, n;
    
    // we need to fill the total DOF vector m_Ut
    // TODO: I need to find an easier way to do this
    FEMesh& mesh = m_fem.GetMesh();
    gather(m_Ut, mesh, m_dofVX);
    gather(m_Ut, mesh, m_dofVY);
    gather(m_Ut, mesh, m_dofVZ);
    gather(m_Ut, mesh, m_dofE);
    
    // TODO: move this somewhere else
    int nsl = m_fem.SurfaceLoads();
    for (int i=0; i<nsl; ++i)
    {
        FESurfaceLoad* psl = m_fem.SurfaceLoad(i);
        FEFluidResistanceBC* pfr = dynamic_cast<FEFluidResistanceBC*>(psl);
        FEFluidNormalVelocity* pnv = dynamic_cast<FEFluidNormalVelocity*>(psl);
        if (pfr && psl->IsActive()) pfr->MarkDilatation();
        else if (pnv && psl->IsActive() && pnv->m_bpv) pnv->MarkVelocity();
    }
    
    return true;
}

//-----------------------------------------------------------------------------
//! Initialize equations
bool FEFluidSolver::InitEquations()
{
    // base class initialization
    FENewtonSolver::InitEquations();
    
    int i;
    
    // determined the nr of velocity and dilatation equations
    FEMesh& mesh = m_fem.GetMesh();
    m_nveq = m_ndeq = 0;
    
    for (i=0; i<mesh.Nodes(); ++i)
    {
        FENode& n = mesh.Node(i);
        if (n.m_ID[m_dofVX] != -1) m_nveq++;
        if (n.m_ID[m_dofVY] != -1) m_nveq++;
        if (n.m_ID[m_dofVZ] != -1) m_nveq++;
        if (n.m_ID[m_dofE ] != -1) m_ndeq++;
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
        nid = n.m_ID[m_dofVX];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            vi[m++] = ui[nid];
            assert(m <= (int) vi.size());
        }
        nid = n.m_ID[m_dofVY];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            vi[m++] = ui[nid];
            assert(m <= (int) vi.size());
        }
        nid = n.m_ID[m_dofVZ];
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
        nid = n.m_ID[m_dofE];
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
        ar << m_Vtol << m_Dtol << m_Etol << m_Rtol << m_Rmin;
        ar << m_bsymm;
        ar << m_nrhs;
        ar << m_niter;
        ar << m_nref << m_ntotref;
        ar << m_naug;
    }
    else
    {
        ar >> m_Vtol >> m_Dtol >> m_Etol >> m_Rtol >> m_Rmin;
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
    int i, n;
    
    // get the mesh
    FEMesh& mesh = m_fem.GetMesh();
    
    // update nodes
    vector<double> U(m_Ut.size());
    for (size_t i=0; i<m_Ut.size(); ++i) U[i] = ui[i] + m_Ui[i] + m_Ut[i];
    
    scatter(U, mesh, m_dofVX);
    scatter(U, mesh, m_dofVY);
    scatter(U, mesh, m_dofVZ);
    scatter(U, mesh, m_dofE);
    
    // make sure the prescribed velocities are fullfilled
    int nvel = m_fem.PrescribedBCs();
    for (i=0; i<nvel; ++i)
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
    
    // update time derivatives of velocity and dilatation
    // for dynamic simulations
    FEAnalysis* pstep = m_fem.GetCurrentStep();
    if (pstep->m_nanalysis == FE_DYNAMIC)
    {
        int N = mesh.Nodes();
        double dt = pstep->m_dt;
        double cgi = 1 - 1.0/m_gamma;
        for (int i=0; i<N; ++i)
        {
            FENode& n = mesh.Node(i);
            
            // velocity time derivative
            vec3d vt = n.get_vec3d(m_dofVX, m_dofVY, m_dofVZ);
            n.m_at = n.m_ap*cgi + (vt - n.m_vp)/(m_gamma*dt);
            
            // dilatation time derivative
            double et = n.get(m_dofE);
            double ep = n.get(m_dofEP);
            double aep = n.get(m_dofAEP);
            double aet = aep*cgi + (et - ep)/(m_gamma*dt);
            n.set(m_dofAE, aet);
        }
    }
}

//-----------------------------------------------------------------------------
//! Updates the current state of the model
void FEFluidSolver::Update(vector<double>& ui)
{
	TimerTracker t(m_UpdateTime);

    // update kinematics
    UpdateKinematics(ui);
    
    // TODO: move this somewhere else
    int nsl = m_fem.SurfaceLoads();
    for (int i=0; i<nsl; ++i)
    {
        FESurfaceLoad* psl = m_fem.SurfaceLoad(i);
        FEFluidResistanceBC* pfr = dynamic_cast<FEFluidResistanceBC*>(psl);
        FEFluidNormalVelocity* pnv = dynamic_cast<FEFluidNormalVelocity*>(psl);
        if (pfr && psl->IsActive()) pfr->SetDilatation();
        else if (pnv && psl->IsActive() && pnv->m_bpv) pnv->SetVelocity();
    }
    
    // update element stresses
    UpdateStresses();
}

//-----------------------------------------------------------------------------
//!  Updates the element stresses
void FEFluidSolver::UpdateStresses()
{
    FEMesh& mesh = m_fem.GetMesh();
	FETimeInfo tp = m_fem.GetTime();
    tp.alpha = m_alphaf;
    tp.beta  = m_alpham;
    tp.gamma = m_gamma;
    
    // update the stresses on all domains
    for (int i=0; i<mesh.Domains(); ++i) mesh.Domain(i).Update(tp);
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
    FETimeInfo tp = m_fem.GetTime();
    tp.alpha = m_alphaf;
    tp.beta  = m_alpham;
    tp.gamma = m_gamma;
    
    // Assume we will pass (can't hurt to be optimistic)
    bool bconv = true;
    
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
    
    // get the time information
    FETimeInfo tp = fem.GetTime();
    tp.alpha = m_alphaf;
    tp.beta  = m_alpham;
    tp.gamma = m_gamma;
    
    // evaluate load curve values at current (or intermediate) time
    double t = tp.currentTime;
    double dt = tp.timeIncrement;
    double ta = (t > 0) ? t - (1-m_alphaf)*dt : m_alphaf*dt;
    
    return FESolver::InitStep(ta);
}

//-----------------------------------------------------------------------------
//! Prepares the data for the first BFGS-iteration.
void FEFluidSolver::PrepStep(const FETimeInfo& timeInfo)
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
    
    // store previous mesh state
    // we need them for strain and acceleration calculations
    FEMesh& mesh = m_fem.GetMesh();
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& ni = mesh.Node(i);
        ni.m_rp = ni.m_rt = ni.m_r0;
        ni.m_vp = ni.get_vec3d(m_dofVX, m_dofVY, m_dofVZ);
        ni.m_ap = ni.m_at;
        ni.set(m_dofEP, ni.get(m_dofE));
        ni.set(m_dofAEP, ni.get(m_dofAE));
        
        switch (m_pred) {
            case 0:
                // initial guess at start of new time step (default)
                ni.m_at = ni.m_ap*(m_gamma-1)/m_gamma;
                ni.set(m_dofAE, ni.get(m_dofAEP)*(m_gamma-1)/m_gamma);
                break;
                
            case 1:
                // initial guess at start of new time step (Zero Ydot)
                ni.m_at = vec3d(0,0,0);
                ni.set(m_dofAE, 0);
                
                ni.set_vec3d(m_dofVX, m_dofVY, m_dofVZ, ni.m_vp + ni.m_ap*dt*(1-m_gamma)*m_alphaf);
                ni.set(m_dofE, ni.get(m_dofEP) + ni.get(m_dofAEP)*dt*(1-m_gamma)*m_alphaf);
                break;
                
            case 2:
                // initial guess at start of new time step (Same Ydot)
                ni.m_at = ni.m_ap;
                ni.set(m_dofAE, ni.get(m_dofAEP));
                
                ni.set_vec3d(m_dofVX, m_dofVY, m_dofVZ, ni.m_vp + ni.m_ap*dt);
                ni.set(m_dofE, ni.get(m_dofEP) + ni.get(m_dofAEP)*dt);
                break;
                
            default:
                break;
        }
    }
    
    FETimeInfo tp = m_fem.GetTime();
    tp.alpha = m_alphaf;
    tp.beta  = m_alpham;
    tp.gamma = m_gamma;
    
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
    
    // intialize material point data
    // NOTE: do this before the stresses are updated
    // TODO: does it matter if the stresses are updated before
    //       the material point data is initialized
	// update domain data
    for (int i=0; i<mesh.Domains(); ++i) mesh.Domain(i).PreSolveUpdate(tp);
    
    // update stresses
    UpdateStresses();
    
    // see if we have to do nonlinear constraint augmentations
    if (m_fem.NonlinearConstraints() != 0) m_baugment = true;
}

//-----------------------------------------------------------------------------
//! Implements the BFGS2 algorithm to solve the nonlinear FE equations.
bool FEFluidSolver::Quasin(double time)
{
    int i;
    
    vector<double> u0(m_neq);
    vector<double> Rold(m_neq);
    
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
    
    // initialize flags
    bool bconv = false;		// convergence flag
    
    // Get the current step
    FEAnalysis* pstep = m_fem.GetCurrentStep();
    
    // prepare for the first iteration
	FETimeInfo tp = m_fem.GetTime();
    tp.alpha = m_alphaf;
    tp.beta  = m_alpham;
    tp.gamma = m_gamma;
    PrepStep(tp);
    
	// calculate initial stiffness matrix
	bool breform = m_breformtimestep;
	if (pstep->m_ntotiter == 0) breform = true;
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
            
            // update velocities
            for (i=0; i<m_nveq; ++i) m_Vi[i] += s*m_vi[i];

            // update dilatations
            for (i=0; i<m_ndeq; ++i) m_Di[i] += s*m_di[i];
            
            // calculate the norms
            normR1 = m_R1*m_R1;
			normv  = (m_vi*m_vi)*(s*s);
            normV  = m_Vi*m_Vi;
            normd  = (m_di*m_di)*(s*s);
            normD  = m_Di*m_Di;
            normE1 = s*fabs(m_ui*m_R1);
        
			// check for nans
			if (ISNAN(normR1)) throw NANDetected();
		}
		m_UpdateTime.stop();
        
        // check residual norm
        if ((m_Rtol > 0) && (normR1 > m_Rtol*normRi)) bconv = false;
        
        // check velocity norm
        if ((m_Vtol > 0) && (normv  > (m_Vtol*m_Vtol)*normV )) bconv = false;
        
        // check dilatation norm
        if ((m_Dtol > 0) && (normd  > (m_Dtol*m_Dtol)*normD )) bconv = false;
        
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
        felog.printf("\t   velocity         %15le %15le %15le \n", normVi, normv ,(m_Vtol*m_Vtol)*normV );
        felog.printf("\t   dilatation       %15le %15le %15le \n", normDi, normd ,(m_Dtol*m_Dtol)*normD );
        
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
                normVi = normv;
                normDi = normd;
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
        m_Ut += m_Ui;
    }
    
    return bconv;
}

//-----------------------------------------------------------------------------
//! Calculates global stiffness matrix.

bool FEFluidSolver::StiffnessMatrix(const FETimeInfo& tp)
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
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        dom.StiffnessMatrix(this, tp);
    }
    
    // calculate the body force stiffness matrix for each domain
    for (i=0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        int NBL = m_fem.BodyLoads();
        for (int j=0; j<NBL; ++j)
        {
            FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(j));
            if (pbf) dom.BodyForceStiffness(this, tp, *pbf);
        }
    }
    
    // calculate stiffness matrix due to surface loads
    int nsl = m_fem.SurfaceLoads();
    for (int i=0; i<nsl; ++i)
    {
        FESurfaceLoad* psl = m_fem.SurfaceLoad(i);
        if (psl->IsActive()) psl->StiffnessMatrix(tp, this);
    }
    
    // Add mass matrix
    // loop over all domains
    for (i=0; i<mesh.Domains(); ++i)
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
    // get nodal DOFS
    DOFS& fedofs = m_fem.GetDOFS();
    int MAX_NDOFS = fedofs.GetTotalDOFS();
    
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
//! calculates the residual vector
//! Note that the concentrated nodal forces are not calculated here.
//! This is because they do not depend on the geometry 
//! so we only calculate them once (in Quasin) and then add them here.

bool FEFluidSolver::Residual(vector<double>& R)
{
	TimerTracker t(m_RHSTime);

    // get the time information
    FETimeInfo tp = m_fem.GetTime();
    tp.alpha = m_alphaf;
    tp.beta  = m_alpham;
    tp.gamma = m_gamma;

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
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        for (int j=0; j<m_fem.BodyLoads(); ++j)
        {
            FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(j));
            dom.BodyForce(RHS, tp, *pbf);
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
        if ((n = -node.m_ID[m_dofVX]-2) >= 0) node.m_Fr.x = -m_Fr[n];
        if ((n = -node.m_ID[m_dofVY]-2) >= 0) node.m_Fr.y = -m_Fr[n];
        if ((n = -node.m_ID[m_dofVZ]-2) >= 0) node.m_Fr.z = -m_Fr[n];
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
