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
#include <FECore/FEModelLoad.h>
#include "FECore/FEAnalysis.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_PARAMETER_LIST(FEFluidSolver, FENewtonSolver)
	ADD_PARAMETER(m_Vtol         , FE_PARAM_DOUBLE, "vtol"        );
    ADD_PARAMETER(m_Etol         , FE_PARAM_DOUBLE, "etol"        );
	ADD_PARAMETER(m_Rtol         , FE_PARAM_DOUBLE, "rtol"        );
	ADD_PARAMETER(m_Rmin         , FE_PARAM_DOUBLE, "min_residual");
	ADD_PARAMETER(m_bdivreform   , FE_PARAM_BOOL  , "diverge_reform");
	ADD_PARAMETER(m_bdoreforms   , FE_PARAM_BOOL  , "do_reforms"  );
	ADD_PARAMETER(m_bsymm        , FE_PARAM_BOOL  , "symmetric_stiffness");
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
    m_Rmin = 1.0e-20;
    
    m_niter = 0;
    
    m_bsymm = false;
    m_bdivreform = true;
    m_bdoreforms = true;

    m_baugment = false;
    
	// a different solution strategy is used here
	m_nqnsolver = QN_BFGS2;

	// turn off checking for a zero diagonal
	CheckZeroDiagonal(false);

	// Allocate degrees of freedom
	DOFS& dofs = pfem->GetDOFS();
	int nV = dofs.AddVariable("fluid velocity", VAR_VEC3);
	int nE = dofs.AddVariable("fluid dilation", VAR_SCALAR);
	dofs.SetDOFName(nV, 0, "vx");
	dofs.SetDOFName(nV, 1, "vy");
	dofs.SetDOFName(nV, 2, "vz");
	dofs.SetDOFName(nE, 0, "e");

	// get the dof indices
	m_dofVX = pfem->GetDOFIndex("vx");
	m_dofVY = pfem->GetDOFIndex("vy");
	m_dofVZ = pfem->GetDOFIndex("vz");
	m_dofE  = pfem->GetDOFIndex("e");
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
    if (m_Etol <  0.0) { felog.printf("Error: etol must be nonnegative.\n"); return false; }
    if (m_Rtol <  0.0) { felog.printf("Error: rtol must be nonnegative.\n"); return false; }
    if (m_Rmin <  0.0) { felog.printf("Error: min_residual must be nonnegative.\n"  ); return false; }
    
    // allocate vectors
    int neq = m_neq;
    m_Fn.assign(neq, 0);
    m_Fd.assign(neq, 0);
    m_Fr.assign(neq, 0);
    m_Vi.assign(neq, 0);
    m_Vt.assign(neq, 0);
    
    int i, n;
    
    // we need to fill the total velocity vector m_Vt
    // TODO: I need to find an easier way to do this
    FEMesh& mesh = m_fem.GetMesh();
    for (i=0; i<mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        
        // velocity dofs
        n = node.m_ID[m_dofVX]; if (n >= 0) m_Vt[n] = node.get(m_dofVX);
        n = node.m_ID[m_dofVY]; if (n >= 0) m_Vt[n] = node.get(m_dofVY);
        n = node.m_ID[m_dofVZ]; if (n >= 0) m_Vt[n] = node.get(m_dofVZ);
        
        // dilatation dofs
        n = node.m_ID[m_dofE]; if (n >= 0) m_Vt[n] = node.get(m_dofE);
    }
    
    return true;
}

//-----------------------------------------------------------------------------
//! Save data to dump file

void FEFluidSolver::Serialize(DumpStream& ar)
{
    // Serialize parameters
    FENewtonSolver::Serialize(ar);
    
    if (ar.IsSaving())
    {
        ar << m_Vtol << m_Etol << m_Rtol << m_Rmin;
        ar << m_bsymm;
        ar << m_nrhs;
        ar << m_niter;
        ar << m_nref << m_ntotref;
        ar << m_naug;
    }
    else
    {
        ar >> m_Vtol >> m_Etol >> m_Rtol >> m_Rmin;
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
void FEFluidSolver::UpdateKinematics(vector<double>& vi)
{
    int i, n;
    
    // get the mesh
    FEMesh& mesh = m_fem.GetMesh();
    
    // update nodes
    for (i=0; i<mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        
        // velocity dofs
        // current velocity = total at prev conv step + total increment so far + current increment
        if ((n = node.m_ID[m_dofVX]) >= 0) node.set(m_dofVX, m_Vt[n] + m_Vi[n] + vi[n]);
        if ((n = node.m_ID[m_dofVY]) >= 0) node.set(m_dofVY, m_Vt[n] + m_Vi[n] + vi[n]);
        if ((n = node.m_ID[m_dofVZ]) >= 0) node.set(m_dofVZ, m_Vt[n] + m_Vi[n] + vi[n]);
        
        // dilatation dofs
        // current dilatation = total at prev conv step + total increment so far + current increment
        if ((n = node.m_ID[m_dofE]) >= 0) node.set(m_dofE, m_Vt[n] + m_Vi[n] + vi[n]);
    }
    
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
    if (m_fem.m_LinC.size() > 0)
    {
        int nlin = (int)m_fem.m_LinC.size();
        list<FELinearConstraint>::iterator it = m_fem.m_LinC.begin();
        double d;
        for (int n=0; n<nlin; ++n, ++it)
        {
            FELinearConstraint& lc = *it;
            FENode& node = mesh.Node(lc.master.node);
            
            d = 0;
            int ns = (int)lc.slave.size();
            list<FELinearConstraint::SlaveDOF>::iterator si = lc.slave.begin();
            for (int i=0; i<ns; ++i, ++si)
            {
                FENode& node = mesh.Node(si->node);
                switch (si->bc)
                {
                    case 0: d += si->val*node.get(m_dofVX); break;
                    case 1: d += si->val*node.get(m_dofVY); break;
                    case 2: d += si->val*node.get(m_dofVZ); break;
                    case 3: d += si->val*node.get(m_dofE); break;
                }
            }
            
            switch (lc.master.bc)
            {
                case 0: node.set(m_dofVX, d); break;
                case 1: node.set(m_dofVY, d); break;
                case 2: node.set(m_dofVZ, d); break;
                case 3: node.set(m_dofE, d); break;
            }
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
    
    // update element stresses
    UpdateStresses();
}

//-----------------------------------------------------------------------------
//!  Updates the element stresses
void FEFluidSolver::UpdateStresses()
{
    FEMesh& mesh = m_fem.GetMesh();
	FETimeInfo tp = m_fem.GetTime();
    
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
//! Prepares the data for the first BFGS-iteration.
void FEFluidSolver::PrepStep(const FETimeInfo& timeInfo)
{
	TimerTracker t(m_UpdateTime);

    // initialize counters
    m_niter = 0;	// nr of iterations
    m_nrhs  = 0;	// nr of RHS evaluations
    m_nref  = 0;	// nr of stiffness reformations
    m_ntotref = 0;
    m_pbfgs->m_nups	= 0;	// nr of stiffness updates between reformations
    m_naug  = 0;	// nr of augmentations
    
    // zero total velocities
    zero(m_Vi);
    
    // store previous mesh state
    // we need them for strain and acceleration calculations
    FEMesh& mesh = m_fem.GetMesh();
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& ni = mesh.Node(i);
        ni.m_rp = ni.m_rt = ni.m_r0;
        ni.m_vp = ni.get_vec3d(m_dofVX, m_dofVY, m_dofVZ);
        ni.m_ap = ni.m_at;
    }
    
    // apply concentrated nodal forces
    // since these forces do not depend on the geometry
    // we can do this once outside the NR loop.
    NodalForces(m_Fn, timeInfo);
    
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
    FEMaterialPoint::dt = timeInfo.timeIncrement;
	FEMaterialPoint::time = timeInfo.currentTime;
    
	// update domain data
    for (int i=0; i<mesh.Domains(); ++i) mesh.Domain(i).PreSolveUpdate(timeInfo);
    
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
    double	normU;		// velocity norm
    double	normu;		// velocity increment norm
    double	normRi = 0;		// initial residual norm
    double	normUi = 0;		// initial velocity norm
    double	normEi;		// initial energy norm
    double	normEm;		// max energy norm
    
    // initialize flags
    bool bconv = false;		// convergence flag
    bool breform = false;	// reformation flag
    
    // Get the current step
    FEAnalysis* pstep = m_fem.GetCurrentStep();
    
    // prepare for the first iteration
	FETimeInfo tp = m_fem.GetTime();
    PrepStep(tp);
    
    // calculate initial stiffness matrix
    if (ReformStiffness(tp) == false) return false;
    
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
            (pstep->GetPrintLevel() != FE_PRINT_NEVER)) felog.SetMode(Logfile::FILE_ONLY);
        
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
        
        // set initial convergence norms
        if (m_niter == 0)
        {
            normRi = fabs(m_R0*m_R0);
            normEi = fabs(m_ui*m_R0);
            normUi = fabs(m_ui*m_ui);
            normEm = normEi;
        }
        
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
        normR1 = m_R1*m_R1;
        normu  = (m_ui*m_ui)*(s*s);
        normE1 = s*fabs(m_ui*m_R1);
        
        // check for nans
        if (ISNAN(normR1) || ISNAN(normu)) throw NANDetected();
        
        // update total velocities
        int neq = (int)m_Vi.size();
        for (i=0; i<neq; ++i) m_Vi[i] += s*m_ui[i];
        normU  = m_Vi*m_Vi;
        
        // check residual norm
        if ((m_Rtol > 0) && (normR1 > m_Rtol*normRi)) bconv = false;
        
        // check velocity norm
        if ((m_Vtol > 0) && (normu  > (m_Vtol*m_Vtol)*normU )) bconv = false;
        
        // check energy norm
        if ((m_Etol > 0) && (normE1 > m_Etol*normEi)) bconv = false;
        
        // check linestep size
        if ((m_LStol > 0) && (s < m_LSmin)) bconv = false;
        
        // check energy divergence
        if (normE1 > normEm) bconv = false;
        
        // print convergence summary
        oldmode = felog.GetMode();
        if ((pstep->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
            (pstep->GetPrintLevel() != FE_PRINT_NEVER)) felog.SetMode(Logfile::FILE_ONLY);
        
        felog.printf(" Nonlinear solution status: time= %lg\n", time);
        felog.printf("\tstiffness updates             = %d\n", m_pbfgs->m_nups);
        felog.printf("\tright hand side evaluations   = %d\n", m_nrhs);
        felog.printf("\tstiffness matrix reformations = %d\n", m_nref);
        if (m_LStol > 0) felog.printf("\tstep from line search         = %lf\n", s);
        felog.printf("\tconvergence norms :     INITIAL         CURRENT         REQUIRED\n");
        felog.printf("\t   residual         %15le %15le %15le \n", normRi, normR1, m_Rtol*normRi);
        felog.printf("\t   energy           %15le %15le %15le \n", normEi, normE1, m_Etol*normEi);
        felog.printf("\t   velocity         %15le %15le %15le \n", normUi, normu ,(m_Vtol*m_Vtol)*normU );
        
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
                        if (m_pbfgs->Update(s, m_ui, m_R0, m_R1) == false)
                        {
                            // Stiffness update has failed.
                            // this might be due a too large condition number
                            // or the update was no longer positive definite.
                            felog.printbox("WARNING", "The BFGS update has failed.\nStiffness matrix will now be reformed.");
                            breform = true;
                        }
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
            zero(m_ui);
            
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
            m_R0 = m_R1;
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
        if (mode != Logfile::NEVER)
        {
            felog.SetMode(Logfile::FILE_ONLY);
            felog.printf("\nconvergence summary\n");
            felog.printf("    number of iterations   : %d\n", m_niter);
            felog.printf("    number of reformations : %d\n", m_nref);
            felog.SetMode(mode);
        }
    }
    
    // if converged we update the total velocities
    if (bconv)
    {
        m_Vt += m_Vi;
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
        dom.StiffnessMatrix(this);
    }
    
    // calculate the body force stiffness matrix for each domain
    for (i=0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        int NBL = m_fem.BodyLoads();
        for (int j=0; j<NBL; ++j)
        {
            FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(j));
            if (pbf) dom.BodyForceStiffness(this, *pbf);
        }
    }
    
    // Add mass matrix
    // loop over all domains
    for (i=0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        dom.MassMatrix(this);
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
    if (n >= 0) R[n] += f;
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
    if (m_fem.m_LinC.size() > 0)
    {
        int i, j, l;
        
        int ndof = ke.rows();
        int ndn = ndof / en.size();
        
        SparseMatrix& K = *m_pK;
        
        // loop over all stiffness components 
        // and correct for linear constraints
        int ni, nj, li, lj, I, J, k;
        double kij;
        for (i=0; i<ndof; ++i)
        {
            ni = MAX_NDOFS*(en[i/ndn]) + i%ndn;
            li = m_fem.m_LCT[ni];
            for (j=0; j<ndof; ++j)
            {
                nj = MAX_NDOFS*(en[j/ndn]) + j%ndn;
                lj = m_fem.m_LCT[nj];
                
                if ((li >= 0) && (lj < 0))
                {
                    // dof i is constrained
                    FELinearConstraint& Li = *m_fem.m_LCA[li];
                    
                    assert(elm[i] == -1);
                    
                    list<FELinearConstraint::SlaveDOF>::iterator is = Li.slave.begin();
                    for (k=0; k<(int)Li.slave.size(); ++k, ++is)
                    {
                        I = is->neq;
                        J = elm[j];
                        kij = is->val*ke[i][j];
                        if ((J>=I) && (I >=0)) K.add(I,J, kij);
                        else
                        {
                            // adjust for prescribed dofs
                            J = -J-2;
                            if ((J>=0) && (I>=0)) m_Fd[I] -= kij*ui[J];
                        }
                    }
                }
                else if ((lj >= 0) && (li < 0))
                {
                    // dof j is constrained
                    FELinearConstraint& Lj = *m_fem.m_LCA[lj];
                    
                    assert(elm[j] == -1);
                    
                    list<FELinearConstraint::SlaveDOF>::iterator js = Lj.slave.begin();
                    
                    for (k=0; k<(int)Lj.slave.size(); ++k, ++js)
                    {
                        I = elm[i];
                        J = js->neq;
                        kij = js->val*ke[i][j];
                        if ((J>=I) && (I >=0)) K.add(I,J, kij);
                        else
                        {
                            // adjust for prescribed dofs
                            J = -J-2;
                            if ((J>=0) && (I>=0)) m_Fd[I] -= kij*ui[J];
                        }
                    }
                }
                else if ((li >= 0) && (lj >= 0))
                {
                    // both dof i and j are constrained
                    FELinearConstraint& Li = *m_fem.m_LCA[li];
                    FELinearConstraint& Lj = *m_fem.m_LCA[lj];
                    
                    list<FELinearConstraint::SlaveDOF>::iterator is = Li.slave.begin();
                    list<FELinearConstraint::SlaveDOF>::iterator js = Lj.slave.begin();
                    
                    assert(elm[i] == -1);
                    assert(elm[j] == -1);
                    
                    for (k=0; k<(int)Li.slave.size(); ++k, ++is)
                    {
                        js = Lj.slave.begin();
                        for  (l=0; l<(int)Lj.slave.size(); ++l, ++js)
                        {
                            I = is->neq;
                            J = js->neq;
                            kij = ke[i][j]*is->val*js->val;
                            
                            if ((J>=I) && (I >=0)) K.add(I,J, kij);
                            else
                            {
                                // adjust for prescribed dofs
                                J = -J-2;
                                if ((J>=0) && (I>=0)) m_Fd[I] -= kij*ui[J];
                            }
                        }
                    }
                }
            }
        }
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
        dom.InternalForces(RHS);
    }
    
    // calculate the body forces
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        for (int j=0; j<m_fem.BodyLoads(); ++j)
        {
            FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(m_fem.GetBodyLoad(j));
            dom.BodyForce(RHS, *pbf);
        }
    }
    
    // calculate inertial forces
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        dom.InertialForces(RHS);
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
