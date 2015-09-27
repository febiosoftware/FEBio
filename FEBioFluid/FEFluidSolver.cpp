#include "stdafx.h"
#include "FEFluidSolver.h"
#include "FEFluidDomain.h"
#include "FEFluidResidualVector.h"
#include "FECore/FENodeReorder.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include "NumCore/NumCore.h"
#include <assert.h>

#ifdef WIN32
#include <float.h>
#define ISNAN(x) _isnan(x)
#endif

#ifdef LINUX
#define ISNAN(x) std::isnan(x)
#endif

#ifdef __APPLE__
#include <math.h>
#define ISNAN(x) isnan(x)
#endif

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_PARAMETER_LIST(FEFluidSolver, FESolver)
ADD_PARAMETER(m_Vtol         , FE_PARAM_DOUBLE, "vtol"        );
ADD_PARAMETER(m_Rtol         , FE_PARAM_DOUBLE, "rtol"        );
ADD_PARAMETER(m_Rmin         , FE_PARAM_DOUBLE, "min_residual");
ADD_PARAMETER(m_LStol , FE_PARAM_DOUBLE, "lstol"       );
ADD_PARAMETER(m_LSmin , FE_PARAM_DOUBLE, "lsmin"       );
ADD_PARAMETER(m_LSiter, FE_PARAM_INT   , "lsiter"      );
ADD_PARAMETER(m_bfgs.m_maxref, FE_PARAM_INT   , "max_refs"    );
ADD_PARAMETER(m_bfgs.m_maxups, FE_PARAM_INT   , "max_ups"     );
ADD_PARAMETER(m_bfgs.m_cmax  , FE_PARAM_DOUBLE, "cmax"        );
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
    m_Rtol = 0;	// deactivate residual convergence
    m_Vtol = 0.001;
    m_Rmin = 1.0e-20;
    
    m_niter = 0;
    
    m_pK = 0;
    m_neq = 0;
    m_plinsolve = 0;
    
    m_bsymm = false;
    m_bdivreform = true;
    m_bdoreforms = true;
    
}

//-----------------------------------------------------------------------------
FEFluidSolver::~FEFluidSolver()
{
    delete m_plinsolve;	// clean up linear solver data
    delete m_pK;		// clean up stiffnes matrix data
}

//-----------------------------------------------------------------------------
//! Clean
//! \todo Why can this not be done in destructor?
void FEFluidSolver::Clean()
{
    if (m_plinsolve) m_plinsolve->Destroy();
}

//-----------------------------------------------------------------------------
//! Allocates and initializes the data structures used by the FEFluidSolver
//
bool FEFluidSolver::Init()
{
    // check parameters
    if (m_Vtol <  0.0) { felog.printf("Error: vtol must be nonnegative.\n"   ); return false; }
    if (m_Rtol <  0.0) { felog.printf("Error: rtol must be nonnegative.\n"); return false; }
    if (m_Rmin <  0.0) { felog.printf("Error: min_residual must be nonnegative.\n"  ); return false; }
    if (m_LStol  < 0.0) { felog.printf("Error: lstol must be nonnegative.\n" ); return false; }
    if (m_LSmin  < 0.0) { felog.printf("Error: lsmin must be nonnegative.\n" ); return false; }
    if (m_LSiter < 0) { felog.printf("Error: lsiter must be nonnegative.\n"  ); return false; }
    if (m_bfgs.m_maxref < 0) { felog.printf("Error: max_refs must be nonnegative.\n"); return false; }
    if (m_bfgs.m_maxups < 0) { felog.printf("Error: max_ups must be nonnegative.\n" ); return false; }
    if (m_bfgs.m_cmax   < 0) { felog.printf("Error: cmax must be nonnegative.\n"    ); return false; }
    
    // Now that we have determined the equation numbers we can continue
    // with creating the stiffness matrix. First we select the linear solver
    // The stiffness matrix is created in CreateStiffness
    // Note that if a particular solver was requested in the input file
    // then the solver might already be allocated. That's way we need to check it.
    if (m_plinsolve == 0)
    {
        m_plinsolve = NumCore::CreateLinearSolver(m_fem.m_nsolver);
        if (m_plinsolve == 0)
        {
            felog.printbox("FATAL ERROR","Unknown solver type selected\n");
            return false;
        }
    }
    
    // allocate storage for the sparse matrix that will hold the stiffness matrix data
    // we let the solver allocate the correct type of matrix format
    SparseMatrix* pS = m_plinsolve->CreateSparseMatrix(m_bsymm? SPARSE_SYMMETRIC : SPARSE_UNSYMMETRIC);
    if (pS == 0)
    {
        felog.printbox("FATAL ERROR", "The selected linear solver does not support the requested\n matrix format.\nPlease select a different linear solver.\n");
        return false;
    }
    
    // clean up the stiffness matrix if we have one
    if (m_pK) delete m_pK; m_pK = 0;
    
    // Create the stiffness matrix.
    // Note that this does not construct the stiffness matrix. This
    // is done later in the StiffnessMatrix routine.
    m_pK = new FEStiffnessMatrix(pS);
    if (m_pK == 0)
    {
        felog.printbox("FATAL ERROR", "Failed allocating stiffness matrix\n\n");
        return false;
    }
    
    // get nr of equations
    int neq = m_neq;
    
    // allocate vectors
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
        n = node.m_ID[DOF_VX]; if (n >= 0) m_Vt[n] = node.m_vt.x;
        n = node.m_ID[DOF_VY]; if (n >= 0) m_Vt[n] = node.m_vt.y;
        n = node.m_ID[DOF_VZ]; if (n >= 0) m_Vt[n] = node.m_vt.z;
        
        // dilatation dofs
        n = node.m_ID[DOF_E]; if (n >= 0) m_Vt[n] = node.m_et;
    }
    
    // initialize BFGS data
    m_bfgs.Init(neq, this, m_plinsolve);
    
    // set the create stiffness matrix flag
    m_breshape = true;
    
    return true;
}

//-----------------------------------------------------------------------------
//! Save data to dump file

void FEFluidSolver::Serialize(DumpFile& ar)
{
    // Serialize parameters
    FESolver::Serialize(ar);
    
    if (ar.IsSaving())
    {
        ar << m_Vtol << m_Rtol << m_Rmin;
        ar << m_bsymm;
        ar << m_nrhs;
        ar << m_niter;
        ar << m_nref << m_ntotref;
        ar << m_naug;
        ar << m_neq;
        
        ar << m_LStol << m_LSiter << m_LSmin;
        ar << m_bfgs.m_maxups;
        ar << m_bfgs.m_maxref;
        ar << m_bfgs.m_cmax;
        ar << m_bfgs.m_nups;
    }
    else
    {
        ar >> m_Vtol >> m_Rtol >> m_Rmin;
        ar >> m_bsymm;
        ar >> m_nrhs;
        ar >> m_niter;
        ar >> m_nref >> m_ntotref;
        ar >> m_naug;
        ar >> m_neq;
        
        ar >> m_LStol >> m_LSiter >> m_LSmin;
        ar >> m_bfgs.m_maxups;
        ar >> m_bfgs.m_maxref;
        ar >> m_bfgs.m_cmax;
        ar >> m_bfgs.m_nups;
    }
}

//-----------------------------------------------------------------------------
//! Determine the number of linear equations and assign equation numbers
//!

//-----------------------------------------------------------------------------
//!	This function initializes the equation system.
//! It is assumed that all free dofs up until now have been given an ID >= 0
//! and the fixed or rigid dofs an ID < 0.
//! After this operation the nodal ID array will contain the equation
//! number assigned to the corresponding degree of freedom. To distinguish
//! between free or unconstrained dofs and constrained ones the following rules
//! apply to the ID array:
//!
//!           /
//!          |  >=  0 --> dof j of node i is a free dof
//! ID[i][j] <  == -1 --> dof j of node i is a fixed (no equation assigned too)
//!          |  <  -1 --> dof j of node i is constrained and has equation nr = -ID[i][j]-2
//!           \
//!
bool FEFluidSolver::InitEquations()
{
    int i, j;
    
    // get the mesh
    FEMesh& mesh = m_fem.GetMesh();
    
    // initialize nr of equations
    int neq = 0;
    
    // see if we need to optimize the bandwidth
    if (m_fem.m_bwopt)
    {
        // reorder the node numbers
        vector<int> P(mesh.Nodes());
        FENodeReorder mod;
        mod.Apply(mesh, P);
        
        // set the equation numbers
        for (i=0; i<mesh.Nodes(); ++i)
        {
            FENode& node = mesh.Node(P[i]);
            for (j=0; j<(int)node.m_ID.size(); ++j)
            {
                if      (node.m_ID[j] == DOF_FIXED     ) { node.m_ID[j] = -1; }
                else if (node.m_ID[j] == DOF_OPEN      ) { node.m_ID[j] =  neq++; }
                else if (node.m_ID[j] == DOF_PRESCRIBED) { node.m_ID[j] = -neq-2; neq++; }
                else { assert(false); return false; }
            }
        }
    }
    else
    {
        // give all free dofs an equation number
        for (i=0; i<mesh.Nodes(); ++i)
        {
            FENode& node = mesh.Node(i);
            for (j=0; j<(int)node.m_ID.size(); ++j)
            {
                if      (node.m_ID[j] == DOF_FIXED     ) { node.m_ID[j] = -1; }
                else if (node.m_ID[j] == DOF_OPEN      ) { node.m_ID[j] =  neq++; }
                else if (node.m_ID[j] == DOF_PRESCRIBED) { node.m_ID[j] = -neq-2; neq++; }
                else { assert(false); return false; }
            }
        }
    }
    
    // store the number of equations
    m_neq = neq;
    
    // All initialization is done
    return true;
}

//-----------------------------------------------------------------------------
//!  Creates the global stiffness matrix
//! \todo Can we move this to the FEStiffnessMatrix::Create function?
bool FEFluidSolver::CreateStiffness(bool breset)
{
    // clean up the solver
    if (m_pK->NonZeroes()) m_plinsolve->Destroy();
    
    // clean up the stiffness matrix
    m_pK->Clear();
    
    // create the stiffness matrix
    felog.printf("===== reforming stiffness matrix:\n");
    if (m_pK->Create(&GetFEModel(), m_neq, breset) == false)
    {
        felog.printf("FATAL ERROR: An error occured while building the stiffness matrix\n\n");
        return false;
    }
    else
    {
        // output some information about the direct linear solver
        int neq = m_pK->Rows();
        int nnz = m_pK->NonZeroes();
        felog.printf("\tNr of equations ........................... : %d\n", neq);
        felog.printf("\tNr of nonzeroes in stiffness matrix ....... : %d\n", nnz);
        felog.printf("\n");
    }
    // let's flush the logfile to make sure the last output will not get lost
    felog.flush();
    
    // Do the preprocessing of the solver
    m_SolverTime.start();
    {
        if (!m_plinsolve->PreProcess()) throw FatalError();
    }
    m_SolverTime.stop();
    
    // done!
    return true;
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
        if ((n = node.m_ID[DOF_VX]) >= 0) node.m_vt.x = m_Vt[n] + m_Vi[n] + vi[n];
        if ((n = node.m_ID[DOF_VY]) >= 0) node.m_vt.y = m_Vt[n] + m_Vi[n] + vi[n];
        if ((n = node.m_ID[DOF_VZ]) >= 0) node.m_vt.z = m_Vt[n] + m_Vi[n] + vi[n];
        
        // dilatation dofs
        // current dilatation = total at prev conv step + total increment so far + current increment
        if ((n = node.m_ID[DOF_E]) >= 0) node.m_et = m_Vt[n] + m_Vi[n] + vi[n];
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
                    case 0: d += si->val*node.m_vt.x; break;
                    case 1: d += si->val*node.m_vt.y; break;
                    case 2: d += si->val*node.m_vt.z; break;
                    case 3: d += si->val*node.m_et;   break;
                }
            }
            
            switch (lc.master.bc)
            {
                case 0: node.m_vt.x = d; break;
                case 1: node.m_vt.y = d; break;
                case 2: node.m_vt.z = d; break;
                case 3: node.m_et   = d; break;
            }
        }
    }
    
}

//-----------------------------------------------------------------------------
//! Updates the current state of the model
void FEFluidSolver::Update(vector<double>& ui)
{
    // update kinematics
    UpdateKinematics(ui);
    
    // update element stresses
    UpdateStresses();
    
    // write the new state
    m_fem.Write(FE_UNCONVERGED);
}

//-----------------------------------------------------------------------------
//!  Updates the element stresses
void FEFluidSolver::UpdateStresses()
{
    FEMesh& mesh = m_fem.GetMesh();
    
    // update the stresses on all domains
    for (int i=0; i<mesh.Domains(); ++i)
    {
        FEFluidDomain& dom = dynamic_cast<FEFluidDomain&>(mesh.Domain(i));
        dom.UpdateStresses(m_fem);
    }
}

//-----------------------------------------------------------------------------
//!  This function mainly calls the Quasin routine
//!  and deals with exceptions that require the immediate termination of
//!	quasi-Newton iterations.
bool FEFluidSolver::SolveStep(double time)
{
    bool bret;
    
    try
    {
        // let's try to call Quasin
        bret = Quasin(time);
    }
    catch (NegativeJacobian e)
    {
        // A negative jacobian was detected
        felog.printbox("ERROR","Negative jacobian was detected at element %d at gauss point %d\njacobian = %lg\n", e.m_iel, e.m_ng+1, e.m_vol);
        return false;
    }
    catch (MaxStiffnessReformations)
    {
        // max nr of reformations is reached
        felog.printbox("ERROR", "Max nr of reformations reached.");
        return false;
    }
    catch (ForceConversion)
    {
        // user forced conversion of problem
        felog.printbox("WARNING", "User forced conversion.\nSolution might not be stable.");
        return true;
    }
    catch (IterationFailure)
    {
        // user caused a forced iteration failure
        felog.printbox("WARNING", "User forced iteration failure.");
        return false;
    }
    catch (ZeroLinestepSize)
    {
        // a zero line step size was detected
        felog.printbox("ERROR", "Zero line step size.");
        return false;
    }
    catch (EnergyDiverging)
    {
        // problem was diverging after stiffness reformation
        felog.printbox("ERROR", "Problem diverging uncontrollably.");
        return false;
    }
    catch (FEMultiScaleException)
    {
        // the RVE problem didn't solve
        felog.printbox("ERROR", "The RVE problem has failed. Aborting macro run.");
        return false;
    }
    catch (DoRunningRestart)
    {
        // a request to fail the iteration and restart the time step
        return false;
    }
    
    return bret;
}

//-----------------------------------------------------------------------------
//! Prepares the data for the first BFGS-iteration.
void FEFluidSolver::PrepStep(double time)
{
    // initialize counters
    m_niter = 0;	// nr of iterations
    m_nrhs  = 0;	// nr of RHS evaluations
    m_nref  = 0;	// nr of stiffness reformations
    m_ntotref = 0;
    m_bfgs.m_nups	= 0;	// nr of stiffness updates between reformations
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
        ni.m_vp = ni.m_vt;
        ni.m_ap = ni.m_at;
    }
    
    // TODO: Pass this parameter to this function instead of time
    FETimePoint tp = m_fem.GetTime();
    
    // apply concentrated nodal forces
    // since these forces do not depend on the geometry
    // we can do this once outside the NR loop.
    NodalForces(m_Fn, tp);
    
    // apply prescribed velocities
    // we save the prescribed velocity increments in the ui vector
    vector<double>& ui = m_bfgs.m_ui;
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
    FEMaterialPoint::dt = m_fem.GetCurrentStep()->m_dt;
    FEMaterialPoint::time = m_fem.m_ftime;
    
    for (int i=0; i<mesh.Domains(); ++i) mesh.Domain(i).InitElements();
    
    // update stresses
    UpdateStresses();
}

//-----------------------------------------------------------------------------
//! Implements the BFGS algorithm to solve the nonlinear FE equations.
//! The details of this implementation of the BFGS method can be found in:
//!   "Finite Element Procedures", K.J. Bathe, p759 and following
bool FEFluidSolver::Quasin(double time)
{
    int i;
    
    vector<double> u0(m_neq);
    vector<double> Rold(m_neq);
    
    // convergence norms
    double	normR1;		// residual norm
    double	normU;		// velocity norm
    double	normu;		// velocity increment norm
    double	normRi = 0;		// initial residual norm
    double	normUi = 0;		// initial velocity norm
    
    // initialize flags
    bool bconv = false;		// convergence flag
    bool breform = false;	// reformation flag
    
    // Get the current step
    FEAnalysis* pstep = m_fem.GetCurrentStep();
    
    // prepare for the first iteration
    PrepStep(time);
    
    // do minor iterations callbacks
    m_fem.DoCallback(CB_MINOR_ITERS);
    
    // calculate initial stiffness matrix
    if (ReformStiffness() == false) return false;
    
    // calculate initial residual
    if (Residual(m_bfgs.m_R0) == false) return false;
    
    m_bfgs.m_R0 += m_Fd;
    
    // TODO: I can check here if the residual is zero.
    // If it is than there is probably no force acting on the system
    // if (m_R0*m_R0 < eps) bconv = true;
    
    //	double r0 = m_R0*m_R0;
    
    felog.printf("\n===== beginning time step %d : %lg =====\n", pstep->m_ntimesteps+1, m_fem.m_ftime);
    
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
            m_bfgs.SolveEquations(m_bfgs.m_ui, m_bfgs.m_R0);
        }
        m_SolverTime.stop();
        
        // set initial convergence norms
        if (m_niter == 0)
        {
            normRi = fabs(m_bfgs.m_R0*m_bfgs.m_R0);
            normUi = fabs(m_bfgs.m_ui*m_bfgs.m_ui);
        }
        
        // perform a linesearch
        // the geometry is also updated in the line search
        if (m_LStol > 0) s = LineSearch(1.0);
        else
        {
            s = 1;
            
            // Update geometry
            Update(m_bfgs.m_ui);
            
            // calculate residual at this point
            Residual(m_bfgs.m_R1);
        }
        
        // calculate norms
        normR1 = m_bfgs.m_R1*m_bfgs.m_R1;
        normu  = (m_bfgs.m_ui*m_bfgs.m_ui)*(s*s);
        
        // check for nans
        if (ISNAN(normR1) || ISNAN(normu)) throw NANDetected();
        
        // update total velocities
        int neq = (int)m_Vi.size();
        for (i=0; i<neq; ++i) m_Vi[i] += s*m_bfgs.m_ui[i];
        normU  = m_Vi*m_Vi;
        
        // check residual norm
        if ((m_Rtol > 0) && (normR1 > m_Rtol*normRi)) bconv = false;
        
        // check velocity norm
        if ((m_Vtol > 0) && (normu  > (m_Vtol*m_Vtol)*normU )) bconv = false;
        
        // check linestep size
        if ((m_LStol > 0) && (s < m_LSmin)) bconv = false;
        
        // print convergence summary
        oldmode = felog.GetMode();
        if ((pstep->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
            (pstep->GetPrintLevel() != FE_PRINT_NEVER)) felog.SetMode(Logfile::FILE_ONLY);
        
        felog.printf(" Nonlinear solution status: time= %lg\n", time);
        felog.printf("\tstiffness updates             = %d\n", m_bfgs.m_nups);
        felog.printf("\tright hand side evaluations   = %d\n", m_nrhs);
        felog.printf("\tstiffness matrix reformations = %d\n", m_nref);
        if (m_LStol > 0) felog.printf("\tstep from line search         = %lf\n", s);
        felog.printf("\tconvergence norms :     INITIAL         CURRENT         REQUIRED\n");
        felog.printf("\t   residual         %15le %15le %15le \n", normRi, normR1, m_Rtol*normRi);
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
            else
            {
                // If we havn't reached max nr of BFGS updates
                // do an update
                if (!breform)
                {
                    if (m_bfgs.m_nups < m_bfgs.m_maxups-1)
                    {
                        if (m_bfgs.Update(s, m_bfgs.m_ui, m_bfgs.m_R0, m_bfgs.m_R1) == false)
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
                        if (m_bfgs.m_maxups > 0)
                            felog.printbox("WARNING", "Max nr of iterations reached.\nStiffness matrix will now be reformed.");
                        
                    }
                }
            }
            
            // zero velocity increments
            // we must set this to zero before the reformation
            // because we assume that the prescribed velocities are stored
            // in the m_ui vector.
            zero(m_bfgs.m_ui);
            
            // reform stiffness matrices if necessary
            if (breform && m_bdoreforms)
            {
                felog.printf("Reforming stiffness matrix: reformation #%d\n\n", m_nref);
                
                // reform the matrix
                if (ReformStiffness() == false) break;
                
                // reset reformation flag
                breform = false;
            }
            
            // copy last calculated residual
            m_bfgs.m_R0 = m_bfgs.m_R1;
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
//! Reforms a stiffness matrix and factorizes it
bool FEFluidSolver::ReformStiffness()
{
    // first, let's make sure we have not reached the max nr of reformations allowed
    if (m_nref >= m_bfgs.m_maxref) throw MaxStiffnessReformations();
    
    // recalculate the shape of the stiffness matrix if necessary
    if (m_breshape)
    {
        // TODO: I don't think I need to update here
        //		if (m_fem.m_bcontact) UpdateContact();
        
        // reshape the stiffness matrix
        if (!CreateStiffness(m_niter == 0)) return false;
        
        // reset reshape flag, except for contact
        m_breshape = (m_fem.SurfacePairInteractions() > 0? true : false);
    }
    
    // calculate the global stiffness matrix
    FETimePoint tp = m_fem.GetTime();
    bool bret = StiffnessMatrix(tp);
    
    if (bret)
    {
        m_SolverTime.start();
        {
            // factorize the stiffness matrix
            m_plinsolve->Factor();
        }
        m_SolverTime.stop();
        
        // increase total nr of reformations
        m_nref++;
        m_ntotref++;
        
        // reset bfgs update counter
        m_bfgs.m_nups = 0;
    }
    
    return bret;
}

//-----------------------------------------------------------------------------
//! Calculates global stiffness matrix.

bool FEFluidSolver::StiffnessMatrix(const FETimePoint& tp)
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
    
    // let's check the stiffness matrix for zero diagonal elements
/*    int neq = K.Size();
    for (i=0; i<neq; ++i)
    {
        if (fabs(K.diag(i)) < 1e-15) throw ZeroDiagonal(-1, -1);
    }*/
    
    return true;
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
    DOFS& fedofs = *DOFS::GetInstance();
    int MAX_NDOFS = fedofs.GetNDOFS();
    
    // assemble into global stiffness matrix
    m_pK->Assemble(ke, elm);
    
    vector<double>& ui = m_bfgs.m_ui;
    
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
        if (psl->IsActive()) psl->Residual(RHS);
    }
    
    // get the time information
    FETimePoint tp = m_fem.GetTime();
    
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
        if ((n = -node.m_ID[DOF_VX]-2) >= 0) node.m_Fr.x = -m_Fr[n];
        if ((n = -node.m_ID[DOF_VY]-2) >= 0) node.m_Fr.y = -m_Fr[n];
        if ((n = -node.m_ID[DOF_VZ]-2) >= 0) node.m_Fr.z = -m_Fr[n];
    }
    
    // increase RHS counter
    m_nrhs++;
    
    return true;
}

//-----------------------------------------------------------------------------
//! calculates the concentrated nodal forces

void FEFluidSolver::NodalForces(vector<double>& F, const FETimePoint& tp)
{
    // zero nodal force vector
    zero(F);
    
    // loop over nodal loads
    int NNL = m_fem.NodalLoads();
    for (int i=0; i<NNL; ++i)
    {
        FENodalLoad& fc = *m_fem.NodalLoad(i);
        if (fc.IsActive())
        {
            int nid = fc.m_node;	// node ID
            int dof = fc.m_bc;		// degree of freedom
            
            // get the nodal load value
            double f = fc.Value();
            
            // assemble into residual
            AssembleResidual(nid, dof, f, F);
        }
    }
}
