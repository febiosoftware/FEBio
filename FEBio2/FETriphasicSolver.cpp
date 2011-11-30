#include "stdafx.h"
#include "FETriphasicSolver.h"
#include "FETriphasicDomain.h"
#include "Interrupt.h"
#include "log.h"

#ifdef WIN32
#include <float.h>
#define ISNAN(x) _isnan(x)
#endif

#ifdef LINUX
#include <math.h>
#define ISNAN(x) isnan(x)
#endif

#ifdef __APPLE__
#include <math.h>
#define ISNAN(x) isnan(x)
#endif

//-----------------------------------------------------------------------------
FETriphasicSolver::FETriphasicSolver(FEM& fem) : FEBiphasicSolver(fem)
{
	m_Ctol = 0.01;
}

//-----------------------------------------------------------------------------
//! Allocates and initializes the data structures.
//
bool FETriphasicSolver::Init()
{
	// initialize base class
	if (FEBiphasicSolver::Init() == false) return false;
	
	// allocate concentration-vectors
	assert ((m_nceq[0] > 0) || (m_nceq[1] > 0));
	m_ci[0].assign(m_nceq[0], 0);
	m_ci[1].assign(m_nceq[1], 0);
	m_Ci[0].assign(m_nceq[0], 0);
	m_Ci[1].assign(m_nceq[1], 0);
	
	// we need to fill the total displacement vector m_Ut
	// TODO: I need to find an easier way to do this
	FEMesh& mesh = m_fem.m_mesh;
	int i, n;
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		
		// concentration dofs
		n = node.m_ID[DOF_C]; if (n >= 0) m_Ut[n] = node.m_ct[0];
		n = node.m_ID[DOF_C+1]; if (n >= 0) m_Ut[n] = node.m_ct[1];
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Prepares the data for the first QN iteration. 
//!
void FETriphasicSolver::PrepStep(double time)
{
	zero(m_Ci[0]);
	zero(m_Ci[1]);
	FEBiphasicSolver::PrepStep(time);
}

//-----------------------------------------------------------------------------
//! Implements the BFGS algorithm to solve the nonlinear FE equations.
//! The details of this implementation of the BFGS method can be found in:
//!   "Finite Element Procedures", K.J. Bathe, p759 and following
//!
bool FETriphasicSolver::Quasin(double time)
{
	int i;
	double s;
	
	// convergence norms
	double	normR1;		// residual norm
	double	normE1;		// energy norm
	double	normU;		// displacement norm
	double	normu;		// displacement increment norm
	double	normRi;		// initial residual norm
	double	normEi;		// initial energy norm
	double	normEm;		// max energy norm
	double	normUi;		// initial displacement norm
	
	// poro convergence norms data
	double	normPi;		// initial pressure norm
	double	normP;		// current pressure norm
	double	normp;		// incremement pressure norm
	
	// solute convergence data
	double	normCi[2];	// initial concentration norm
	double	normC[2];	// current concentration norm
	double	normc[2];	// incremement concentration norm
	
	// initialize flags
	bool bconv = false;		// convergence flag
	bool breform = false;	// reformation flag

	// get the current step
	FEAnalysisStep* pstep = dynamic_cast<FEAnalysisStep*>(m_fem.m_pStep);
	
	// make sure this is poro-solute problem
	assert(pstep->m_nModule == FE_TRIPHASIC);
	
	// prepare for the first iteration
	PrepStep(time);
	
	// check for CTRL+C interruption before we do any work
	if (m_fem.m_bInterruptable)
	{
		Interruption itr;
		if (itr.m_bsig)
		{
			itr.m_bsig = false;
			itr.interrupt();
		}
	}
	
	// calculate initial stiffness matrix
	if (ReformStiffness() == false) return false;
	
	// calculate initial residual
	if (Residual(m_bfgs.m_R0) == false) return false;
	
	m_bfgs.m_R0 += m_Fd;
	
	// TODO: I can check here if the residual is zero.
	// If it is than there is probably no force acting on the system
	// if (m_R0*m_R0 < eps) bconv = true;
	
	//	double r0 = m_R0*m_R0;
	
	Logfile::MODE oldmode;
	
	clog.printf("\n===== beginning time step %d : %lg =====\n", m_fem.m_pStep->m_ntimesteps+1, m_fem.m_ftime);
	
	// loop until converged or when max nr of reformations reached
	do
	{
		oldmode = clog.GetMode();
		if ((m_fem.m_pStep->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
			(m_fem.m_pStep->GetPrintLevel() != FE_PRINT_NEVER)) clog.SetMode(Logfile::FILE_ONLY);
		
		clog.printf(" %d\n", m_niter+1);
		clog.SetMode(oldmode);
		
		// assume we'll converge. 
		bconv = true;
		
		// solve the equations
		m_SolverTime.start();
		{
			m_bfgs.SolveEquations(m_bfgs.m_ui, m_bfgs.m_R0);
		}
		m_SolverTime.stop();
		
		// check for nans
		if (m_fem.GetDebugFlag())
		{
			double du = m_bfgs.m_ui*m_bfgs.m_ui;
			if (ISNAN(du)) throw NANDetected();
		}
		
		// set initial convergence norms
		if (m_niter == 0)
		{
			normRi = fabs(m_bfgs.m_R0*m_bfgs.m_R0);
			normEi = fabs(m_bfgs.m_ui*m_bfgs.m_R0);
			normUi = fabs(m_bfgs.m_ui*m_bfgs.m_ui);
			normEm = normEi;
		}
		
		// perform a linesearch
		// the geometry is also updated in the line search
		if (m_bfgs.m_LStol > 0) s = m_bfgs.LineSearch();
		else
		{
			s = 1;
			
			// Update geometry
			Update(m_bfgs.m_ui);
			
			// calculate residual at this point
			Residual(m_bfgs.m_R1);
		}
		
		// update total displacements
		for (i=0; i<m_neq; ++i) m_Ui[i] += s*m_bfgs.m_ui[i];
		
		// calculate norms
		normR1 = m_bfgs.m_R1*m_bfgs.m_R1;
		normu  = (m_bfgs.m_ui*m_bfgs.m_ui)*(s*s);
		normU  = m_Ui*m_Ui;
		normE1 = s*fabs(m_bfgs.m_ui*m_bfgs.m_R1);
		
		// check residual norm
		if ((m_Rtol > 0) && (normR1 > m_Rtol*normRi)) bconv = false;	
		
		// check displacement norm
		if ((m_Dtol > 0) && (normu  > (m_Dtol*m_Dtol)*normU )) bconv = false;
		
		// check energy norm
		if ((m_Etol > 0) && (normE1 > m_Etol*normEi)) bconv = false;
		
		// check linestep size
		if ((m_bfgs.m_LStol > 0) && (s < m_bfgs.m_LSmin)) bconv = false;
		
		// check energy divergence
		if (normE1 > normEm) bconv = false;
		
		// check triphasic convergence
		{
			// extract the pressure increments
			GetPressureData(m_pi, m_bfgs.m_ui);
			
			// set initial norm
			if (m_niter == 0) normPi = fabs(m_pi*m_pi);
			
			// update total pressure
			for (i=0; i<m_npeq; ++i) m_Pi[i] += s*m_pi[i];
			
			// calculate norms
			normP = m_Pi*m_Pi;
			normp = (m_pi*m_pi)*(s*s);
			
			// check convergence
			if ((m_Ptol > 0) && (normp > (m_Ptol*m_Ptol)*normP)) bconv = false;
		}
		
		// check solute convergence
		{
			// extract the concentration increments
			GetConcentrationData(m_ci[0], m_bfgs.m_ui,0);
			GetConcentrationData(m_ci[1], m_bfgs.m_ui,1);
			
			// set initial norm
			if (m_niter == 0) {
				normCi[0] = fabs(m_ci[0]*m_ci[0]);
				normCi[1] = fabs(m_ci[1]*m_ci[1]);
			}
			
			// update total concentration
			for (i=0; i<m_nceq[0]; ++i) m_Ci[0][i] += s*m_ci[0][i];
			for (i=0; i<m_nceq[1]; ++i) m_Ci[1][i] += s*m_ci[1][i];
			
			// calculate norms
			normC[0] = m_Ci[0]*m_Ci[0]; normC[1] = m_Ci[1]*m_Ci[1];
			normc[0] = (m_ci[0]*m_ci[0])*(s*s); normc[1] = (m_ci[1]*m_ci[1])*(s*s);
			
			// check convergence
			if ((m_Ctol > 0) && 
				((normc[0] > (m_Ctol*m_Ctol)*normC[0]) || 
				 (normc[1] > (m_Ctol*m_Ctol)*normC[1]))) bconv = false;
		}
		
		// print convergence summary
		oldmode = clog.GetMode();
		if ((m_fem.m_pStep->GetPrintLevel() <= FE_PRINT_MAJOR_ITRS) &&
			(m_fem.m_pStep->GetPrintLevel() != FE_PRINT_NEVER)) clog.SetMode(Logfile::FILE_ONLY);
		
		clog.printf(" Nonlinear solution status: time= %lg\n", time); 
		clog.printf("\tstiffness updates             = %d\n", m_bfgs.m_nups);
		clog.printf("\tright hand side evaluations   = %d\n", m_nrhs);
		clog.printf("\tstiffness matrix reformations = %d\n", m_nref);
		if (m_bfgs.m_LStol > 0) clog.printf("\tstep from line search         = %lf\n", s);
		clog.printf("\tconvergence norms :        INITIAL         CURRENT         REQUIRED\n");
		clog.printf("\t residual             %15le %15le %15le \n", normRi, normR1, m_Rtol*normRi);
		clog.printf("\t energy               %15le %15le %15le \n", normEi, normE1, m_Etol*normEi);
		clog.printf("\t displacement         %15le %15le %15le \n", normUi, normu ,(m_Dtol*m_Dtol)*normU );
		clog.printf("\t fluid pressure       %15le %15le %15le \n", normPi, normp ,(m_Ptol*m_Ptol)*normP );
		clog.printf("\t cation concentration %15le %15le %15le \n", normCi[0], normc[0] ,(m_Ctol*m_Ctol)*normC[0] );
		clog.printf("\t anion concentration  %15le %15le %15le \n", normCi[1], normc[1] ,(m_Ctol*m_Ctol)*normC[1] );
		
		clog.SetMode(oldmode);
		
		// check if we have converged. 
		// If not, calculate the BFGS update vectors
		if (bconv == false)
		{
			if ((normR1 < m_Rmin))
			{
				// check for almost zero-residual on the first iteration
				// this might be an indication that there is no force on the system
				clog.printbox("WARNING", "No force acting on the system.");
				bconv = true;
			}
			else if (s < m_bfgs.m_LSmin)
			{
				// check for zero linestep size
				clog.printbox("WARNING", "Zero linestep size. Stiffness matrix will now be reformed");
				breform = true;
			}
			else if (normE1 > normEm)
			{
				// check for diverging
				clog.printbox("WARNING", "Problem is diverging. Stiffness matrix will now be reformed");
				normEm = normE1;
				normEi = normE1;
				normRi = normR1;
				normPi = normp;
				normCi[0] = normc[0]; normCi[1] = normc[1];
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
							clog.printbox("WARNING", "The BFGS update has failed.\nStiffness matrix will now be reformed.");
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
							clog.printbox("WARNING", "Max nr of iterations reached.\nStiffness matrix will now be reformed.");
						
					}
				}
			}	
			
			// zero displacement increments
			// we must set this to zero before the reformation
			// because we assume that the prescribed displacements are stored 
			// in the m_ui vector.
			zero(m_bfgs.m_ui);
			
			// reform stiffness matrices if necessary
			if (breform)
			{
				clog.printf("Reforming stiffness matrix: reformation #%d\n\n", m_nref);
				
				// reform the matrix
				if (ReformStiffness() == false) break;
				
				// reset reformation flag
				breform = false;
			}
			
			// copy last calculated residual
			m_bfgs.m_R0 = m_bfgs.m_R1;
		}
		else if (pstep->m_baugment)
		{
			// we have converged, so let's see if the augmentations have converged as well
			
			clog.printf("\n........................ augmentation # %d\n", m_naug+1);
			
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
				Residual(m_bfgs.m_R0);
				
				// reform the matrix if we are using full-Newton
				FEAnalysisStep* pstep = dynamic_cast<FEAnalysisStep*>(m_fem.m_pStep);
				if (pstep->m_psolver->m_bfgs.m_maxups == 0)
				{
					clog.printf("Reforming stiffness matrix: reformation #%d\n\n", m_nref);
					if (ReformStiffness() == false) break;
				}
			}
		}
		
		// increase iteration number
		m_niter++;
		
		// let's flush the logfile to make sure the last output will not get lost
		clog.flush();
		
		// check for CTRL+C interruption
		if (m_fem.m_bInterruptable)
		{
			Interruption itr;
			if (itr.m_bsig)
			{
				itr.m_bsig = false;
				itr.interrupt();
			}
		}
	}
	while (bconv == false);
	
	// when converged, 
	// print a convergence summary to the clog file
	if (bconv)
	{
		Logfile::MODE mode = clog.SetMode(Logfile::FILE_ONLY);
		if (mode != Logfile::NEVER)
		{
			clog.printf("\nconvergence summary\n");
			clog.printf("    number of iterations   : %d\n", m_niter);
			clog.printf("    number of reformations : %d\n", m_nref);
		}
		clog.SetMode(mode);
	}
	
	// if converged we update the total displacements
	if (bconv)
	{
		m_Ut += m_Ui;
	}
	
	return bconv;
}

//-----------------------------------------------------------------------------
void FETriphasicSolver::GetConcentrationData(vector<double> &ci, vector<double> &ui, const int ion)
{
	int N = m_fem.m_mesh.Nodes(), nid, m = 0;
	zero(ci);
	for (int i=0; i<N; ++i)
	{
		FENode& n = m_fem.m_mesh.Node(i);
		nid = n.m_ID[DOF_C+ion];
		if (nid != -1)
		{
			nid = (nid < -1 ? -nid-2 : nid);
			ci[m++] = ui[nid];
			assert(m <= (int) ci.size());
		}
	}
}


//-----------------------------------------------------------------------------
//! Save data to dump file

void FETriphasicSolver::Serialize(DumpFile& ar)
{
	FEBiphasicSolver::Serialize(ar);
	
	if (ar.IsSaving())
	{
		ar << m_Ctol;
	}
	else
	{
		ar >> m_Ctol;
	}
}
