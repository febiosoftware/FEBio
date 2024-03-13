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
#include "FEBiphasicSolver.h"
#include "FEBiphasicDomain.h"
#include <FEBioMech/FESlidingElasticInterface.h>
#include "FESlidingInterface2.h"
#include "FESlidingInterface3.h"
#include "FESlidingInterfaceBiphasic.h"
#include "FESlidingInterfaceBiphasicMixed.h"
#include <FEBioMech/FEElasticDomain.h>
#include <FEBioMech/FEPressureLoad.h>
#include <FEBioMech/FEResidualVector.h>
#include <FEBioMech/FESolidLinearSystem.h>
#include <FEBioMech/FESSIShellDomain.h>
#include <FEBioMech/FERigidConnector.h>
#include <FECore/log.h>
#include <FECore/sys.h>
#include <FECore/FEModel.h>
#include <FECore/FEModelLoad.h>
#include <FECore/FENodalLoad.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FEBoundaryCondition.h>
#include <FECore/FENLConstraint.h>
#include <FECore/FELinearConstraintManager.h>
#include "FEBiphasicAnalysis.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_FECORE_CLASS(FEBiphasicSolver, FENewtonSolver)
	BEGIN_PARAM_GROUP("Nonlinear solver");	// make sure this matches FENewtonSolver. 
		ADD_PARAMETER(m_Dtol      , FE_RANGE_GREATER_OR_EQUAL(0.0), "dtol"        );
		ADD_PARAMETER(m_Etol      , FE_RANGE_GREATER_OR_EQUAL(0.0), "etol");
		ADD_PARAMETER(m_Rtol      , FE_RANGE_GREATER_OR_EQUAL(0.0), "rtol");
		ADD_PARAMETER(m_Ptol, "ptol"        );
		ADD_PARAMETER(m_biphasicFormulation, "mixed_formulation");
	END_PARAM_GROUP();

	// obsolete parameters that used to be inherited from FESolidSolver2
	ADD_PARAMETER(m_rhoi      , "rhoi"            )->SetFlags(FE_PARAM_HIDDEN);
	ADD_PARAMETER(m_alpha     , "alpha"           )->SetFlags(FE_PARAM_HIDDEN);
	ADD_PARAMETER(m_beta      , "beta"            )->SetFlags(FE_PARAM_HIDDEN);
	ADD_PARAMETER(m_gamma     , "gamma"           )->SetFlags(FE_PARAM_HIDDEN);
	ADD_PARAMETER(m_logSolve  , "logSolve"        )->SetFlags(FE_PARAM_HIDDEN);
	ADD_PARAMETER(m_arcLength , "arc_length"      )->SetFlags(FE_PARAM_HIDDEN);
	ADD_PARAMETER(m_al_scale  , "arc_length_scale")->SetFlags(FE_PARAM_HIDDEN);

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEBiphasicSolver::FEBiphasicSolver(FEModel* pfem) : FENewtonSolver(pfem), 
	m_dofU(pfem), m_dofV(pfem), m_dofRQ(pfem),
	m_dofSU(pfem), m_dofSV(pfem), m_dofSA(pfem),
	m_dofP(pfem), m_dofSP(pfem),
	m_rigidSolver(pfem)
{
	// default values
	m_Rtol = 0;	// deactivate residual convergence 
	m_Dtol = 0.001;
	m_Etol = 0.01;
	m_Ptol = 0.01;
	m_Rmin = 1.0e-20;
	m_Rmax = 0;	// not used if zero

	m_ndeq = 0;
	m_npeq = 0;
	m_nreq = 0;
	m_niter = 0;

	m_msymm = REAL_UNSYMMETRIC; // assume non-symmetric stiffness matrix by default

	// set default formulation (full shape functions)
	m_biphasicFormulation = 0;

	m_solutionNorm.push_back(ConvergenceInfo());

	// get pressure dof
	if (pfem)
	{
		m_dofP.AddDof(pfem->GetDOFIndex("p"));
		m_dofSP.AddDof(pfem->GetDOFIndex("q"));
		m_dofU.AddVariable("displacement");
		m_dofRQ.AddVariable("rigid rotation");
		m_dofV.AddVariable("velocity");
		m_dofSU.AddVariable("shell displacement");
		m_dofSV.AddVariable("shell velocity");
		m_dofSA.AddVariable("shell acceleration");
	}
}

//-----------------------------------------------------------------------------
FEBiphasicSolver::~FEBiphasicSolver() {}

//-----------------------------------------------------------------------------
//! Allocates and initializes the data structures.
//
bool FEBiphasicSolver::Init()
{
	// initialize base class
	if (FENewtonSolver::Init() == false) return false;

	FEModel& fem = *GetFEModel();

	// allocate vectors
	//	m_Fn.assign(m_neq, 0);
	m_Fr.assign(m_neq, 0);
	m_Ui.assign(m_neq, 0);
	m_Ut.assign(m_neq, 0);
	m_Uip.assign(m_neq, 0);

	// we need to fill the total displacement vector m_Ut
	FEMesh& mesh = fem.GetMesh();
	gather(m_Ut, mesh, m_dofU[0]);
	gather(m_Ut, mesh, m_dofU[1]);
	gather(m_Ut, mesh, m_dofU[2]);
	gather(m_Ut, mesh, m_dofSU[0]);
	gather(m_Ut, mesh, m_dofSU[1]);
	gather(m_Ut, mesh, m_dofSU[2]);

	SolverWarnings();

	// allocate poro-vectors
    assert((m_ndeq > 0) || (m_npeq > 0));
    m_di.assign(m_ndeq, 0);
	m_Di.assign(m_ndeq, 0);

	if (m_npeq > 0) {
		m_pi.assign(m_npeq, 0);
		m_Pi.assign(m_npeq, 0);

		// we need to fill the total displacement vector m_Ut
		// (displacements are already handled in base class)
		FEMesh& mesh = GetFEModel()->GetMesh();
		gather(m_Ut, mesh, m_dofP[0]);
        gather(m_Ut, mesh, m_dofSP[0]);
    }

	return true;
}

//! Generate warnings if needed
void FEBiphasicSolver::SolverWarnings()
{
	FEModel& fem = *GetFEModel();

	// Generate warning if rigid connectors are used with symmetric stiffness
	if (m_msymm == REAL_SYMMETRIC) {
		for (int i = 0; i < fem.NonlinearConstraints(); ++i)
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
			for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
			{
				FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));
				FESlidingElasticInterface* pbw = dynamic_cast<FESlidingElasticInterface*>(pci);
				if (pbw) {
					feLogWarning("The sliding-elastic contact algorithm runs better with a non-symmetric stiffness matrix.\nYou may set symmetric_stiffness flag to 0 in Control section.");
					break;
				}
			}
		}
	}
}

//! Initialize equations
bool FEBiphasicSolver::InitEquations()
{
	// define the solution variables for the Newton solver
	// Do this before calling base class!
	// TODO: Maybe I can get default values from the domains?
	int pressureOrder = (m_biphasicFormulation == 1 ? 1 : -1);
	AddSolutionVariable(&m_dofU, -1, "displacement", m_Dtol);
	AddSolutionVariable(&m_dofSU, -1, "shell displacement", m_Dtol);
	AddSolutionVariable(&m_dofP, pressureOrder, "pressure", m_Ptol);
	AddSolutionVariable(&m_dofSP, pressureOrder, "shell fluid pressure", m_Ptol);

	// set the interpolation orders
	DOFS& dofs = GetFEModel()->GetDOFS();
	int var_u = dofs.GetVariableIndex("displacement");
	int var_p = dofs.GetVariableIndex("fluid pressure");
	dofs.SetVariableInterpolationOrder(var_u, -1);
	dofs.SetVariableInterpolationOrder(var_p, pressureOrder);

	// First call the base class.
	// This will initialize all equation numbers, except the rigid body equation numbers
	if (FENewtonSolver::InitEquations2() == false) return false;

	// store the number of equations we currently have
	m_nreq = m_neq;

	// Next, we assign equation numbers to the rigid body degrees of freedom
	int neq = m_rigidSolver.InitEquations(m_neq);
	if (neq == -1) return false;
	else m_neq = neq;

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
	
	// determined the nr of pressure and concentration equations
	m_ndeq = m_npeq = 0;
	
	FEMesh& mesh = fem.GetMesh();
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& n = mesh.Node(i);
		if (n.m_ID[m_dofU[0]] != -1) m_ndeq++;
		if (n.m_ID[m_dofU[1]] != -1) m_ndeq++;
		if (n.m_ID[m_dofU[2]] != -1) m_ndeq++;
        if (n.m_ID[m_dofSU[0]] != -1) m_ndeq++;
        if (n.m_ID[m_dofSU[1]] != -1) m_ndeq++;
        if (n.m_ID[m_dofSU[2]] != -1) m_ndeq++;
        if (n.m_ID[m_dofP [0]] != -1) m_npeq++;
        if (n.m_ID[m_dofSP[0]] != -1) m_npeq++;
    }

	return true;
}

//-----------------------------------------------------------------------------
//! Prepares the data for the first QN iteration. 
//!
//! \todo There is some more stuff in the base method that 
//!       I need to move to this method, but since it will
//!       change the order of some operations I need to make
//!       sure it won't break anything
void FEBiphasicSolver::PrepStep()
{
	zero(m_Pi);
	zero(m_Di);

	// for pressure nodal loads we need to multiply the time step size
	FEModel& fem = *GetFEModel();
	for (int i = 0; i < fem.ModelLoads(); ++i)
	{
		FENodalDOFLoad* pl = dynamic_cast<FENodalDOFLoad*>(fem.ModelLoad(i));
		if (pl && pl->IsActive())
		{
			if ((pl->GetDOF() == m_dofP[0]) || (pl->GetDOF() == m_dofSP[0]))
			{
				pl->SetDtScale(true);
			}
		}
	}

	FETimeInfo& tp = fem.GetTime();
	double dt = tp.timeIncrement;
	tp.augmentation = 0;

	// zero total displacements
	zero(m_Ui);

	// store previous mesh state
	// we need them for velocity and acceleration calculations
	FEMesh& mesh = fem.GetMesh();
	for (int i = 0; i < mesh.Nodes(); ++i)
	{
		FENode& ni = mesh.Node(i);
		vec3d vs = (ni.m_rt - ni.m_rp)/dt;
		vec3d vq = (ni.m_dt - ni.m_dp)/dt;
		ni.m_rp = ni.m_rt;
		ni.m_vp = ni.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);
		ni.m_dp = ni.m_dt;
		ni.UpdateValues();

		ni.set_vec3d(m_dofV[0], m_dofV[1], m_dofV[2], vs);

		// solid shell
		ni.set_vec3d(m_dofSV[0], m_dofSV[1], m_dofSV[2], vs - vq);
	}

	// apply concentrated nodal forces
	// since these forces do not depend on the geometry
	// we can do this once outside the NR loop.
//	vector<double> dummy(m_neq, 0.0);
//	zero(m_Fn);
//	FEResidualVector Fn(*GetFEModel(), m_Fn, dummy);
//	NodalLoads(Fn, tp);

	// apply boundary conditions
	// we save the prescribed displacements increments in the ui vector
	vector<double>& ui = m_ui;
	zero(ui);
	int nbc = fem.BoundaryConditions();
	for (int i = 0; i < nbc; ++i)
	{
		FEBoundaryCondition& dc = *fem.BoundaryCondition(i);
		if (dc.IsActive()) dc.PrepStep(ui);
	}

	// do the linear constraints
	fem.GetLinearConstraintManager().PrepStep();

	// initialize rigid bodies
	m_rigidSolver.PrepStep(tp, ui);

	// intialize material point data
	// NOTE: do this before the stresses are updated
	// TODO: does it matter if the stresses are updated before
	//       the material point data is initialized
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		if (dom.IsActive()) dom.PreSolveUpdate(tp);
	}

	// update model state
	UpdateModel();

	for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	{
		FENLConstraint* plc = fem.NonlinearConstraint(i);
		if (plc && plc->IsActive()) plc->PrepStep();
	}

	// see if we need to do contact augmentations
	m_baugment = false;
	for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
	{
		FEContactInterface& ci = dynamic_cast<FEContactInterface&>(*fem.SurfacePairConstraint(i));
		if (ci.IsActive() && (ci.m_laugon == 1)) m_baugment = true;
	}

	// see if we have to do nonlinear constraint augmentations
	if (fem.NonlinearConstraints() != 0) m_baugment = true;
}

//-----------------------------------------------------------------------------
//! Implements the BFGS algorithm to solve the nonlinear FE equations.
bool FEBiphasicSolver::Quasin()
{
	// convergence norms
	double	normR1;		// residual norm
	double	normE1;		// energy norm
	double	normD;		// displacement norm
	double	normd;		// displacement increment norm
	double	normRi = 0; // initial residual norm
	double	normEi = 0; // initial energy norm
	double	normEm = 0;	// max energy norm
	double	normDi = 0;	// initial displacement norm

	// poro convergence norms data
	double	normPi = 0;	// initial pressure norm
	double	normP;		// current pressure norm
	double	normp;		// incremement pressure norm

	// prepare for the first iteration
	FEModel& fem = *GetFEModel();
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

		// extract the pressure increments
		GetDisplacementData(m_di, m_ui);

		// set initial convergence norms
		if (m_niter == 0)
		{
			normRi = fabs(m_R0*m_R0);
			normEi = fabs(m_ui*m_R0);
			normDi = fabs(m_di*m_di);
			normEm = normEi;

			m_residuNorm.norm0 = normRi;
			m_energyNorm.norm0 = normEi;
			m_solutionNorm[0].norm0 = normDi;
		}

		// update all degrees of freedom
		for (int i = 0; i<m_neq; ++i) m_Ui[i] += s*m_ui[i];

		// update displacements
		for (int i = 0; i<m_ndeq; ++i) m_Di[i] += s*m_di[i];

		// calculate norms
		normR1 = m_R1*m_R1;
		normd  = (m_di*m_di)*(s*s);
		normD  = m_Di*m_Di;
		normE1 = s*fabs(m_ui*m_R1);

		m_residuNorm.norm = normR1;
		m_energyNorm.norm = normR1;
		m_solutionNorm[0].norm = normd;

		// check residual norm
		if ((m_Rtol > 0) && (normR1 > m_Rtol*normRi)) bconv = false;	

		// check displacement norm
		if ((m_Dtol > 0) && (normd  > (m_Dtol*m_Dtol)*normD )) bconv = false;

		// check energy norm
		if ((m_Etol > 0) && (normE1 > m_Etol*normEi)) bconv = false;

		// check linestep size
		if ((m_lineSearch->m_LStol > 0) && (s < m_lineSearch->m_LSmin)) bconv = false;

		// check energy divergence
		if (m_bdivreform)
		{
			if (normE1 > normEm) bconv = false;
		}

		// check poroelastic convergence
		{
			// extract the pressure increments
			GetPressureData(m_pi, m_ui);

			// set initial norm
			if (m_niter == 0) normPi = fabs(m_pi*m_pi);

			// update total pressure
			for (int i = 0; i<m_npeq; ++i) m_Pi[i] += s*m_pi[i];

			// calculate norms
			normP = m_Pi*m_Pi;
			normp = (m_pi*m_pi)*(s*s);

			// check convergence
			if ((m_Ptol > 0) && (normp > (m_Ptol*m_Ptol)*normP)) bconv = false;
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
		feLog("\t   displacement     %15le %15le %15le \n", normDi, normd ,(m_Dtol*m_Dtol)*normD );
		feLog("\t   fluid pressure   %15le %15le %15le \n", normPi, normp ,(m_Ptol*m_Ptol)*normP );

		if ((bconv == false) && (normR1 < m_Rmin))
		{
			// check for almost zero-residual on the first iteration
			// this might be an indication that there is no force on the system
			feLogWarning("No force acting on the system.");
			bconv = true;
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
				normPi = normp;
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

	// if converged we update the total displacements
	if (bconv)
	{
		m_Ut += m_Ui;
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! calculates the residual vector
//! Note that the concentrated nodal forces are not calculated here.
//! This is because they do not depend on the geometry 
//! so we only calculate them once (in Quasin) and then add them here.

bool FEBiphasicSolver::Residual(vector<double>& R)
{
	// get the time information
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();

	// zero nodal reaction forces
	zero(m_Fr);

	// setup global RHS vector
    zero(R);
	FEResidualVector RHS(fem, R, m_Fr);

	// zero rigid body reaction forces
	m_rigidSolver.Residual();

	// calculate internal stress force
    InternalForces(RHS);

    // calculate nodal reaction forces
    for (int i = 0; i < m_neq; ++i) m_Fr[i] -= R[i];
    
    // calculate external forces
    ExternalForces(RHS);
    
	// increase RHS counter
	m_nrhs++;

	return true;
}

//-----------------------------------------------------------------------------
//! Calculates global stiffness matrix.

bool FEBiphasicSolver::StiffnessMatrix()
{
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();

	// get the mesh
	FEMesh& mesh = fem.GetMesh();

	// setup the linear system of equations
	FESolidLinearSystem LS(this, &m_rigidSolver, *m_pK, m_Fd, m_ui, (m_msymm == REAL_SYMMETRIC), m_alpha, m_nreq);

	// calculate the stiffness matrix for each domain
	FEAnalysis* pstep = fem.GetCurrentStep();
	bool bsymm = (m_msymm == REAL_SYMMETRIC);
	if (pstep->m_nanalysis == FEBiphasicAnalysis::STEADY_STATE)
	{
		for (int i=0; i<mesh.Domains(); ++i) 
		{
            // Biphasic analyses may include biphasic and elastic domains
			FEBiphasicDomain* pbdom = dynamic_cast<FEBiphasicDomain*>(&mesh.Domain(i));
			if (pbdom) pbdom->StiffnessMatrixSS(LS, bsymm);
            else
			{
				FEElasticDomain* pedom = dynamic_cast<FEElasticDomain*>(&mesh.Domain(i));
				if (pedom) pedom->StiffnessMatrix(LS);
			}
		}
	}
	else
	{
		for (int i=0; i<mesh.Domains(); ++i) 
		{
            // Biphasic analyses may include biphasic and elastic domains
			FEBiphasicDomain* pbdom = dynamic_cast<FEBiphasicDomain*>(&mesh.Domain(i));
			if (pbdom) pbdom->StiffnessMatrix(LS, bsymm);
            else 
			{
				FEElasticDomain* pedom = dynamic_cast<FEElasticDomain*>(&mesh.Domain(i));
				if (pedom) pedom->StiffnessMatrix(LS);
			}
		}
	}

	// calculate contact stiffness
	ContactStiffness(LS);

	// calculate stiffness matrices for surface loads
	int nml = fem.ModelLoads();
	for (int i=0; i<nml; ++i)
	{
		FEModelLoad* pml = fem.ModelLoad(i);
		if (pml->IsActive()) pml->StiffnessMatrix(LS);
	}

	// calculate nonlinear constraint stiffness
	// note that this is the contribution of the 
	// constrainst enforced with augmented lagrangian
	NonLinearConstraintStiffness(LS, tp);

	// add contributions from rigid bodies
	m_rigidSolver.StiffnessMatrix(*m_pK, tp);

	return true;
}

//! Update the model's kinematic data. 
void FEBiphasicSolver::UpdateKinematics(vector<double>& ui)
{
	// first update all solid-mechanics kinematics
	FEModel& fem = *GetFEModel();

	// get the mesh
	FEMesh& mesh = fem.GetMesh();

	// update rigid bodies
	m_rigidSolver.UpdateRigidBodies(m_Ui, ui);

	// total displacements
	vector<double> U(m_Ut.size());
	int U_size = (int)U.size();
#pragma omp parallel for
	for (int i = 0; i < U_size; ++i)
	{
		U[i] = ui[i] + m_Ui[i] + m_Ut[i];
	}

	// update flexible nodes
	// translational dofs
	scatter3(U, mesh, m_dofU[0], m_dofU[1], m_dofU[2]);
	// shell dofs
	scatter3(U, mesh, m_dofSU[0], m_dofSU[1], m_dofSU[2]);

	// make sure the boundary conditions are fullfilled
	int nbcs = fem.BoundaryConditions();
	for (int i = 0; i < nbcs; ++i)
	{
		FEBoundaryCondition& bc = *fem.BoundaryCondition(i);
		if (bc.IsActive()) bc.Update();
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
	int NN = mesh.Nodes();
	{
		for (int i = 0; i < NN; ++i)
		{
			FENode& node = mesh.Node(i);
			if (node.m_rid == -1) {
				node.m_rt = node.m_r0 + node.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2]);
			}
			node.m_dt = node.m_d0 + node.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2])
				- node.get_vec3d(m_dofSU[0], m_dofSU[1], m_dofSU[2]);
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

	// update poroelastic data
	UpdatePoro(ui);
}

//-----------------------------------------------------------------------------
//! Updates the poroelastic data
void FEBiphasicSolver::UpdatePoro(vector<double>& ui)
{
    int i, n;
    
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	double dt = fem.GetTime().timeIncrement;

	// update poro-elasticity data
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		// update nodal pressures
		n = node.m_ID[m_dofP[0]];
		if (n >= 0) node.set(m_dofP[0], 0 + m_Ut[n] + m_Ui[n] + ui[n]);
        n = node.m_ID[m_dofSP[0]];
        if (n >= 0) node.set(m_dofSP[0], 0 + m_Ut[n] + m_Ui[n] + ui[n]);
    }

	// update poro-elasticity data
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		// update velocities
		vec3d vt = (node.m_rt - node.m_rp) / dt;
		node.set_vec3d(m_dofV[0], m_dofV[1], m_dofV[2], vt);
	}
}

//! Updates the current state of the model
void FEBiphasicSolver::Update(vector<double>& ui)
{
	FEModel& fem = *GetFEModel();
	FETimeInfo& tp = fem.GetTime();
	tp.currentIteration = m_niter;

	// update EAS
	UpdateEAS(ui);
	UpdateIncrementsEAS(ui, true);

	// update kinematics
	UpdateKinematics(ui);

	// update domains 
	FEMesh& mesh = fem.GetMesh();
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		dom.IncrementalUpdate(ui, false);
	}

	// update model state
	UpdateModel();
}

//! Update EAS
void FEBiphasicSolver::UpdateEAS(vector<double>& ui)
{
	FEModel& fem = *GetFEModel();

	FEMesh& mesh = fem.GetMesh();

	// update EAS on shell domains
	for (int i = 0; i < mesh.Domains(); ++i) {
		FESSIShellDomain* sdom = dynamic_cast<FESSIShellDomain*>(&mesh.Domain(i));
		if (sdom && sdom->IsActive()) sdom->UpdateEAS(ui);
	}
}

//! Update EAS
void FEBiphasicSolver::UpdateIncrementsEAS(vector<double>& ui, const bool binc)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// update EAS on shell domains
	for (int i = 0; i < mesh.Domains(); ++i) {
		FESSIShellDomain* sdom = dynamic_cast<FESSIShellDomain*>(&mesh.Domain(i));
		if (sdom && sdom->IsActive()) sdom->UpdateIncrementsEAS(ui, binc);
	}
}

void FEBiphasicSolver::UpdateModel()
{
	// mark all free-draining surfaces
	FEModel& fem = *GetFEModel();
	for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));

		FESlidingInterface2* psi2 = dynamic_cast<FESlidingInterface2*>(pci);
		if (psi2) psi2->MarkFreeDraining();
		FESlidingInterface3* psi3 = dynamic_cast<FESlidingInterface3*>(pci);
		if (psi3) psi3->MarkAmbient();
        FESlidingInterfaceBiphasic* psib = dynamic_cast<FESlidingInterfaceBiphasic*>(pci);
        if (psib) psib->MarkFreeDraining();
		FESlidingInterfaceBiphasicMixed* psbm = dynamic_cast<FESlidingInterfaceBiphasicMixed*>(pci);
		if (psbm) psbm->MarkFreeDraining();
	}

	// Update all contact interfaces
	// NOTE: note that we call the base class version here, not the overridden one!!
	FENewtonSolver::UpdateModel();

	// set free-draining boundary conditions
	for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));

		FESlidingInterface2* psi2 = dynamic_cast<FESlidingInterface2*>(pci);
		if (psi2) psi2->SetFreeDraining();
		FESlidingInterface3* psi3 = dynamic_cast<FESlidingInterface3*>(pci);
		if (psi3) psi3->SetAmbient();
        FESlidingInterfaceBiphasic* psib = dynamic_cast<FESlidingInterfaceBiphasic*>(pci);
        if (psib) psib->SetFreeDraining();
		FESlidingInterfaceBiphasicMixed* psbm = dynamic_cast<FESlidingInterfaceBiphasicMixed*>(pci);
		if (psbm) psbm->SetFreeDraining();
	}
    
    // make sure the prescribed BCs (fluid pressure) are fullfilled
    int nbcs = fem.BoundaryConditions();
    for (int i = 0; i<nbcs; ++i)
    {
        FEBoundaryCondition& bc = *fem.BoundaryCondition(i);
        if (bc.IsActive()) bc.Repair();
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicSolver::GetDisplacementData(vector<double> &di, vector<double> &ui)
{
	FEModel& fem = *GetFEModel();
	int N = fem.GetMesh().Nodes(), nid, m = 0;
	zero(di);
	for (int i=0; i<N; ++i)
	{
		FENode& n = fem.GetMesh().Node(i);
		nid = n.m_ID[m_dofU[0]];
		if (nid != -1)
		{
			nid = (nid < -1 ? -nid-2 : nid);
			di[m++] = ui[nid];
			assert(m <= (int) di.size());
		}
		nid = n.m_ID[m_dofU[1]];
		if (nid != -1)
		{
			nid = (nid < -1 ? -nid-2 : nid);
			di[m++] = ui[nid];
			assert(m <= (int) di.size());
		}
		nid = n.m_ID[m_dofU[2]];
		if (nid != -1)
		{
			nid = (nid < -1 ? -nid-2 : nid);
			di[m++] = ui[nid];
			assert(m <= (int) di.size());
		}
        nid = n.m_ID[m_dofSU[0]];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            di[m++] = ui[nid];
            assert(m <= (int) di.size());
        }
        nid = n.m_ID[m_dofSU[1]];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            di[m++] = ui[nid];
            assert(m <= (int) di.size());
        }
        nid = n.m_ID[m_dofSU[2]];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            di[m++] = ui[nid];
            assert(m <= (int) di.size());
        }
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicSolver::GetPressureData(vector<double> &pi, vector<double> &ui)
{
	FEModel& fem = *GetFEModel();
	int N = fem.GetMesh().Nodes(), nid, m = 0;
	zero(pi);
	for (int i=0; i<N; ++i)
	{
		FENode& n = fem.GetMesh().Node(i);
		nid = n.m_ID[m_dofP[0]];
		if (nid != -1)
		{
			nid = (nid < -1 ? -nid-2 : nid);
			pi[m++] = ui[nid];
			assert(m <= (int) pi.size());
		}
        nid = n.m_ID[m_dofSP[0]];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            pi[m++] = ui[nid];
            assert(m <= (int) pi.size());
        }
    }
}

//-----------------------------------------------------------------------------
//! Save data to dump file

void FEBiphasicSolver::Serialize(DumpStream& ar)
{
	// Serialize parameters
	FENewtonSolver::Serialize(ar);

	ar& m_nrhs;
	ar& m_niter;
	ar& m_nref& m_ntotref;
	ar& m_naug;
	ar& m_nreq;

	ar& m_Ut& m_Ui;

	ar& m_arcLength;
	ar& m_al_scale;

	if (ar.IsLoading())
	{
//		m_Fn.assign(m_neq, 0);
		m_Fr.assign(m_neq, 0);
//		m_Ui.assign(m_neq, 0);
	}

	// serialize rigid solver
	m_rigidSolver.Serialize(ar);

	if (ar.IsShallow()) return;
	ar & m_Ptol & m_ndeq & m_npeq;
	ar & m_nceq;
	ar & m_di & m_Di & m_pi & m_Pi;
}

//-----------------------------------------------------------------------------
//! Internal forces
void FEBiphasicSolver::InternalForces(FEGlobalVector& RHS)
{
    // get the time information
    FEModel& fem = *GetFEModel();
    const FETimeInfo& tp = fem.GetTime();
    
    // get the mesh
    FEMesh& mesh = fem.GetMesh();
    
    // calculate internal stress force
    if (fem.GetCurrentStep()->m_nanalysis == FEBiphasicAnalysis::STEADY_STATE)
    {
        for (int i=0; i<mesh.Domains(); ++i)
        {
            FEBiphasicDomain* pdom = dynamic_cast<FEBiphasicDomain*>(&mesh.Domain(i));
            if (pdom) pdom->InternalForcesSS(RHS);
            else
            {
                FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
                dom.InternalForces(RHS);
            }
        }
    }
    else
    {
        for (int i=0; i<mesh.Domains(); ++i)
        {
            FEBiphasicDomain* pdom = dynamic_cast<FEBiphasicDomain*>(&mesh.Domain(i));
            if (pdom) pdom->InternalForces(RHS);
            else
            {
                FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
                dom.InternalForces(RHS);
            }
        }
    }
    
}

//-----------------------------------------------------------------------------
//! External forces
void FEBiphasicSolver::ExternalForces(FEGlobalVector& RHS)
{
    // get the time information
    FEModel& fem = *GetFEModel();
    const FETimeInfo& tp = fem.GetTime();
    
    // get the mesh
    FEMesh& mesh = fem.GetMesh();
    
    // add model loads
    int NML = fem.ModelLoads();
    for (int i=0; i<NML; ++i)
    {
        FEModelLoad& mli = *fem.ModelLoad(i);
        if (mli.IsActive()) mli.LoadVector(RHS);
    }
    
    // calculate contact forces
    ContactForces(RHS);
    
    // calculate nonlinear constraint forces
    // note that these are the linear constraints
    // enforced using the augmented lagrangian
    NonLinearConstraintForces(RHS, tp);
    
    // set the nodal reaction forces
    // TODO: Is this a good place to do this?
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        node.set_load(m_dofU[0], 0);
        node.set_load(m_dofU[1], 0);
        node.set_load(m_dofU[2], 0);
        node.set_load(m_dofP[0], 0);

        int n;
        if ((n = node.m_ID[m_dofU[0]]) >= 0) node.set_load(m_dofU[0], -m_Fr[n]);
        if ((n = -node.m_ID[m_dofU[0]] - 2) >= 0) node.set_load(m_dofU[0], -m_Fr[n]);
        
        if ((n = node.m_ID[m_dofU[1]]) >= 0) node.set_load(m_dofU[1], -m_Fr[n]);
        if ((n = -node.m_ID[m_dofU[1]] - 2) >= 0) node.set_load(m_dofU[1], -m_Fr[n]);
        
        if ((n = node.m_ID[m_dofU[2]]) >= 0) node.set_load(m_dofU[2], -m_Fr[n]);
        if ((n = -node.m_ID[m_dofU[2]] - 2) >= 0) node.set_load(m_dofU[2], -m_Fr[n]);
        
        if ((n = node.m_ID[m_dofP[0]]) >= 0) node.set_load(m_dofP[0], -m_Fr[n]);
        if ((n = -node.m_ID[m_dofP[0]] - 2) >= 0) node.set_load(m_dofP[0], -m_Fr[n]);
    }
}

//! Calculates the contact forces
void FEBiphasicSolver::ContactForces(FEGlobalVector& R)
{
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();
	for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));
		if (pci->IsActive()) pci->LoadVector(R, tp);
	}
}

//! This function calculates the contact stiffness matrix
void FEBiphasicSolver::ContactStiffness(FELinearSystem& LS)
{
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();
	for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));
		if (pci->IsActive()) pci->StiffnessMatrix(LS, tp);
	}
}

//! calculate the nonlinear constraint forces 
void FEBiphasicSolver::NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp)
{
	FEModel& fem = *GetFEModel();
	int N = fem.NonlinearConstraints();
	for (int i = 0; i < N; ++i)
	{
		FENLConstraint* plc = fem.NonlinearConstraint(i);
		if (plc->IsActive()) plc->LoadVector(R, tp);
	}
}

//! Calculate the stiffness contribution due to nonlinear constraints
void FEBiphasicSolver::NonLinearConstraintStiffness(FELinearSystem& LS, const FETimeInfo& tp)
{
	FEModel& fem = *GetFEModel();
	int N = fem.NonlinearConstraints();
	for (int i = 0; i < N; ++i)
	{
		FENLConstraint* plc = fem.NonlinearConstraint(i);
		if (plc->IsActive()) plc->StiffnessMatrix(LS, tp);
	}
}
