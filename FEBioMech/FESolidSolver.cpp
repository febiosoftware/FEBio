/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FESolidSolver.h"
#include "FE3FieldElasticSolidDomain.h"
#include "FEResidualVector.h"
#include "FEUncoupledMaterial.h"
#include "FEContactInterface.h"
#include <FECore/sys.h>
#include <FECore/log.h>
#include <FECore/DOFS.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FEBoundaryCondition.h>
#include <FECore/FENodalLoad.h>
#include <FECore/FEModelLoad.h>
#include <FECore/FELinearConstraintManager.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEBodyLoad.h>
#include <assert.h>
#include "FEBioMech.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_FECORE_CLASS(FESolidSolver, FENewtonSolver)
	ADD_PARAMETER(m_Dtol        , FE_RANGE_GREATER_OR_EQUAL(0.0), "dtol"        );
	ADD_PARAMETER(m_Etol        , FE_RANGE_GREATER_OR_EQUAL(0.0), "etol"        );
	ADD_PARAMETER(m_Rtol        , FE_RANGE_GREATER_OR_EQUAL(0.0), "rtol"        );
	ADD_PARAMETER(m_Rmin        , FE_RANGE_GREATER_OR_EQUAL(0.0), "min_residual");
	ADD_PARAMETER(m_beta        , "beta"        );
	ADD_PARAMETER(m_gamma       , "gamma"       );
	ADD_PARAMETER(m_bnew_update , "use_new_rigid_update");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! FESolidSolver Construction
//
FESolidSolver::FESolidSolver(FEModel* pfem) : FENewtonSolver(pfem), m_rigidSolver(pfem),\
m_dofU(pfem), m_dofV(pfem), m_dofSQ(pfem), m_dofRQ(pfem)
{
	// default values
	m_Rtol = 0;	// deactivate residual convergence 
	m_Dtol = 0.001;
	m_Etol = 0.01;
	m_Rmin = 1.0e-20;

	m_niter = 0;
	m_nreq = 0;

	// default Newmark parameters for unconditionally stable time integration
	m_beta = 0.25;
	m_gamma = 0.5;

	m_bnew_update = false;

	m_rigidSolver.AllowMixedBCs(true);

	// Allocate degrees of freedom
	DOFS& dofs = pfem->GetDOFS();
	int varU = dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT), VAR_VEC3);
	dofs.SetDOFName(varU, 0, "x");
	dofs.SetDOFName(varU, 1, "y");
	dofs.SetDOFName(varU, 2, "z");
	int varSQ = dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_ROTATION), VAR_VEC3);
	dofs.SetDOFName(varSQ, 0, "u");
	dofs.SetDOFName(varSQ, 1, "v");
	dofs.SetDOFName(varSQ, 2, "w");
	int varRQ = dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::RIGID_ROTATION), VAR_VEC3);
	dofs.SetDOFName(varRQ, 0, "Ru");
	dofs.SetDOFName(varRQ, 1, "Rv");
	dofs.SetDOFName(varRQ, 2, "Rw");
	int varV = dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::VELOCTIY), VAR_VEC3);
	dofs.SetDOFName(varV, 0, "vx");
	dofs.SetDOFName(varV, 1, "vy");
	dofs.SetDOFName(varV, 2, "vz");

	// get the DOF indices
	m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
	m_dofV.AddVariable(FEBioMech::GetVariableName(FEBioMech::VELOCTIY));
	m_dofSQ.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_ROTATION));
	m_dofRQ.AddVariable(FEBioMech::GetVariableName(FEBioMech::RIGID_ROTATION));
}

//-----------------------------------------------------------------------------
FESolidSolver::~FESolidSolver()
{
}

//-----------------------------------------------------------------------------
//! Allocates and initializes the data structures used by the FESolidSolver
//
bool FESolidSolver::Init()
{
	// initialize base class
	if (FENewtonSolver::Init() == false) return false;

	// allocate vectors
	int neq = m_neq;
	m_Fn.assign(neq, 0);
	m_Fr.assign(neq, 0);
	m_Ui.assign(neq, 0);
	m_Ut.assign(neq, 0);

	// we need to fill the total displacement vector m_Ut
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	gather(m_Ut, mesh, m_dofU[0]);
	gather(m_Ut, mesh, m_dofU[1]);
	gather(m_Ut, mesh, m_dofU[2]);
	gather(m_Ut, mesh, m_dofSQ[0]);
	gather(m_Ut, mesh, m_dofSQ[1]);
	gather(m_Ut, mesh, m_dofSQ[2]);

	// set the dynamic update flag only if we are running a dynamic analysis
	bool b = (fem.GetCurrentStep()->m_nanalysis == FE_DYNAMIC ? true : false);
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEElasticSolidDomain* d = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(i));
		if (d) d->SetDynamicUpdateFlag(b);
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Save data to dump file

void FESolidSolver::Serialize(DumpStream& ar)
{
	// Serialize parameters
	FENewtonSolver::Serialize(ar);
	
	ar & m_nrhs;
	ar & m_niter;
	ar & m_nref & m_ntotref;
	ar & m_naug;
	ar & m_nreq;

	if (ar.IsLoading())
	{
		// re-initialize data
		Init();
	}

	m_rigidSolver.Serialize(ar);
}

//-----------------------------------------------------------------------------
//! Determine the number of linear equations and assign equation numbers
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
bool FESolidSolver::InitEquations()
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
bool FESolidSolver::Augment()
{
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();

	// Assume we will pass (can't hurt to be optimistic)
	bool bconv = true;

	// Do contact augmentations
	if (fem.SurfacePairConstraints() > 0)
	{
		// loop over all contact interfaces
		for (int i=0; i<fem.SurfacePairConstraints(); ++i)
		{
			FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));
			if (pci->IsActive()) bconv = (pci->Augment(m_naug, tp) && bconv);
		}
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
//! Update the kinematics of the model, such as nodal positions, velocities,
//! accelerations, etc.
void FESolidSolver::UpdateKinematics(vector<double>& ui)
{
	// get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// update rigid bodies
	// (this also updates the kinematics of rigid nodes)
	m_rigidSolver.UpdateRigidBodies(m_Ui, ui, m_bnew_update);

	// total displacements
	vector<double> U(m_Ut.size());
	for (size_t i=0; i<m_Ut.size(); ++i) U[i] = ui[i] + m_Ui[i] + m_Ut[i];

	// update flexible nodes
	// translational dofs
	scatter(U, mesh, m_dofU[0]);
	scatter(U, mesh, m_dofU[1]);
	scatter(U, mesh, m_dofU[2]);
	// rotational dofs
	scatter(U, mesh, m_dofSQ[0]);
	scatter(U, mesh, m_dofSQ[1]);
	scatter(U, mesh, m_dofSQ[2]);

	// make sure the boundary conditions are fullfilled
	int nbcs = fem.BoundaryConditions();
	for (int i = 0; i<nbcs; ++i)
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
	for (int i = 0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		if (node.m_rid == -1)
			node.m_rt = node.m_r0 + node.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2]);
	}

	// update velocity and accelerations
	// for dynamic simulations
	FEAnalysis* pstep = fem.GetCurrentStep();
	if (pstep->m_nanalysis == FE_DYNAMIC)
	{
		int N = mesh.Nodes();
		double dt = fem.GetTime().timeIncrement;
		double a = 1.0 / (m_beta*dt);
		double b = a / dt;
		double c = 1.0 - 0.5/m_beta;
		for (int i = 0; i<N; ++i)
		{
			FENode& n = mesh.Node(i);
			n.m_at = (n.m_rt - n.m_rp)*b - n.m_vp*a + n.m_ap*c;
			vec3d vt = n.m_vp + (n.m_ap*(1.0 - m_gamma) + n.m_at*m_gamma)*dt;
			n.set_vec3d(m_dofV[0], m_dofV[1], m_dofV[2], vt);
		}

		// update the rigid body kinematics
		m_rigidSolver.UpdateRigidKinematics();
	}
}

//-----------------------------------------------------------------------------
//! Updates the current state of the model
void FESolidSolver::Update(vector<double>& ui)
{
	FEModel& fem = *GetFEModel();

	// update kinematics
	UpdateKinematics(ui);

	// update element stresses
	fem.Update();
}

//-----------------------------------------------------------------------------
//! Prepares the data for the first BFGS-iteration. 
void FESolidSolver::PrepStep()
{
	// zero total displacements
	zero(m_Ui);

	// store previous mesh state
	// we need them for velocity and acceleration calculations
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& ni = mesh.Node(i);
		ni.m_rp = ni.m_rt;
		ni.m_vp = ni.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);
		ni.m_ap = ni.m_at;
	}

	const FETimeInfo& tp = fem.GetTime();

	// apply concentrated nodal forces
	// since these forces do not depend on the geometry
	// we can do this once outside the NR loop.
	vector<double> dummy(m_neq, 0.0);
	zero(m_Fn);
	FEResidualVector Fn(*GetFEModel(), m_Fn, dummy);
	NodalForces(Fn, tp);

	// apply boundary conditions
	// we save the prescribed displacements increments in the ui vector
	vector<double>& ui = m_ui;
	zero(ui);
	int nbc = fem.BoundaryConditions();
	for (int i = 0; i<nbc; ++i)
	{
		FEBoundaryCondition& dc = *fem.BoundaryCondition(i);
		if (dc.IsActive()) dc.PrepStep(ui);
	}

	// initialize rigid bodies
	m_rigidSolver.PrepStep(tp, ui);

	// intialize material point data
	// NOTE: do this before the stresses are updated
	// TODO: does it matter if the stresses are updated before
	//       the material point data is initialized
	for (int i=0; i<mesh.Domains(); ++i) mesh.Domain(i).PreSolveUpdate(tp);

	// update stresses
	fem.Update();

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
	for (int i=0; i<fem.NonlinearConstraints(); ++i)
	{
		FENLConstraint& ci = *fem.NonlinearConstraint(i);
		if (ci.IsActive()) m_baugment = true;
	}
}

//-----------------------------------------------------------------------------
//! Implements the BFGS algorithm to solve the nonlinear FE equations.
bool FESolidSolver::Quasin()
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
	FEModel& fem = *GetFEModel();
	FEAnalysis* pstep = fem.GetCurrentStep();
	const FETimeInfo& tp = fem.GetTime();

	// prepare for the first iteration
	PrepStep();

	// Init QN method
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

		// set initial convergence norms
		if (m_niter == 0)
		{
			normRi = fabs(m_R0*m_R0);
			normEi = fabs(m_ui*m_R0);
			normUi = fabs(m_ui*m_ui);
			normEm = normEi;
		}

		// calculate norms
		normR1 = m_R1*m_R1;
		normu  = (m_ui*m_ui)*(s*s);
		normE1 = s*fabs(m_ui*m_R1);

		// check for nans
		if (ISNAN(normR1) || ISNAN(normu)) throw NANDetected();

		// update total displacements
		int neq = (int)m_Ui.size();
		for (int i = 0; i<neq; ++i) m_Ui[i] += s*m_ui[i];
		normU  = m_Ui*m_Ui;

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
		feLog(" Nonlinear solution status: time= %lg\n", tp.currentTime); 
		feLog("\tstiffness updates             = %d\n", m_qnstrategy->m_nups);
		feLog("\tright hand side evaluations   = %d\n", m_nrhs);
		feLog("\tstiffness matrix reformations = %d\n", m_nref);
		if (m_lineSearch->m_LStol > 0) feLog("\tstep from line search         = %lf\n", s);
		feLog("\tconvergence norms :     INITIAL         CURRENT         REQUIRED\n");
		feLog("\t   residual         %15le %15le %15le \n", normRi, normR1, m_Rtol*normRi);
		feLog("\t   energy           %15le %15le %15le \n", normEi, normE1, m_Etol*normEi);
		feLog("\t   displacement     %15le %15le %15le \n", normUi, normu ,(m_Dtol*m_Dtol)*normU );

		// see if we may have a small residual
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
//! Calculates global stiffness matrix.
bool FESolidSolver::StiffnessMatrix()
{
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();

	// get the mesh
	FEMesh& mesh = fem.GetMesh();

	// calculate the stiffness matrix for each domain
	for (int i=0; i<mesh.Domains(); ++i) 
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
		dom.StiffnessMatrix(this);
	}

	// calculate the body force stiffness matrix for each domain
	int NBL = fem.BodyLoads();
	for (int j = 0; j<NBL; ++j)
	{
		FEBodyLoad* pbl = fem.GetBodyLoad(j);
		if (pbl->IsActive()) pbl->StiffnessMatrix(this, tp);
	}

	// Add mass matrix for dynamic problems
	FEAnalysis* pstep = fem.GetCurrentStep();
	if (pstep->m_nanalysis == FE_DYNAMIC)
	{
		// scale factor
		double dt = tp.timeIncrement;
		double a = 1.0 / (m_beta*dt*dt);

		// loop over all domains
		for (int i = 0; i<mesh.Domains(); ++i)
		{
			FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
			dom.MassMatrix(this, a);
		}
	}

	// calculate contact stiffness
	if (fem.SurfacePairConstraints() > 0)
	{
		ContactStiffness();
	}

	// calculate stiffness matrices for surface loads
	int nsl = fem.SurfaceLoads();
	for (int i = 0; i<nsl; ++i)
	{
		FESurfaceLoad* psl = fem.SurfaceLoad(i);
		if (psl->IsActive())
		{
			psl->StiffnessMatrix(this, tp); 
		}
	}

	// calculate nonlinear constraint stiffness
	// note that this is the contribution of the 
	// constrainst enforced with augmented lagrangian
	NonLinearConstraintStiffness(tp);

	// calculate the stiffness contributions for the rigid forces
	for (int i = 0; i<fem.ModelLoads(); ++i) fem.ModelLoad(i)->StiffnessMatrix(this, tp);

	// we still need to set the diagonal elements to 1
	// for the prescribed rigid body dofs.
	m_rigidSolver.StiffnessMatrix(*m_pK, tp);

	return true;
}

//-----------------------------------------------------------------------------
//! Calculate the stiffness contribution due to nonlinear constraints
void FESolidSolver::NonLinearConstraintStiffness(const FETimeInfo& tp)
{
	FEModel& fem = *GetFEModel();
	int N = fem.NonlinearConstraints();
	for (int i=0; i<N; ++i) 
	{
		FENLConstraint* plc = fem.NonlinearConstraint(i);
		if (plc->IsActive()) plc->StiffnessMatrix(this, tp);
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the contact stiffness matrix

void FESolidSolver::ContactStiffness()
{
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();
	for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));
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
void FESolidSolver::AssembleResidual(int node_id, int dof, double f, vector<double>& R)
{
	// get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

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
void FESolidSolver::AssembleStiffness(std::vector<int>& lm, matrix& ke)
{
	m_pK->Assemble(ke, lm);

	// adjust for prescribed dofs
	vector<double>& ui = m_ui;

	SparseMatrix& K = *m_pK;

	// loop over columns
	int cols = ke.columns();
	int rows = ke.rows();
	for (int j = 0; j<cols; ++j)
	{
		int J = -lm[j] - 2;
		if ((J >= 0) && (J<m_neq))
		{
			// dof j is a prescribed degree of freedom

			// loop over rows
			for (int i = 0; i<rows; ++i)
			{
				int I = lm[i];
				if (I >= 0)
				{
					// dof i is not a prescribed degree of freedom
					m_Fd[I] -= ke[i][j] * ui[J];
				}
			}

			// set the diagonal element of K to 1
			K.set(J, J, 1);
		}
	}
}

//-----------------------------------------------------------------------------
//!  Assembles the element stiffness matrix into the global stiffness matrix.
//!  Also adjusts the global stiffness matrix and residual to take the 
//!  prescribed displacements into account.

//! \todo In stead of changing the global stiffness matrix to accomodate for 
//!       the rigid bodies and linear constraints, can I modify the element stiffness
//!       matrix prior to assembly? I might have to change the elm vector as well as 
//!       the element matrix size.

void FESolidSolver::AssembleStiffness(vector<int>& en, vector<int>& elmi, vector<int>& elmj, matrix& ke)
{
	// assemble into global stiffness matrix
	m_pK->Assemble(ke, elmi, elmj);

	vector<double>& ui = m_ui;

	// adjust for linear constraints
	FEModel& fem = *GetFEModel();
	FELinearConstraintManager& LCM = fem.GetLinearConstraintManager();
	if (LCM.LinearConstraints() > 0)
	{
		LCM.AssembleStiffness(*m_pK, m_Fd, m_ui, en, elmi, elmj, ke);
	}

	// adjust stiffness matrix for prescribed degrees of freedom
	// NOTE: I had to comment this if statement out since otherwise
	//       poroelastic DOF's that are set as free-draining in the
	//       sliding2 contact code are skipt and zeroes will appear
	//       on the diagonal of the stiffness matrix.
//	if (fem.m_DC.size() > 0)
	{
		int i, j;
		int I, J;

		SparseMatrix& K = *m_pK;

		int N = ke.rows();

		// loop over columns
		for (j=0; j<N; ++j)
		{
			J = -elmj[j]-2;
			if ((J >= 0) && (J<m_nreq))
			{
				// dof j is a prescribed degree of freedom

				// loop over rows
				for (i=0; i<N; ++i)
				{
					I = elmi[i];
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
	m_rigidSolver.RigidStiffness(*m_pK, m_ui, m_Fd, en, elmi, elmj, ke, 1.0);
}

//-----------------------------------------------------------------------------
//! Calculates the contact forces
void FESolidSolver::ContactForces(FEGlobalVector& R)
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

bool FESolidSolver::Residual(vector<double>& R)
{
	// get the time information
	FEModel& fem = *GetFEModel();
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

	// calculate the internal (stress) forces
	for (int i=0; i<mesh.Domains(); ++i)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
		dom.InternalForces(RHS);
	}

	// calculate the body forces
	for (int j = 0; j<fem.BodyLoads(); ++j)
	{
		FEBodyLoad* pbl = fem.GetBodyLoad(j);
		if (pbl->IsActive()) pbl->Residual(RHS, tp);
	}

	// calculate inertial forces for dynamic problems
	if (fem.GetCurrentStep()->m_nanalysis == FE_DYNAMIC) InertialForces(RHS);

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

	// forces due to point constraints
//	for (i=0; i<(int) fem.m_PC.size(); ++i) fem.m_PC[i]->Residual(this, R);

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
void FESolidSolver::NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp)
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

void FESolidSolver::NodalForces(FEGlobalVector& R, const FETimeInfo& tp)
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

//-----------------------------------------------------------------------------
//! This function calculates the inertial forces for dynamic problems

void FESolidSolver::InertialForces(FEGlobalVector& R)
{
	// get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// allocate F
	vector<double> F(3*mesh.Nodes());
	zero(F);

	// get the time information
	const FETimeInfo& tp = fem.GetTime();

	// calculate F
	double dt = tp.timeIncrement;
	double a = 1.0 / (m_beta*dt);
	double b = a / dt;
	double c = 1.0 - 0.5/m_beta;
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		vec3d& rt = node.m_rt;
		vec3d& rp = node.m_rp;
		vec3d& vp = node.m_vp;
		vec3d& ap = node.m_ap;

		F[3*i  ] = b*(rt.x - rp.x) - a*vp.x + c * ap.x;
		F[3*i+1] = b*(rt.y - rp.y) - a*vp.y + c * ap.y;
		F[3*i+2] = b*(rt.z - rp.z) - a*vp.z + c * ap.z;
	}

	// now multiply F with the mass matrix
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(nd));
		dom.InertialForces(R, F);
	}
}
