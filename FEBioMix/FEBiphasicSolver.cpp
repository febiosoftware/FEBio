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
#include "FEBiphasicSolver.h"
#include "FEBiphasicDomain.h"
#include "FESlidingInterface2.h"
#include "FESlidingInterface3.h"
#include "FESlidingInterfaceBiphasic.h"
#include <FEBioMech/FEElasticDomain.h>
#include <FEBioMech/FEElasticDomain.h>
#include <FEBioMech/FEPressureLoad.h>
#include <FEBioMech/FEResidualVector.h>
#include <FECore/log.h>
#include <FECore/sys.h>
#include <FECore/FEModel.h>
#include <FECore/FEModelLoad.h>
#include <FECore/FENodalLoad.h>
#include <FECore/FEAnalysis.h>

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_FECORE_CLASS(FEBiphasicSolver, FESolidSolver2)
	ADD_PARAMETER(m_Ptol, "ptol"        );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEBiphasicSolver::FEBiphasicSolver(FEModel* pfem) : FESolidSolver2(pfem)
{
	m_Ptol = 0.01;
	m_ndeq = 0;
	m_npeq = 0;

    m_msymm = REAL_UNSYMMETRIC; // assume non-symmetric stiffness matrix by default
    
	// Allocate degrees of freedom
	DOFS& dofs = pfem->GetDOFS();
	int varP = dofs.AddVariable("fluid pressure");
	dofs.SetDOFName(varP, 0, "p");
    int varQ = dofs.AddVariable("shell fluid pressure");
    dofs.SetDOFName(varQ, 0, "q");
    
	// get pressure dof
	m_dofP = pfem->GetDOFIndex("p");
    m_dofQ = pfem->GetDOFIndex("q");
}

//-----------------------------------------------------------------------------
//! Allocates and initializes the data structures.
//
bool FEBiphasicSolver::Init()
{
	// initialize base class
	if (FESolidSolver2::Init() == false) return false;

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
		gather(m_Ut, mesh, m_dofP);
        gather(m_Ut, mesh, m_dofQ);
    }

	return true;
}


//-----------------------------------------------------------------------------
//! Initialize equations
bool FEBiphasicSolver::InitEquations()
{
	// base class does most of the work
	FESolidSolver2::InitEquations();

	int i;

	// determined the nr of pressure and concentration equations
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	m_ndeq = m_npeq = 0;
	
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& n = mesh.Node(i);
		if (n.m_ID[m_dofX] != -1) m_ndeq++;
		if (n.m_ID[m_dofY] != -1) m_ndeq++;
		if (n.m_ID[m_dofZ] != -1) m_ndeq++;
        if (n.m_ID[m_dofSX] != -1) m_ndeq++;
        if (n.m_ID[m_dofSY] != -1) m_ndeq++;
        if (n.m_ID[m_dofSZ] != -1) m_ndeq++;
        if (n.m_ID[m_dofP] != -1) m_npeq++;
        if (n.m_ID[m_dofQ] != -1) m_npeq++;
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

	FESolidSolver2::PrepStep();
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
	double	normRi;		// initial residual norm
	double	normEi;		// initial energy norm
	double	normEm;		// max energy norm
	double	normDi;		// initial displacement norm

	// poro convergence norms data
	double	normPi;		// initial pressure norm
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
//! calculates the concentrated nodal forces
void FEBiphasicSolver::NodalForces(FEGlobalVector& R, const FETimeInfo& tp)
{
	// loop over nodal loads
	FEModel& fem = *GetFEModel();
	int NNL = fem.NodalLoads();
	for (int i=0; i<NNL; ++i)
	{
		FENodalDOFLoad& fc = dynamic_cast<FENodalDOFLoad&>(*fem.NodalLoad(i));
		if (fc.IsActive())
		{
			int dof = fc.GetDOF();

			FENodeSet& nset = *fc.GetNodeSet();
			int N = nset.Size();
			for (int j=0; j<N; ++j)
			{
				int nid = nset[j];	// node ID

				// get the nodal load value
				double f = fc.NodeValue(j);
			
				// For pressure and concentration loads, multiply by dt
				// for consistency with evaluation of residual and stiffness matrix
				if ((dof == m_dofP) || (dof == m_dofQ)) f *= tp.timeIncrement;

				// assemble into residual
				AssembleResidual(nid, dof, f, R);
			}
		}
	}
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

	// initialize residual with concentrated nodal loads
	R = m_Fn;

	// zero nodal reaction forces
	zero(m_Fr);

	// setup global RHS vector
	FEResidualVector RHS(fem, R, m_Fr);

	// zero rigid body reaction forces
	m_rigidSolver.Residual();

	// get the mesh
	FEMesh& mesh = fem.GetMesh();

	// calculate internal stress force
	if (fem.GetCurrentStep()->m_nanalysis == FE_STEADY_STATE)
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

    // calculate the body forces
	for (int j = 0; j<fem.BodyLoads(); ++j)
	{
		FEBodyLoad* pbl = fem.GetBodyLoad(j);
		if (pbl->IsActive()) pbl->Residual(RHS, tp);
    }
    
	// calculate forces due to surface loads
	for (int i=0; i<fem.SurfaceLoads(); ++i)
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
		if ((n = -node.m_ID[m_dofX]-2) >= 0) node.m_Fr.x = -m_Fr[n];
		if ((n = -node.m_ID[m_dofY]-2) >= 0) node.m_Fr.y = -m_Fr[n];
		if ((n = -node.m_ID[m_dofZ]-2) >= 0) node.m_Fr.z = -m_Fr[n];
	}

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

	// calculate the stiffness matrix for each domain
	FEAnalysis* pstep = fem.GetCurrentStep();
	bool bsymm = (m_msymm == REAL_SYMMETRIC);
	if (pstep->m_nanalysis == FE_STEADY_STATE)
	{
		for (int i=0; i<mesh.Domains(); ++i) 
		{
            // Biphasic analyses may include biphasic and elastic domains
			FEBiphasicDomain* pbdom = dynamic_cast<FEBiphasicDomain*>(&mesh.Domain(i));
			if (pbdom) pbdom->StiffnessMatrixSS(this, bsymm);
            else
			{
				FEElasticDomain* pedom = dynamic_cast<FEElasticDomain*>(&mesh.Domain(i));
				if (pedom) pedom->StiffnessMatrix(this);
			}
		}
	}
	else
	{
		for (int i=0; i<mesh.Domains(); ++i) 
		{
            // Biphasic analyses may include biphasic and elastic domains
			FEBiphasicDomain* pbdom = dynamic_cast<FEBiphasicDomain*>(&mesh.Domain(i));
			if (pbdom) pbdom->StiffnessMatrix(this, bsymm);
            else 
			{
				FEElasticDomain* pedom = dynamic_cast<FEElasticDomain*>(&mesh.Domain(i));
				if (pedom) pedom->StiffnessMatrix(this);
			}
		}
	}

    // calculate the body force stiffness matrix for each domain
	// TODO: This is not  going to work with FEDiscreteSpringDomain
	int NBL = fem.BodyLoads();
	for (int j = 0; j<NBL; ++j)
	{
		FEBodyLoad* pbl = fem.GetBodyLoad(j);
		if (pbl->IsActive()) pbl->StiffnessMatrix(this, tp);
    }
    
	// calculate contact stiffness
	ContactStiffness();

	// calculate stiffness matrices for surface loads
	int nsl = fem.SurfaceLoads();
	for (int i=0; i<nsl; ++i)
	{
		FESurfaceLoad* psl = fem.SurfaceLoad(i);
		if (psl->IsActive()) psl->StiffnessMatrix(this, tp);
	}

	// calculate nonlinear constraint stiffness
	// note that this is the contribution of the 
	// constrainst enforced with augmented lagrangian
	NonLinearConstraintStiffness(tp);

	// add contributions from rigid bodies
	m_rigidSolver.StiffnessMatrix(*m_pK, tp);

	return true;
}

//-----------------------------------------------------------------------------
//! Update the model's kinematic data. This is overriden from FESolidSolver2 so
//! that biphasic data is updated
void FEBiphasicSolver::UpdateKinematics(vector<double>& ui)
{
	// first update all solid-mechanics kinematics
	FESolidSolver2::UpdateKinematics(ui);

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
		n = node.m_ID[m_dofP];
		if (n >= 0) node.set(m_dofP, 0 + m_Ut[n] + m_Ui[n] + ui[n]);
        n = node.m_ID[m_dofQ];
        if (n >= 0) node.set(m_dofQ, 0 + m_Ut[n] + m_Ui[n] + ui[n]);
    }

	// update poro-elasticity data
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		// update velocities
		vec3d vt = (node.m_rt - node.m_rp) / dt;
		node.set_vec3d(m_dofVX, m_dofVY, m_dofVZ, vt); 
	}
}

//-----------------------------------------------------------------------------
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
	}

	// Update all contact interfaces
	FESolidSolver2::UpdateModel();

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
		nid = n.m_ID[m_dofX];
		if (nid != -1)
		{
			nid = (nid < -1 ? -nid-2 : nid);
			di[m++] = ui[nid];
			assert(m <= (int) di.size());
		}
		nid = n.m_ID[m_dofY];
		if (nid != -1)
		{
			nid = (nid < -1 ? -nid-2 : nid);
			di[m++] = ui[nid];
			assert(m <= (int) di.size());
		}
		nid = n.m_ID[m_dofZ];
		if (nid != -1)
		{
			nid = (nid < -1 ? -nid-2 : nid);
			di[m++] = ui[nid];
			assert(m <= (int) di.size());
		}
        nid = n.m_ID[m_dofSX];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            di[m++] = ui[nid];
            assert(m <= (int) di.size());
        }
        nid = n.m_ID[m_dofSY];
        if (nid != -1)
        {
            nid = (nid < -1 ? -nid-2 : nid);
            di[m++] = ui[nid];
            assert(m <= (int) di.size());
        }
        nid = n.m_ID[m_dofSZ];
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
		nid = n.m_ID[m_dofP];
		if (nid != -1)
		{
			nid = (nid < -1 ? -nid-2 : nid);
			pi[m++] = ui[nid];
			assert(m <= (int) pi.size());
		}
        nid = n.m_ID[m_dofQ];
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
	FESolidSolver2::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_Ptol & m_ndeq & m_npeq;
	ar & m_nceq;
	ar & m_di & m_Di & m_pi & m_Pi;
}
