#include "stdafx.h"
#include "fem.h"
#include "FESolidSolver.h"
#include "FEMicroMaterial.h"
#include "FEBioLib/FETrussMaterial.h"
#include "FEBioLib/FEPointBodyForce.h"
#include "FEBioLib/FEElasticDomain.h"

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : FESolidSolver::Update
// Updates the current nodal positions based on the displacement increment
// and line search factor. If bfinal is true, it also updates the current
// displacements
//

void FESolidSolver::Update(vector<double>& ui)
{
	int i, n;

	// get the mesh
	FEM& fem = dynamic_cast<FEM&>(m_fem);
	FEMesh& mesh = m_fem.m_mesh;

	// update rigid bodies
	UpdateRigidBodies(ui);

	// update flexible nodes
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		// displacement dofs
		// current position = initial + total at prev conv step + total increment so far + current increment  
		if ((n = node.m_ID[DOF_X]) >= 0) node.m_rt.x = node.m_r0.x + m_Ut[n] + m_Ui[n] + ui[n];
		if ((n = node.m_ID[DOF_Y]) >= 0) node.m_rt.y = node.m_r0.y + m_Ut[n] + m_Ui[n] + ui[n];
		if ((n = node.m_ID[DOF_Z]) >= 0) node.m_rt.z = node.m_r0.z + m_Ut[n] + m_Ui[n] + ui[n];

		// rotational dofs
		if ((n = node.m_ID[DOF_U]) >= 0) node.m_Dt.x = node.m_D0.x + m_Ut[n] + m_Ui[n] + ui[n];
		if ((n = node.m_ID[DOF_V]) >= 0) node.m_Dt.y = node.m_D0.y + m_Ut[n] + m_Ui[n] + ui[n];
		if ((n = node.m_ID[DOF_W]) >= 0) node.m_Dt.z = node.m_D0.z + m_Ut[n] + m_Ui[n] + ui[n];
	}

	// make sure the prescribed displacements are fullfilled
	int ndis = m_fem.m_DC.size();
	for (i=0; i<ndis; ++i)
	{
		FEPrescribedBC& dc = *m_fem.m_DC[i];
		if (dc.IsActive())
		{
			int n    = dc.node;
			int lc   = dc.lc;
			int bc   = dc.bc;
			double s = dc.s;
			double r = dc.r;	// GAA

			FENode& node = mesh.Node(n);

			double g = r + s*m_fem.GetLoadCurve(lc)->Value(); // GAA

			switch (bc)
			{
			case 0:
				node.m_rt.x = node.m_r0.x + g;
				break;
			case 1:
				node.m_rt.y = node.m_r0.y + g;
				break;
			case 2:
				node.m_rt.z = node.m_r0.z + g;
				break;
			case 20:
				{
					vec3d dr = node.m_r0;
					dr.x = 0; dr.unit(); dr *= g;

					node.m_rt.y = node.m_r0.y + dr.y;
					node.m_rt.z = node.m_r0.z + dr.z;
				}
				break;
			}
		}
	}

	// enforce the linear constraints
	// TODO: do we really have to do this? Shouldn't the algorithm
	// already guarantee that the linear constraints are satisfied?
	if (fem.m_LinC.size() > 0)
	{
		int nlin = fem.m_LinC.size();
		list<FELinearConstraint>::iterator it = fem.m_LinC.begin();
		double d;
		for (int n=0; n<nlin; ++n, ++it)
		{
			FELinearConstraint& lc = *it;
			FENode& node = mesh.Node(lc.master.node);

			d = 0;
			int ns = lc.slave.size();
			list<FELinearConstraint::SlaveDOF>::iterator si = lc.slave.begin();
			for (int i=0; i<ns; ++i, ++si)
			{
				FENode& node = mesh.Node(si->node);
				switch (si->bc)
				{
				case 0: d += si->val*(node.m_rt.x - node.m_r0.x); break;
				case 1: d += si->val*(node.m_rt.y - node.m_r0.y); break;
				case 2: d += si->val*(node.m_rt.z - node.m_r0.z); break;
				}
			}

			switch (lc.master.bc)
			{
			case 0: node.m_rt.x = node.m_r0.x + d; break;
			case 1: node.m_rt.y = node.m_r0.y + d; break;
			case 2: node.m_rt.z = node.m_r0.z + d; break;
			}
		}
	}


	// update velocity and accelerations
	// for dynamic simulations
	FEAnalysis* pstep = fem.GetCurrentStep();
	if (pstep->m_nanalysis == FE_DYNAMIC)
	{
		int N = mesh.Nodes();
		double dt = pstep->m_dt;
		double a = 4.0 / dt;
		double b = a / dt;
		for (i=0; i<N; ++i)
		{
			FENode& n = mesh.Node(i);
			n.m_at = (n.m_rt - n.m_rp)*b - n.m_vp*a - n.m_ap;
			n.m_vt = n.m_vp + (n.m_ap + n.m_at)*dt*0.5;
		}
	}

	// update poroelastic data
	if (pstep->GetType() == FE_BIPHASIC) UpdatePoro(ui);

	// update solute-poroelastic data
	if (pstep->GetType() == FE_POROSOLUTE) { UpdatePoro(ui); UpdateSolute(ui); }

	// update triphasic data
	if (pstep->GetType() == FE_TRIPHASIC) { UpdatePoro(ui); UpdateTriphasic(ui); }

	// update contact
	if (fem.ContactInterfaces() > 0) UpdateContact();

	// update element stresses
	UpdateStresses();

	// update other stuff that may depend on the deformation
	int NBF = m_fem.BodyForces();
	for (i=0; i<NBF; ++i)
	{
		FEPointBodyForce* pbf = dynamic_cast<FEPointBodyForce*>(m_fem.GetBodyForce(i));
		if (pbf) pbf->Update();
	}

	// dump all states to the plot file
	// when requested
	if (pstep->m_nplot == FE_PLOT_MINOR_ITRS) fem.m_plot->Write(m_fem);
}

///////////////////////////////////////////////////////////////////////////////
//! Updates the poroelastic data

void FESolidSolver::UpdatePoro(vector<double>& ui)
{
	int i, n;

	FEMesh& mesh = m_fem.m_mesh;
	FEAnalysis* pstep = m_fem.GetCurrentStep();

	// update poro-elasticity data
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		// update nodal pressures
		n = node.m_ID[DOF_P];
		if (n >= 0) node.m_pt = 0 + m_Ut[n] + m_Ui[n] + ui[n];
	}

	// update poro-elasticity data
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		// update velocities
		node.m_vt  = (node.m_rt - node.m_rp) / pstep->m_dt;
	}

	// make sure the prescribed pressures are fullfilled
	int ndis = m_fem.m_DC.size();
	for (i=0; i<ndis; ++i)
	{
		FEPrescribedBC& dc = *m_fem.m_DC[i];
		if (dc.IsActive())
		{
			int n    = dc.node;
			int lc   = dc.lc;
			int bc   = dc.bc;
			double s = dc.s;
			double r = dc.r;	// GAA

			FENode& node = mesh.Node(n);

			if (bc == DOF_P) node.m_pt = r + s*m_fem.GetLoadCurve(lc)->Value(); // GAA
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
//! Updates the rigid body data

void FESolidSolver::UpdateRigidBodies(vector<double>& ui)
{
	FEM& fem = dynamic_cast<FEM&>(m_fem);
	FEMesh& mesh = m_fem.m_mesh;

	// update rigid bodies
	int nrb = fem.m_Obj.size();
	for (int i=0; i<nrb; ++i)
	{
		// get the rigid body
		FERigidBody& RB = dynamic_cast<FERigidBody&>(*fem.m_Obj[i]);
		if (RB.IsActive()) RB.Update(m_Ui, ui);
	}

	// update rigid joints
	int NC = fem.NonlinearConstraints();
	for (int i=0; i<NC; ++i)
	{
		FENLConstraint* plc = fem.NonlinearConstraint(i);
		plc->Update();
	}
}

///////////////////////////////////////////////////////////////////////////////
//! Updates the solute data

void FESolidSolver::UpdateSolute(vector<double>& ui)
{
	int i, n;
	
	FEMesh& mesh = m_fem.m_mesh;
	FEAnalysis* pstep = m_fem.GetCurrentStep();
	
	// update solute data
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		
		// update nodal concentration
		n = node.m_ID[DOF_C];
		if (n >= 0) node.m_ct[0] = 0 + m_Ut[n] + m_Ui[n] + ui[n];
	}
	
	// update solute data
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		
		// update velocities
		node.m_vt  = (node.m_rt - node.m_rp) / pstep->m_dt;
	}
	
	// make sure the prescribed concentrations are fullfilled
	int ndis = m_fem.m_DC.size();
	for (i=0; i<ndis; ++i)
	{
		FEPrescribedBC& dc = *m_fem.m_DC[i];
		if (dc.IsActive())
		{
			int n    = dc.node;
			int lc   = dc.lc;
			int bc   = dc.bc;
			double s = dc.s;
			double r = dc.r;	// GAA
			
			FENode& node = mesh.Node(n);
			
			if (bc == DOF_C) node.m_ct[0] = r + s*m_fem.GetLoadCurve(lc)->Value();
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
//! Updates the triphasic solute data

void FESolidSolver::UpdateTriphasic(vector<double>& ui)
{
	int i, n;
	
	FEMesh& mesh = m_fem.m_mesh;
	FEAnalysis* pstep = m_fem.GetCurrentStep();
	
	// update solute data
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		
		// update nodal concentration
		n = node.m_ID[DOF_C];
		if (n >= 0) node.m_ct[0] = 0 + m_Ut[n] + m_Ui[n] + ui[n];
		n = node.m_ID[DOF_C+1];
		if (n >= 0) node.m_ct[1] = 0 + m_Ut[n] + m_Ui[n] + ui[n];
	}
	
	// update solute data
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		
		// update velocities
		node.m_vt  = (node.m_rt - node.m_rp) / pstep->m_dt;
	}
	
	// make sure the prescribed concentrations are fullfilled
	int ndis = m_fem.m_DC.size();
	for (i=0; i<ndis; ++i)
	{
		FEPrescribedBC& dc = *m_fem.m_DC[i];
		if (dc.IsActive())
		{
			int n    = dc.node;
			int lc   = dc.lc;
			int bc   = dc.bc;
			double s = dc.s;
			double r = dc.r;	// GAA
			
			FENode& node = mesh.Node(n);
			
			if (bc == DOF_C) node.m_ct[0] = r + s*m_fem.GetLoadCurve(lc)->Value();	// GAA
			if (bc == DOF_C+1) node.m_ct[1] = r + s*m_fem.GetLoadCurve(lc)->Value();	// GAA
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
//  Updates the element stresses
//

void FESolidSolver::UpdateStresses()
{
	FEMesh& mesh = m_fem.m_mesh;

	// update the stresses on all domains
	for (int i=0; i<mesh.Domains(); ++i)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
		dom.UpdateStresses(m_fem);
	}
}

///////////////////////////////////////////////////////////////////////////////
//  Update contact data
//
void FESolidSolver::UpdateContact()
{
	// Update all contact interfaces
	for (int i=0; i<m_fem.ContactInterfaces(); ++i) m_fem.ContactInterface(i)->Update(m_niter);
}
