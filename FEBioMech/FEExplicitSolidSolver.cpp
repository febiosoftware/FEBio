/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "FEExplicitSolidSolver.h"
#include "FEElasticSolidDomain.h"
#include "FERigidMaterial.h"
#include "FEBodyForce.h"
#include "FEContactInterface.h"
#include "FERigidBody.h"
#include "FERigidSystem.h"
#include "RigidBC.h"
#include "FEMechModel.h"
#include <FECore/FENodalLoad.h>
#include <FECore/FESurfaceLoad.h>
#include "FECore/log.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FELinearConstraintManager.h>
#include "FEResidualVector.h"
#include "FEBioMech.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_FECORE_CLASS(FEExplicitSolidSolver, FESolver)
	ADD_PARAMETER(m_dyn_damping, "dyn_damping");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEExplicitSolidSolver::FEExplicitSolidSolver(FEModel* pfem) : FESolver(pfem), m_dofU(pfem), m_dofV(pfem), m_dofSQ(pfem), m_dofRQ(pfem)
{
	m_dyn_damping = 0.99;
	m_niter = 0;
	m_nreq = 0;

	// Allocate degrees of freedom
	DOFS& dofs = pfem->GetDOFS();
	int varD = dofs.AddVariable("displacement", VAR_VEC3);
	dofs.SetDOFName(varD, 0, "x");
	dofs.SetDOFName(varD, 1, "y");
	dofs.SetDOFName(varD, 2, "z");
	int varQ = dofs.AddVariable("shell rotation", VAR_VEC3);
	dofs.SetDOFName(varQ, 0, "u");
	dofs.SetDOFName(varQ, 1, "v");
	dofs.SetDOFName(varQ, 2, "w");
	int varQR = dofs.AddVariable("rigid rotation", VAR_VEC3);
	dofs.SetDOFName(varQR, 0, "Ru");
	dofs.SetDOFName(varQR, 1, "Rv");
	dofs.SetDOFName(varQR, 2, "Rw");
	int varV = dofs.AddVariable("velocity", VAR_VEC3);
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
void FEExplicitSolidSolver::Clean()
{
}

//-----------------------------------------------------------------------------
bool FEExplicitSolidSolver::Init()
{
	if (FESolver::Init() == false) return false;

	// get nr of equations
	int neq = m_neq;

	// allocate vectors
	m_Fn.assign(neq, 0);
	m_Fd.assign(neq, 0);
	m_Fr.assign(neq, 0);
	m_Ui.assign(neq, 0);
	m_ui.assign(neq, 0);
	m_Ut.assign(neq, 0);
	m_inv_mass.assign(neq, 1);

	// we need to fill the total displacement vector m_Ut
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	gather(m_Ut, mesh, m_dofU[0]);
	gather(m_Ut, mesh, m_dofU[1]);
	gather(m_Ut, mesh, m_dofU[2]);
	gather(m_Ut, mesh, m_dofSQ[0]);
	gather(m_Ut, mesh, m_dofSQ[1]);
	gather(m_Ut, mesh, m_dofSQ[2]);

	// calculate the inverse mass vector for the explicit analysis
	vector<double> dummy(m_inv_mass);
	FEGlobalVector Mi(fem, m_inv_mass, dummy);
	matrix ke;
	int j, iel;
	int nint, neln;
	double *H, kab;
	vector <int> lm;
	vector <double> el_lumped_mass;
	// Data structure to store element mass data for dynamic damping:-
	// Define an overall dynamic array of pointers to the array that points to the element data
	// domain_mass is a list of pointers to the data for each domain
	domain_mass = new double ** [mesh.Domains()];

	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		// check whether it is a solid domain
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(nd));
		if (pbd)  // it is an elastic solid domain
		{
			// for each domain define an array of pointers to the individual element_mass records
			double ** elmasses; 
			elmasses = new double * [pbd->Elements()];
			// and set a pointer in domain_mass to the new element array
			domain_mass[nd] = elmasses;

			FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(pbd->GetMaterial());
			FEParamDouble& density = pme->Density();

			for (iel=0; iel<pbd->Elements(); ++iel)
			{
				FESolidElement& el = pbd->Element(iel);
				pbd->UnpackLM(el, lm);

				nint = el.GaussPoints();
				neln = el.Nodes();

				ke.resize(3*neln, 3*neln);
				ke.zero();
				el_lumped_mass.resize(3*neln);
				for (int i=0; i<3*neln; ++i) el_lumped_mass[i]=0.0;

				// create the element mass matrix
				for (int n=0; n<nint; ++n)
				{
					FEMaterialPoint& mp = *el.GetMaterialPoint(n);
					double d = density(mp);
					double detJ0 = pbd->detJ0(el, n)*el.GaussWeights()[n];

					H = el.H(n);
					for (int i=0; i<neln; ++i)
						for (j=0; j<neln; ++j)
						{
							kab = H[i]*H[j]*detJ0*d;
							ke[3*i  ][3*j  ] += kab;
							ke[3*i+1][3*j+1] += kab;
							ke[3*i+2][3*j+2] += kab;
						}	
				}
				// reduce to a lumped mass vector and add up the total
				double total_mass = 0.0;
				for (int i=0; i<3*neln; ++i)
				{
						for (j=0; j<neln; ++j)
						{
							el_lumped_mass[i]+=ke[i][j];
						}
						total_mass += el_lumped_mass[i];
				}	
				total_mass /= 3.0; // because each mass is represented three times for each direction
				// define an element mass record
				double * thiselement;
				thiselement = new double [neln+1];
				// and set a pointer to the element data
				elmasses[iel] = thiselement;
				thiselement[0] = total_mass; // total mass of the element first, followed by the fraction at each node
				// for each node, store the fraction of the element mass associated with it
				for (int i=0; i<neln; ++i)
				{
					thiselement[i+1] = (el_lumped_mass[3*i]+el_lumped_mass[3*i+1]+el_lumped_mass[3*i+2])/(3*total_mass);
				} // loop over nodes within element
				// invert and assemble element matrix into inv_mass vector 
				for (int i=0; i<3*neln; ++i)
				{
					el_lumped_mass[i] = 1.0/el_lumped_mass[i];
				}
				Mi.Assemble(el.m_node, lm, el_lumped_mass);
				// formerly AssembleResidual(el.m_node, lm, el_lumped_mass, m_inv_mass);
			} // loop over elements
		} // was an elastic solid domain
		else domain_mass[nd] = 0;  // no masses stored for other types of domain
	}

	// Calculate initial residual to be used on the first time step
	if (Residual(m_R1) == false) return false;
	m_R1 += m_Fd;

	return true;
}

//-----------------------------------------------------------------------------
//! Updates the current state of the model
void FEExplicitSolidSolver::Update(vector<double>& ui)
{
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();

	// update kinematics
	UpdateKinematics(ui);

	// update element stresses
	fem.Update();
}

//-----------------------------------------------------------------------------
//! Update the kinematics of the model, such as nodal positions, velocities,
//! accelerations, etc.
void FEExplicitSolidSolver::UpdateKinematics(vector<double>& ui)
{
	// get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// update rigid bodies
	UpdateRigidBodies(ui);

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

	// make sure the prescribed displacements are fullfilled
	int ndis = fem.BoundaryConditions();
	for (int i=0; i<ndis; ++i)
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
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		if (node.m_rid == -1)
			node.m_rt = node.m_r0 + node.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2]);
	}
}

//-----------------------------------------------------------------------------
//! Updates the rigid body data
void FEExplicitSolidSolver::UpdateRigidBodies(vector<double>& ui)
{
	// get the number of rigid bodies
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidSystem& rigid = *fem.GetRigidSystem();
	const int NRB = rigid.Objects();

	// first calculate the rigid body displacement increments
	for (int i=0; i<NRB; ++i)
	{
		// get the rigid body
		FERigidBody& RB = *rigid.Object(i);
		int *lm = RB.m_LM;
		double* du = RB.m_du;

		if (RB.m_prb == 0)
		{
			for (int j=0; j<6; ++j)
			{
				du[j] = (lm[j] >=0 ? m_Ui[lm[j]] + ui[lm[j]] : 0);
			}
		}
	}

	// for prescribed displacements, the displacement increments are evaluated differently
	// TODO: Is this really necessary? Why can't the ui vector contain the correct values?
	const int NRD = rigid.PrescribedBCs();
	for (int i=0; i<NRD; ++i)
	{
		FERigidBodyDisplacement& dc = *rigid.PrescribedBC(i);
		if (dc.IsActive())
		{
			FERigidBody& RB = *rigid.Object(dc.GetID());
			int I = dc.GetBC();
			RB.m_du[I] = dc.Value() - RB.m_Up[I];
		}
	}

	// update the rigid bodies
	for (int i=0; i<NRB; ++i)
	{
		// get the rigid body
		FERigidBody& RB = *rigid.Object(i);
		double* du = RB.m_du;

		// This is the "old" update algorithm which has some issues. It does not produce the correct
		// rigid body orientation when the rotational degrees of freedom are prescribed.
		RB.m_rt.x = RB.m_rp.x + du[0];
		RB.m_rt.y = RB.m_rp.y + du[1];
		RB.m_rt.z = RB.m_rp.z + du[2];

		vec3d r = vec3d(du[3], du[4], du[5]);
		double w = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
		quatd dq = quatd(w, r);

		quatd Q = dq*RB.m_qp;
		Q.MakeUnit();
		RB.SetRotation(Q);

		if (RB.m_prb) du = RB.m_dul;
		RB.m_Ut[0] = RB.m_Up[0] + du[0];
		RB.m_Ut[1] = RB.m_Up[1] + du[1];
		RB.m_Ut[2] = RB.m_Up[2] + du[2];
		RB.m_Ut[3] = RB.m_Up[3] + du[3];
		RB.m_Ut[4] = RB.m_Up[4] + du[4];
		RB.m_Ut[5] = RB.m_Up[5] + du[5];
	}

	// we need to update the position of rigid nodes
	rigid.UpdateMesh();

	// Since the rigid nodes are repositioned we need to update the displacement DOFS
	FEMesh& mesh = fem.GetMesh();
	int N = mesh.Nodes();
	for (int i=0; i<N; ++i)
	{
		FENode& node = mesh.Node(i);
		if (node.m_rid >= 0)
		{
			vec3d ut = node.m_rt - node.m_r0;
			node.set_vec3d(m_dofU[0], m_dofU[1], m_dofU[2], ut);
		}
	}
}

//-----------------------------------------------------------------------------
//! Save data to dump file

void FEExplicitSolidSolver::Serialize(DumpStream& ar)
{
	FESolver::Serialize(ar);
	ar & m_nrhs & m_niter & m_nref & m_ntotref & m_naug & m_neq & m_nreq;
}

//-----------------------------------------------------------------------------
//!  This function mainly calls the DoSolve routine 
//!  and deals with exceptions that require the immediate termination of
//!	 the solution eg negative Jacobians.
bool FEExplicitSolidSolver::SolveStep()
{
	bool bret;

	try
	{
		// let's try to solve the step
		bret = DoSolve();
	}
	catch (NegativeJacobian e)
	{
		// A negative jacobian was detected
		feLogError("Negative jacobian was detected at element %d at gauss point %d\njacobian = %lg\n", e.m_iel, e.m_ng+1, e.m_vol);
		return false;
	}
	catch (MaxStiffnessReformations) // shouldn't happen for an explicit analysis!
	{
		// max nr of reformations is reached
		feLogError("Max nr of reformations reached.");
		return false;
	}
	catch (ForceConversion)
	{
		// user forced conversion of problem
		feLogWarning("User forced conversion.\nSolution might not be stable.");
		return true;
	}
	catch (IterationFailure)
	{
		// user caused a forced iteration failure
		feLogWarning("User forced iteration failure.");
		return false;
	}
	catch (ZeroLinestepSize) // shouldn't happen for an explicit analysis!
	{
		// a zero line step size was detected
		feLogError("Zero line step size.");
		return false;
	}
	catch (EnergyDiverging) // shouldn't happen for an explicit analysis!
	{
		// problem was diverging after stiffness reformation
		feLogError("Problem diverging uncontrollably.");
		return false;
	}
	catch (FEMultiScaleException)
	{
		// the RVE problem didn't solve
		feLogError("The RVE problem has failed. Aborting macro run.");
		return false;
	}

	return bret;
}


//-----------------------------------------------------------------------------
//! Prepares the data for the time step. 
void FEExplicitSolidSolver::PrepStep()
{
	int i, j;

	// initialize counters
	m_niter = 0;	// nr of iterations
	m_nrhs  = 0;	// nr of RHS evaluations
	m_nref  = 0;	// nr of stiffness reformations
	m_ntotref = 0;
	m_naug  = 0;	// nr of augmentations

	// zero total displacements
	zero(m_Ui);

	// store previous mesh state
	// we need them for velocity and acceleration calculations
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FEMesh& mesh = fem.GetMesh();
	for (i=0; i<mesh.Nodes(); ++i)
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

	// apply prescribed displacements
	// we save the prescribed displacements increments in the ui vector
	vector<double>& ui = m_ui;
	zero(ui);
	int neq = m_neq;
	int nbc = fem.BoundaryConditions();
	for (i=0; i<nbc; ++i)
	{
		FEBoundaryCondition& bc = *fem.BoundaryCondition(i);
		if (bc.IsActive()) bc.PrepStep(ui);
	}

	// initialize rigid bodies
	FERigidSystem& rigid = *fem.GetRigidSystem();
	int NO = rigid.Objects();
	for (i=0; i<NO; ++i) rigid.Object(i)->Init();

	// calculate local rigid displacements
	for (i=0; i<(int) rigid.PrescribedBCs(); ++i)
	{
		FERigidBodyDisplacement& DC = *rigid.PrescribedBC(i);
		FERigidBody& RB = *rigid.Object(DC.GetID());
		if (DC.IsActive())
		{
			int I = DC.GetBC();
			RB.m_dul[I] = DC.Value() - RB.m_Ut[I];
		}
	}

	// calculate global rigid displacements
	for (i=0; i<NO; ++i)
	{
		FERigidBody* prb = rigid.Object(i);
		if (prb)
		{
			FERigidBody& RB = *prb;
			if (RB.m_prb == 0)
			{
				for (j=0; j<6; ++j) RB.m_du[j] = RB.m_dul[j];
			}
			else
			{
				double* dul = RB.m_dul;
				vec3d dr = vec3d(dul[0], dul[1], dul[2]);
				
				vec3d v = vec3d(dul[3], dul[4], dul[5]);
				double w = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
				quatd dq = quatd(w, v);

				FERigidBody* pprb = RB.m_prb;

				vec3d r0 = RB.m_rt;
				quatd Q0 = RB.GetRotation();

				dr = Q0*dr;
				dq = Q0*dq*Q0.Inverse();

				while (pprb)
				{
					vec3d r1 = pprb->m_rt;
					dul = pprb->m_dul;

					quatd Q1 = pprb->GetRotation();
					
					dr = r0 + dr - r1;

					// grab the parent's local displacements
					vec3d dR = vec3d(dul[0], dul[1], dul[2]);
					v = vec3d(dul[3], dul[4], dul[5]);
					w = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
					quatd dQ = quatd(w, v);

					dQ = Q1*dQ*Q1.Inverse();

					// update global displacements
					quatd Qi = Q1.Inverse();
					dr = dR + r1 + dQ*dr - r0;
					dq = dQ*dq;

					// move up in the chain
					pprb = pprb->m_prb;
					Q0 = Q1;
				}

				// set global displacements
				double* du = RB.m_du;

				du[0] = dr.x;
				du[1] = dr.y;
				du[2] = dr.z;

				v = dq.GetVector();
				w = dq.GetAngle();
				du[3] = w*v.x;
				du[4] = w*v.y;
				du[5] = w*v.z;
			}
		}
	}

	// store rigid displacements in Ui vector
	for (i=0; i<NO; ++i)
	{
		FERigidBody& RB = *rigid.Object(i);
		for (j=0; j<6; ++j)
		{
			int I = -RB.m_LM[j]-2;
			if (I >= 0) ui[I] = RB.m_du[j];
		}
	}

	// intialize material point data
	for (i=0; i<mesh.Domains(); ++i) mesh.Domain(i).PreSolveUpdate(tp);

	fem.Update();
}

//-----------------------------------------------------------------------------
//! calculates the concentrated nodal forces

void FEExplicitSolidSolver::NodalForces(FEGlobalVector& R, const FETimeInfo& tp)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// loop over nodal forces
	int ncnf = fem.NodalLoads();
	for (int i=0; i<ncnf; ++i)
	{
		FENodalLoad& fc = *fem.NodalLoad(i);
		if (fc.IsActive()) fc.Residual(R, tp);
	}
}

//-----------------------------------------------------------------------------
bool FEExplicitSolidSolver::DoSolve()
{
	int i, n;

	vector<double> u0(m_neq);
	vector<double> Rold(m_neq);

	// Get the current step
	FEModel& fem = *GetFEModel();
	FEAnalysis* pstep = fem.GetCurrentStep();

	// prepare for the first iteration
	PrepStep();

	feLog(" %d\n", m_niter+1);

	// get the mesh
	FEMesh& mesh = fem.GetMesh();
	int N = mesh.Nodes(); // this is the total number of nodes in the mesh
	int j,iel;
    double dt = fem.GetTime().timeIncrement;
	double avx = 0.0;  // average element velocity in each direction
	double avy = 0.0;
	double avz = 0.0;
	double mass_at_node;

	for (i=0; i<N; ++i) // zero the new acceleration vector ready to add in the damping components
	{
		FENode& node = mesh.Node(i);
		node.m_at.x = 0.0;
		node.m_at.y = 0.0;
		node.m_at.z = 0.0;
	}

	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(nd));
		if (pbd)  // it is an elastic solid domain
		{
			double ** emass = domain_mass[nd]; // array of pointers to the element mass records for this domain
			for (iel=0; iel<pbd->Elements(); ++iel)
			{
				FESolidElement& el = pbd->Element(iel);

				// zero the three average velocity components
				avx = 0.0;
				avy = 0.0;
				avz = 0.0;

				// will use previously calculated element mass data for weighted averaging of velocities

				// loop over each element to find the average velocity
				// then calculate the weighted velocity change for each node
				// add each velocity change into node.m_vt
				double * this_element = emass[iel]; // pointer to the array of fractional nodal masses for this element
				for (j=0; j<el.Nodes(); j++) // loop over each node in the element
				{
					FENode& node = mesh.Node(el.m_node[j]);  // get the node 
					avx += node.m_vp.x*this_element[j+1];  // add each of the three components to the averages
					avy += node.m_vp.y*this_element[j+1];  // weighted by the fractional mass of the node
					avz += node.m_vp.z*this_element[j+1];  // remembering that this_element[0] is the total mass
				}
				for (j=0; j<el.Nodes(); j++) // loop over each node in the element again
				// and calculate and add in the velocity change contribution to each dof
				{
					FENode& node = mesh.Node(el.m_node[j]);  // get the node 
					//	need to find node.m_vt.x += (avx-node.m_vp.x)*dt*m_dyn_damping*element_mass_at_node/total_mass at node;
					// should be t* = dt/(h/c) not dt
					// put this into the accelerations as (avx-node.m_vp.x)*m_dyn_damping*element_mass_at_node
					// then it will be multiplied by dt and divided by m_inv_mass later 
					mass_at_node = this_element[j+1]*this_element[0];
					node.m_at.x += (avx-node.m_vp.x)*mass_at_node*m_dyn_damping;
					node.m_at.y += (avy-node.m_vp.y)*mass_at_node*m_dyn_damping;
					node.m_at.z += (avz-node.m_vp.z)*mass_at_node*m_dyn_damping;
				}
			}  // loop over elements
		}  // if (pbd)
	}  // loop over domains

	for (i=0; i<N; ++i)
	{
		FENode& node = mesh.Node(i);
		//  calculate acceleration using F=ma and update - note m_inv_mass is 1/m so multiply not divide
		n=(int)m_R1[0];
		if ((n = node.m_ID[m_dofU[0]]) >= 0) node.m_at.x = (node.m_at.x+m_R1[n])*m_inv_mass[n];
		if ((n = node.m_ID[m_dofU[1]]) >= 0) node.m_at.y = (node.m_at.y+m_R1[n])*m_inv_mass[n];
		if ((n = node.m_ID[m_dofU[2]]) >= 0) node.m_at.z = (node.m_at.z+m_R1[n])*m_inv_mass[n];
		// and update the velocities using the accelerations
		// which are added to the previously calculated velocity changes from damping
		vec3d vt = node.m_vp + node.m_at*dt;
		node.set_vec3d(m_dofV[0], m_dofV[1], m_dofV[2], vt);	//  update velocity using acceleration m_at
		//	calculate incremental displacement using the velocity
		if ((n = node.m_ID[m_dofU[0]]) >= 0) m_ui[n] = vt.x*dt;
		if ((n = node.m_ID[m_dofU[1]]) >= 0) m_ui[n] = vt.y*dt;
		if ((n = node.m_ID[m_dofU[2]]) >= 0) m_ui[n] = vt.z*dt;
	}

	// need to update everything for the explicit solver
	// Update geometry
	Update(m_ui);

	// calculate new residual at this point - which will be used on the next step to find the acceleration
	Residual(m_R1);

	// update total displacements
	int neq = (int)m_Ui.size();
	for (i=0; i<neq; ++i) m_Ui[i] += m_ui[i];

	// increase iteration number
	m_niter++;

	// do minor iterations callbacks
	fem.DoCallback(CB_MINOR_ITERS);

	// if converged we update the total displacements
	m_Ut += m_Ui;

	return true;
}

//-----------------------------------------------------------------------------
//! calculates the residual vector
//! Note that the concentrated nodal forces are not calculated here.
//! This is because they do not depend on the geometry 
//! so we only calculate them once (in Quasin) and then add them here.

bool FEExplicitSolidSolver::Residual(vector<double>& R)
{
	int i;

	// get the time information
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	const FETimeInfo& tp = fem.GetTime();

	// initialize residual with concentrated nodal loads
	R = m_Fn;

	// zero nodal reaction forces
	zero(m_Fr);

	// setup the global vector
	FEGlobalVector RHS(fem, R, m_Fr);

	// zero rigid body reaction forces
	FERigidSystem& rigid = *fem.GetRigidSystem();
	int NRB = rigid.Objects();
	for (i=0; i<NRB; ++i)
	{
		FERigidBody& RB = *rigid.Object(i);
		RB.m_Fr = RB.m_Mr = vec3d(0,0,0);
	}

	// get the mesh
	FEMesh& mesh = fem.GetMesh();

	// calculate the internal (stress) forces
	for (i=0; i<mesh.Domains(); ++i)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
		dom.InternalForces(RHS);
	}

	// calculate the body forces
	for (int j = 0; j<fem.BodyLoads(); ++j)
	{
		FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.GetBodyLoad(j));
		if (pbf && pbf->IsActive())
		{
			for (i = 0; i<pbf->Domains(); ++i)
			{
				FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(*pbf->Domain(i));
				if (pbf) dom.BodyForce(RHS, *pbf);
			}
		}
	}

	// calculate inertial forces for dynamic problems
	if (fem.GetCurrentStep()->m_nanalysis == FE_DYNAMIC) InertialForces(RHS);

	// calculate forces due to surface loads
	int nsl = fem.SurfaceLoads();
	for (i=0; i<nsl; ++i)
	{
		FESurfaceLoad* psl = fem.SurfaceLoad(i);
		if (psl->IsActive()) psl->Residual(RHS, tp);
	}

	// calculate contact forces
	if (fem.SurfacePairConstraints() > 0)
	{
		ContactForces(RHS);
	}

	// calculate nonlinear constraint forces
	// note that these are the linear constraints
	// enforced using the augmented lagrangian
	NonLinearConstraintForces(RHS, tp);

	// forces due to point constraints
//	for (i=0; i<(int) fem.m_PC.size(); ++i) fem.m_PC[i]->Residual(this, R);

	// set the nodal reaction forces
	// TODO: Is this a good place to do this?
	for (i=0; i<mesh.Nodes(); ++i)
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
//! Calculates the contact forces
void FEExplicitSolidSolver::ContactForces(FEGlobalVector& R)
{
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();
	for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));
		pci->Residual(R, tp);
	}
}


//-----------------------------------------------------------------------------
//! calculate the nonlinear constraint forces 
void FEExplicitSolidSolver::NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp)
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
//! This function calculates the inertial forces for dynamic problems

void FEExplicitSolidSolver::InertialForces(FEGlobalVector& R)
{
	// get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// allocate F
	vector<double> F(3*mesh.Nodes());
	zero(F);

	// calculate F
    double dt = fem.GetTime().timeIncrement;
	double a = 4.0 / dt;
	double b = a / dt;
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		vec3d& rt = node.m_rt;
		vec3d& rp = node.m_rp;
		vec3d& vp = node.m_vp;
		vec3d& ap = node.m_ap;

		F[3*i  ] = b*(rt.x - rp.x) - a*vp.x - ap.x;
		F[3*i+1] = b*(rt.y - rp.y) - a*vp.y - ap.y;
		F[3*i+2] = b*(rt.z - rp.z) - a*vp.z - ap.z;
	}

	// now multiply F with the mass matrix
	matrix ke;
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(nd));
		dom.InertialForces(R, F);
	}
}
