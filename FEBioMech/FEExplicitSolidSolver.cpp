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
#include "FEExplicitSolidSolver.h"
#include "FEElasticSolidDomain.h"
#include "FEElasticShellDomain.h"
#include "FELinearTrussDomain.h"
#include "FEElasticTrussDomain.h"
#include "FERigidMaterial.h"
#include "FEBodyForce.h"
#include "FEContactInterface.h"
#include "FERigidBody.h"
#include "RigidBC.h"
#include "FEMechModel.h"
#include <FECore/FENodalLoad.h>
#include <FECore/FESurfaceLoad.h>
#include "FECore/log.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FELinearConstraintManager.h>
#include <FECore/FENLConstraint.h>
#include "FEResidualVector.h"
#include "FEBioMech.h"
#include "FESolidAnalysis.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_FECORE_CLASS(FEExplicitSolidSolver, FESolver)
	ADD_PARAMETER(m_mass_lumping, "mass_lumping");
	ADD_PARAMETER(m_dyn_damping, "dyn_damping");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEExplicitSolidSolver::FEExplicitSolidSolver(FEModel* pfem) : 
	FESolver(pfem), 
	m_dofU(pfem), m_dofV(pfem), m_dofQ(pfem), m_dofRQ(pfem),
	m_dofSU(pfem), m_dofSV(pfem), m_dofSA(pfem),
	m_rigidSolver(pfem)
{
	m_dyn_damping = 1;
	m_niter = 0;
	m_nreq = 0;

	m_mass_lumping = HRZ_LUMPING;

	// Allocate degrees of freedom
	// TODO: Can this be done in Init, since there is no error checking
	if (pfem)
	{
		DOFS& dofs = pfem->GetDOFS();
		int varD = dofs.AddVariable("displacement", VAR_VEC3);
		dofs.SetDOFName(varD, 0, "x");
		dofs.SetDOFName(varD, 1, "y");
		dofs.SetDOFName(varD, 2, "z");
		int varQ = dofs.AddVariable("rotation", VAR_VEC3);
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
		int varSU = dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_DISPLACEMENT), VAR_VEC3);
		dofs.SetDOFName(varSU, 0, "sx");
		dofs.SetDOFName(varSU, 1, "sy");
		dofs.SetDOFName(varSU, 2, "sz");
		int varSV = dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_VELOCITY), VAR_VEC3);
		dofs.SetDOFName(varSV, 0, "svx");
		dofs.SetDOFName(varSV, 1, "svy");
		dofs.SetDOFName(varSV, 2, "svz");
		int varSA = dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_ACCELERATION), VAR_VEC3);
		dofs.SetDOFName(varSA, 0, "sax");
		dofs.SetDOFName(varSA, 1, "say");
		dofs.SetDOFName(varSA, 2, "saz");

		// get the DOF indices
		m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
		m_dofV.AddVariable(FEBioMech::GetVariableName(FEBioMech::VELOCITY));
		m_dofQ.AddVariable(FEBioMech::GetVariableName(FEBioMech::ROTATION));
		m_dofRQ.AddVariable(FEBioMech::GetVariableName(FEBioMech::RIGID_ROTATION));
		m_dofSU.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_DISPLACEMENT));
		m_dofSV.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_VELOCITY));
		m_dofSA.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_ACCELERATION));
	}
}

//-----------------------------------------------------------------------------
void FEExplicitSolidSolver::Clean()
{
}

//-----------------------------------------------------------------------------
bool FEExplicitSolidSolver::CalculateMassMatrix()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	vector<double> massVector(m_neq), dummy(m_neq);
	FEGlobalVector M(fem, massVector, dummy);
	matrix me;
	vector <int> lm;
	vector <double> el_lumped_mass;

	// loop over all domains
	if (m_mass_lumping == NO_MASS_LUMPING)
	{
		// use consistent mass matrix.
		// TODO: implement this
		assert(false);
		return false;
	}
	else if (m_mass_lumping == ROW_SUM_LUMPING)
	{
		for (int nd = 0; nd < mesh.Domains(); ++nd)
		{
			// check whether it is a solid domain
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(nd));
			if (pbd)  // it is an elastic solid domain
			{
				FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(pbd->GetMaterial());

				// loop over all the elements
				for (int iel = 0; iel < pbd->Elements(); ++iel)
				{
					FESolidElement& el = pbd->Element(iel);
					pbd->UnpackLM(el, lm);

					int nint = el.GaussPoints();
					int neln = el.Nodes();

					me.resize(neln, neln);
					me.zero();

					// create the element mass matrix
					for (int n = 0; n < nint; ++n)
					{
						FEMaterialPoint& mp = *el.GetMaterialPoint(n);
						double d = pme->Density(mp);
						double detJ0 = pbd->detJ0(el, n)*el.GaussWeights()[n];

						double* H = el.H(n);
						for (int i = 0; i < neln; ++i)
							for (int j = 0; j < neln; ++j)
							{
								double kab = H[i] * H[j] * detJ0*d;
								me[i][j] += kab;
							}
					}

					// reduce to a lumped mass vector and add up the total
					el_lumped_mass.assign(3 * neln, 0.0);
					for (int i = 0; i < neln; ++i)
					{
						for (int j = 0; j < neln; ++j)
						{
							double kab = me[i][j];
							el_lumped_mass[3 * i] += kab;
							el_lumped_mass[3 * i + 1] += kab;
							el_lumped_mass[3 * i + 2] += kab;
						}
					}

					// assemble element matrix into inv_mass vector 
					M.Assemble(el.m_node, lm, el_lumped_mass);
				} // loop over elements
			}
			else if (dynamic_cast<FEElasticShellDomain*>(&mesh.Domain(nd)))
			{
				FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&mesh.Domain(nd));
				FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(psd->GetMaterial());
				// loop over all the elements
				for (int iel = 0; iel < psd->Elements(); ++iel)
				{
					FEShellElement& el = psd->Element(iel);
					psd->UnpackLM(el, lm);

					// create the element's stiffness matrix
					FEElementMatrix ke(el);
					int neln = el.Nodes();
					int ndof = 6 * el.Nodes();
					ke.resize(ndof, ndof);
					ke.zero();

					// calculate inertial stiffness
					psd->ElementMassMatrix(el, ke, 1.0);

					// reduce to a lumped mass vector and add up the total
					el_lumped_mass.assign(ndof, 0.0);
					for (int i = 0; i < ndof; ++i)
					{
						for (int j = 0; j < ndof; ++j)
						{
							double kab = ke[i][j];
							el_lumped_mass[i] += kab;
						}
					}
					// assemble element matrix into inv_mass vector 
					M.Assemble(el.m_node, lm, el_lumped_mass);
				}
			}
			else if (dynamic_cast<FELinearTrussDomain*>(&mesh.Domain(nd)))
			{
				FELinearTrussDomain* ptd = dynamic_cast<FELinearTrussDomain*>(&mesh.Domain(nd));

				// loop over all the elements
				for (int iel = 0; iel < ptd->Elements(); ++iel)
				{
					FETrussElement& el = ptd->Element(iel);
					ptd->UnpackLM(el, lm);

					// create the element's stiffness matrix
					FEElementMatrix ke(el);
					int neln = el.Nodes();
					ke.resize(neln, neln);
					ke.zero();

					// calculate inertial stiffness
					ptd->ElementMassMatrix(el, ke);

					// reduce to a lumped mass vector and add up the total
					el_lumped_mass.assign(3 * neln, 0.0);
					for (int i = 0; i < neln; ++i)
					{
						for (int j = 0; j < neln; ++j)
						{
							double kab = ke[i][j];
							el_lumped_mass[3 * i    ] += kab;
							el_lumped_mass[3 * i + 1] += kab;
							el_lumped_mass[3 * i + 2] += kab;
						}
					}

					// assemble element matrix into inv_mass vector 
					M.Assemble(el.m_node, lm, el_lumped_mass);
				}
			}
			else if (dynamic_cast<FEElasticTrussDomain*>(&mesh.Domain(nd)))
			{
				FEElasticTrussDomain* ptd = dynamic_cast<FEElasticTrussDomain*>(&mesh.Domain(nd));

				// loop over all the elements
				for (int iel = 0; iel < ptd->Elements(); ++iel)
				{
					FETrussElement& el = ptd->Element(iel);
					ptd->UnpackLM(el, lm);

					// create the element's stiffness matrix
					FEElementMatrix ke(el);
					int neln = el.Nodes();
					ke.resize(neln, neln);
					ke.zero();

					// calculate inertial stiffness
					ptd->ElementMassMatrix(el, ke);

					// reduce to a lumped mass vector and add up the total
					el_lumped_mass.assign(3 * neln, 0.0);
					for (int i = 0; i < neln; ++i)
					{
						for (int j = 0; j < neln; ++j)
						{
							double kab = ke[i][j];
							el_lumped_mass[3 * i    ] += kab;
							el_lumped_mass[3 * i + 1] += kab;
							el_lumped_mass[3 * i + 2] += kab;
						}
					}

					// assemble element matrix into inv_mass vector 
					M.Assemble(el.m_node, lm, el_lumped_mass);
				}
			}
			else
			{
//				return false;
			}
		}
	}
	else if (m_mass_lumping == HRZ_LUMPING)
	{
		for (int nd = 0; nd < mesh.Domains(); ++nd)
		{
			// check whether it is a solid domain
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(nd));
			if (pbd)  // it is an elastic solid domain
			{
				FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(pbd->GetMaterial());

				// loop over all the elements
				for (int iel = 0; iel < pbd->Elements(); ++iel)
				{
					FESolidElement& el = pbd->Element(iel);
					pbd->UnpackLM(el, lm);

					int nint = el.GaussPoints();
					int neln = el.Nodes();

					me.resize(neln, neln);
					me.zero();

					// calculate the element mass matrix (and element mass).
					double Me = 0.0;
					double* w = el.GaussWeights();
					for (int n = 0; n < nint; ++n)
					{
						FEMaterialPoint& mp = *el.GetMaterialPoint(n);
						double d = pme->Density(mp);
						double detJ0 = pbd->detJ0(el, n)*el.GaussWeights()[n];
						Me += d * detJ0 * w[n];

						double* H = el.H(n);
						for (int i = 0; i < neln; ++i)
							for (int j = 0; j < neln; ++j)
							{
								double kab = H[i] * H[j] * detJ0*d;
								me[i][j] += kab;
							}
					}

					// calculate sum of diagonals
					double S = 0.0;
					for (int i = 0; i < neln; ++i) S += me[i][i];

					// reduce to a lumped mass vector and add up the total
					el_lumped_mass.assign(3 * neln, 0.0);
					for (int i = 0; i < neln; ++i)
					{
						double mab = me[i][i] * Me / S;
						el_lumped_mass[3 * i    ] = mab;
						el_lumped_mass[3 * i + 1] = mab;
						el_lumped_mass[3 * i + 2] = mab;
					}

					// assemble element matrix into inv_mass vector 
					M.Assemble(el.m_node, lm, el_lumped_mass);
				} // loop over elements
			}
			else if(dynamic_cast<FEElasticShellDomain*>(&mesh.Domain(nd)))
			{
				FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&mesh.Domain(nd));
				FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(psd->GetMaterial());
				// loop over all the elements
				for (int iel = 0; iel < psd->Elements(); ++iel)
				{
					FEShellElement& el = psd->Element(iel);
					psd->UnpackLM(el, lm);

					// create the element's stiffness matrix
					FEElementMatrix ke(el);
					int neln = el.Nodes();
					int ndof = 6 * el.Nodes();
					ke.resize(ndof, ndof);
					ke.zero();

					// calculate inertial stiffness
					psd->ElementMassMatrix(el, ke, 1.0);

					// calculate the element mass
					double Me = 0.0;
					double* w = el.GaussWeights();
					for (int n = 0; n < el.GaussPoints(); ++n)
					{
						FEMaterialPoint& mp = *el.GetMaterialPoint(n);
						double d = pme->Density(mp);
						double detJ0 = psd->detJ0(el, n) * el.GaussWeights()[n];
						Me += d * detJ0 * w[n];
					}

					// calculate sum of diagonals
					double S = 0.0;
					for (int i = 0; i < ndof; ++i) S += ke[i][i] / 3.0;

					// reduce to a lumped mass vector and add up the total
					el_lumped_mass.assign(ndof, 0.0);
					for (int i = 0; i < ndof; ++i)
					{
						double mab = ke[i][i] * Me / S;
						el_lumped_mass[i] = mab;
					}

					// assemble element matrix into inv_mass vector 
					M.Assemble(el.m_node, lm, el_lumped_mass);
				}
			}
			else
			{
				// TODO: we can only do solid domains right now.
				return false;
			}
		}
	}
	else
	{
		assert(false);
		return false;
	}

	// we need the inverse of the lumped masses later
	// Also, make sure the lumped masses are positive.
	for (int i = 0; i < massVector.size(); ++i)
	{
		m_data[i].mi = 0;
		if (massVector[i] != 0.0) m_data[i].mi = 1.0 / massVector[i];
	}

	return true;
}

//-----------------------------------------------------------------------------
//! initialize equations
bool FEExplicitSolidSolver::InitEquations()
{
	if (FESolver::InitEquations() == false) return false;

	// allocate rigid body equation numbers
	m_neq = m_rigidSolver.InitEquations(m_neq);

	return true;
}

//-----------------------------------------------------------------------------
bool FEExplicitSolidSolver::Init()
{
	if (FESolver::Init() == false) return false;

	FEMechModel& fem = dynamic_cast<FEMechModel&>(*GetFEModel());
	FEMesh& mesh = fem.GetMesh();

	// set the dynamic update flag only if we are running a dynamic analysis
	// NOTE: I don't think we need to set the dynamic update flag, in fact, we should 
	//		 turn it off, since it's on by default, and it incurs a significant performance overhead
//	bool b = (fem.GetCurrentStep()->m_nanalysis == FESolidAnalysis::DYNAMIC ? true : false);
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEElasticSolidDomain* d = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(i));
		FEElasticShellDomain* s = dynamic_cast<FEElasticShellDomain*>(&mesh.Domain(i));
		if (d) d->SetDynamicUpdateFlag(false);
		if (s) s->SetDynamicUpdateFlag(false);
	}

	// get nr of equations
	int neq = m_neq;

	// allocate vectors
	m_data.resize(neq);
	m_Rt.assign(neq, 0);
	m_Fr.assign(neq, 0);
	m_Ut.assign(neq, 0);
	m_ui.assign(neq, 0);

	fem.Update();

	// we need to fill the total displacement vector m_Ut
	gather(m_Ut, mesh, m_dofU[0]);
	gather(m_Ut, mesh, m_dofU[1]);
	gather(m_Ut, mesh, m_dofU[2]);
	gather(m_Ut, mesh, m_dofQ[0]);
	gather(m_Ut, mesh, m_dofQ[1]);
	gather(m_Ut, mesh, m_dofQ[2]);
	gather(m_Ut, mesh, m_dofSU[0]);
	gather(m_Ut, mesh, m_dofSU[1]);
	gather(m_Ut, mesh, m_dofSU[2]);

	// calculate the inverse mass vector for the explicit analysis
	if (CalculateMassMatrix() == false)
	{
		feLogError("Failed building mass matrix.");
		return false;
	}

	// calculate the initial acceleration
	// (Only when the totiter == 0, in case of a restart)
	if (fem.GetCurrentStep()->m_ntotiter == 0)
	{
		// Calculate initial residual to be used on the first time step
		if (Residual(m_Rt) == false) return false;

		for (int i = 0; i < mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(i);
			int n;
			if ((n = node.m_ID[m_dofU[0]]) >= 0) { m_data[n].a = node.m_at.x = m_Rt[n] * m_data[n].mi; }
			if ((n = node.m_ID[m_dofU[1]]) >= 0) { m_data[n].a = node.m_at.y = m_Rt[n] * m_data[n].mi; }
			if ((n = node.m_ID[m_dofU[2]]) >= 0) { m_data[n].a = node.m_at.z = m_Rt[n] * m_data[n].mi; }

			if ((n = node.m_ID[m_dofSU[0]]) >= 0) { m_data[n].a = m_Rt[n] * m_data[n].mi; node.set(m_dofSA[0], m_data[n].a);}
			if ((n = node.m_ID[m_dofSU[1]]) >= 0) { m_data[n].a = m_Rt[n] * m_data[n].mi; node.set(m_dofSA[1], m_data[n].a);}
			if ((n = node.m_ID[m_dofSU[2]]) >= 0) { m_data[n].a = m_Rt[n] * m_data[n].mi; node.set(m_dofSA[2], m_data[n].a);}
		}

		// do rigid bodies
		for (int i = 0; i < fem.RigidBodies(); ++i)
		{
			FERigidBody& rb = *fem.GetRigidBody(i);
			vec3d Fn, Mn;
			int n;
			if ((n = rb.m_LM[0]) >= 0) { Fn.x = m_Rt[n]; }
			if ((n = rb.m_LM[1]) >= 0) { Fn.y = m_Rt[n]; }
			if ((n = rb.m_LM[2]) >= 0) { Fn.z = m_Rt[n]; }
			if ((n = rb.m_LM[3]) >= 0) { Mn.x = m_Rt[n]; }
			if ((n = rb.m_LM[4]) >= 0) { Mn.y = m_Rt[n]; }
			if ((n = rb.m_LM[5]) >= 0) { Mn.z = m_Rt[n]; }

			rb.m_at = Fn / rb.m_mass;

			// TODO: angular acceleration
		}
	}

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
#pragma omp parallel for
	for (int i=0; i<m_Ut.size(); ++i) U[i] = ui[i] + m_Ut[i];

	// update flexible nodes
	// translational dofs
	scatter3(U, mesh, m_dofU[0], m_dofU[1], m_dofU[2]);

	// rotational dofs
	// TODO: Commenting this out, since this is only needed for the old shells, which I'm not sure
	//       if they would work with the explicit solver anyways.
//	scatter(U, mesh, m_dofQ[0]);
//	scatter(U, mesh, m_dofQ[1]);
//	scatter(U, mesh, m_dofQ[2]);

	// shell displacement
	scatter3(U, mesh, m_dofSU[0], m_dofSU[1], m_dofSU[2]);

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
#pragma omp parallel for
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
	const int NRB = fem.RigidBodies();
	if (NRB == 0) return;

	// first calculate the rigid body displacement increments
	for (int i=0; i<NRB; ++i)
	{
		// get the rigid body
		FERigidBody& RB = *fem.GetRigidBody(i);
		int *lm = RB.m_LM;
		double* du = RB.m_du;

		if (RB.m_prb == 0)
		{
			for (int j=0; j<6; ++j)
			{
				du[j] = (lm[j] >=0 ? ui[lm[j]] : 0);
			}
		}
	}

	// for prescribed displacements, the displacement increments are evaluated differently
	// TODO: Is this really necessary? Why can't the ui vector contain the correct values?
	for (int i = 0; i < NRB; ++i)
	{
		FERigidBody& RB = *fem.GetRigidBody(i);
		for (int j = 0; j < 6; ++j)
		{
			FERigidPrescribedBC* dc = RB.m_pDC[j];
			if (dc)
			{
				int I = dc->GetBC();
				RB.m_du[I] = dc->Value() - RB.m_Up[I];
			}
		}
	}

	// update the rigid bodies
	for (int i=0; i<NRB; ++i)
	{
		// get the rigid body
		FERigidBody& RB = *fem.GetRigidBody(i);
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
	fem.UpdateRigidMesh();

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
	if (ar.IsShallow() == false)
	{
		ar & m_Ut & m_Rt;

		if (ar.IsSaving())
		{
			int N = (int)m_data.size();
			ar << N;
			for (int i = 0; i < N; ++i)
			{
				Data& d = m_data[i];
				ar << d.v << d.a << d.mi;
			}
		}
		else if (ar.IsLoading())
		{
			int N = 0;
			ar >> N;
			m_data.resize(N);
			for (int i = 0; i < N; ++i)
			{
				Data& d = m_data[i];
				ar >> d.v >> d.a >> d.mi;
			}
		}
	}
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

	// store previous mesh state
	// we need them for velocity and acceleration calculations
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FEMesh& mesh = fem.GetMesh();
#pragma omp parallel for
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& ni = mesh.Node(i);
		ni.m_rp = ni.m_rt;
		ni.m_vp = ni.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);
		ni.m_ap = ni.m_at;
		ni.UpdateValues();
	}

	const FETimeInfo& tp = fem.GetTime();

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
	int NO = fem.RigidBodies();
	for (i=0; i<NO; ++i) fem.GetRigidBody(i)->Init();

	// calculate local rigid displacements
	for (i = 0; i < NO; ++i)
	{
		FERigidBody& rb = *fem.GetRigidBody(i);
		for (int j = 0; j < 6; ++j)
		{
			FERigidPrescribedBC* dc = rb.m_pDC[j];
			if (dc)
			{
				int I = dc->GetBC();
				rb.m_dul[I] = dc->Value() - rb.m_Ut[I];
			}
		}
	}

	// calculate global rigid displacements
	for (i=0; i<NO; ++i)
	{
		FERigidBody* prb = fem.GetRigidBody(i);
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
		FERigidBody& RB = *fem.GetRigidBody(i);
		for (j=0; j<6; ++j)
		{
			int I = -RB.m_LM[j]-2;
			if (I >= 0) ui[I] = RB.m_du[j];
		}
	}

	// intialize material point data
	for (i=0; i<mesh.Domains(); ++i) mesh.Domain(i).PreSolveUpdate(tp);

	// NOTE: Commenting this out since I don't think anything needs to be updated here. 
	//       This also halves the number of update calls, so gives a performance boost. 
//	fem.Update();
}

//-----------------------------------------------------------------------------
bool FEExplicitSolidSolver::DoSolve()
{
	// Get the current step
	FEMechModel& fem = dynamic_cast<FEMechModel&>(*GetFEModel());
	FEAnalysis* pstep = fem.GetCurrentStep();

	// prepare for solve
	PrepStep();

//	feLog(" %d\n", m_niter+1);

	// get the mesh
	FEMesh& mesh = fem.GetMesh();
	int N = mesh.Nodes(); // this is the total number of nodes in the mesh
	double dt = fem.GetTime().timeIncrement;

	// collect accelerations, velocities, displacements
	// NOTE: I don't think this is necessary
/*
#pragma omp parallel for shared(mesh)
	for (int i = 0; i < mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		vec3d vt = node.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);
		int n;
		if ((n = node.m_ID[m_dofU[0]]) >= 0) { m_vn[n] = vt.x; m_an[n] = node.m_at.x; }
		if ((n = node.m_ID[m_dofU[1]]) >= 0) { m_vn[n] = vt.y; m_an[n] = node.m_at.y; }
		if ((n = node.m_ID[m_dofU[2]]) >= 0) { m_vn[n] = vt.z; m_an[n] = node.m_at.z; }

		if ((n = node.m_ID[m_dofSU[0]]) >= 0) { m_vn[n] = node.get(m_dofSV[0]); m_an[n] = node.get(m_dofSA[0]); }
		if ((n = node.m_ID[m_dofSU[1]]) >= 0) { m_vn[n] = node.get(m_dofSV[1]); m_an[n] = node.get(m_dofSA[1]); }
		if ((n = node.m_ID[m_dofSU[2]]) >= 0) { m_vn[n] = node.get(m_dofSV[2]); m_an[n] = node.get(m_dofSA[2]); }
	}
*/

	// do rigid bodies
	for (int i = 0; i < fem.RigidBodies(); ++i)
	{
		FERigidBody& rb = *fem.GetRigidBody(i);
		int n;
		if ((n = rb.m_LM[0]) >= 0) { m_data[n].v = rb.m_vt.x; m_data[n].a = rb.m_at.x; }
		if ((n = rb.m_LM[1]) >= 0) { m_data[n].v = rb.m_vt.y; m_data[n].a = rb.m_at.y; }
		if ((n = rb.m_LM[2]) >= 0) { m_data[n].v = rb.m_vt.z; m_data[n].a = rb.m_at.z; }
		
		// convert to rigid frame
		quatd Q = rb.GetRotation();
		quatd Qi = Q.Inverse();
		vec3d Wn = Qi * rb.m_wt;
		vec3d An = Qi * rb.m_alt;
		if ((n = rb.m_LM[3]) >= 0) { m_data[n].v = Wn.x; m_data[n].a = An.x; }
		if ((n = rb.m_LM[4]) >= 0) { m_data[n].v = Wn.y; m_data[n].a = An.y; }
		if ((n = rb.m_LM[5]) >= 0) { m_data[n].v = Wn.z; m_data[n].a = An.z; }
	}

	double Dnorm = 0.0;
#pragma omp parallel for reduction(+: Dnorm)
	for (int i = 0; i < m_neq; ++i)
	{
		// velocity predictor
		m_data[i].v += m_data[i].a * dt*0.5;

		// update displacements
		m_ui[i] = dt * m_data[i].v;

		// update norm
		Dnorm += m_ui[i] * m_ui[i];
	}
	Dnorm = sqrt(Dnorm);
	feLog("\t displacement norm : %lg\n", Dnorm);

	// the update is done in the spatial frame, so we need to update
	// rigid body rotation increment
	for (int i = 0; i < fem.RigidBodies(); ++i)
	{
		FERigidBody& rb = *fem.GetRigidBody(i);
		quatd Q = rb.GetRotation();
		vec3d dq(0, 0, 0);
		int n;
		if ((n = rb.m_LM[3]) >= 0) { dq.x = m_ui[n]; }
		if ((n = rb.m_LM[4]) >= 0) { dq.y = m_ui[n]; }
		if ((n = rb.m_LM[5]) >= 0) { dq.z = m_ui[n]; }
		vec3d qt = Q * dq;
		if ((n = rb.m_LM[3]) >= 0) { m_ui[n] = qt.x; }
		if ((n = rb.m_LM[4]) >= 0) { m_ui[n] = qt.y; }
		if ((n = rb.m_LM[5]) >= 0) { m_ui[n] = qt.z; }
	}

	Update(m_ui);

	// evaluate acceleration
	Residual(m_Rt);

	double Rnorm = 0.0;
#pragma omp parallel shared(Rnorm)
	{
#pragma omp for reduction(+: Rnorm) nowait
		for (int i = 0; i < m_neq; ++i)
		{
			Rnorm += m_Rt[i] * m_Rt[i];
		}

#pragma omp for nowait
		for (int i = 0; i < m_neq; ++i)
		{
			// update total displacement
			m_Ut[i] += m_ui[i];
		}

#pragma omp for
		for (int i = 0; i < m_neq; ++i)
		{
			m_data[i].a = m_Rt[i] * m_data[i].mi;

			// update velocity
			m_data[i].v = m_dyn_damping * (m_data[i].v + m_data[i].a * dt * 0.5);
		}

		// scatter velocity and accelerations
#pragma omp for nowait
		for (int i = 0; i < mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(i);
			int n;
			if ((n = node.m_ID[m_dofU[0]]) >= 0) { node.set(m_dofV[0], m_data[n].v); node.m_at.x = m_data[n].a; }
			if ((n = node.m_ID[m_dofU[1]]) >= 0) { node.set(m_dofV[1], m_data[n].v); node.m_at.y = m_data[n].a; }
			if ((n = node.m_ID[m_dofU[2]]) >= 0) { node.set(m_dofV[2], m_data[n].v); node.m_at.z = m_data[n].a; }

			if ((n = node.m_ID[m_dofSU[0]]) >= 0) { node.set(m_dofSV[0], m_data[n].v); node.set(m_dofSA[0], m_data[n].a); }
			if ((n = node.m_ID[m_dofSU[1]]) >= 0) { node.set(m_dofSV[1], m_data[n].v); node.set(m_dofSA[1], m_data[n].a); }
			if ((n = node.m_ID[m_dofSU[2]]) >= 0) { node.set(m_dofSV[2], m_data[n].v); node.set(m_dofSA[2], m_data[n].a); }
		}
	}

	Rnorm = sqrt(Rnorm);
	feLog("\t force vector norm : %lg\n", Rnorm);

	// do rigid bodies
	for (int i = 0; i < fem.RigidBodies(); ++i)
	{
		FERigidBody& rb = *fem.GetRigidBody(i);
		quatd Q = rb.GetRotation();
		quatd Qi = Q.Inverse();

		// get the force and moment
		vec3d Fn(0, 0, 0), Mn;
		int n;
		if ((n = rb.m_LM[0]) >= 0) Fn.x = m_Rt[n];
		if ((n = rb.m_LM[1]) >= 0) Fn.y = m_Rt[n];
		if ((n = rb.m_LM[2]) >= 0) Fn.z = m_Rt[n];
		if ((n = rb.m_LM[3]) >= 0) Mn.x = m_Rt[n];
		if ((n = rb.m_LM[4]) >= 0) Mn.y = m_Rt[n];
		if ((n = rb.m_LM[5]) >= 0) Mn.z = m_Rt[n];

		// convert to rigid frame
		Mn = Qi * Mn;

		// linear momentum update
		vec3d An = Fn / rb.m_mass;
		rb.m_at = Q*An;

		vec3d Vn(0,0,0);
		if ((n = rb.m_LM[0]) >= 0) Vn.x = m_dyn_damping * (m_data[n].v + An.x * dt * 0.5);
		if ((n = rb.m_LM[1]) >= 0) Vn.y = m_dyn_damping * (m_data[n].v + An.y * dt * 0.5);
		if ((n = rb.m_LM[2]) >= 0) Vn.z = m_dyn_damping * (m_data[n].v + An.z * dt * 0.5);
		rb.m_vt = Q * Vn;

		// angular momentum update
		// NOTE: This currently adds a_{n}, not a_{n+1}, because in order
		//       to evaluate a_{n+1}, I need W_{n+1}. This looks like a nonlinear problem
		//       so probably need to do something else here. 
		vec3d Wn(0,0,0);
		if ((n = rb.m_LM[3]) >= 0) Wn.x = m_dyn_damping * (m_data[n].v + m_data[n].a * dt*0.5);
		if ((n = rb.m_LM[4]) >= 0) Wn.y = m_dyn_damping * (m_data[n].v + m_data[n].a * dt*0.5);
		if ((n = rb.m_LM[5]) >= 0) Wn.z = m_dyn_damping * (m_data[n].v + m_data[n].a * dt*0.5);
		rb.m_wt = Q * Wn;

		mat3ds I0 = rb.m_moi;
		mat3ds I0inv = I0.inverse();
		vec3d Pn = I0inv * (Mn - (Wn ^ (I0 * Wn)));
		rb.m_alt = Q * Pn;

		// update some other stuff
		mat3d R = Q.RotationMatrix();
		mat3d It = R * I0 * R.transpose();
		rb.m_ht = It * rb.m_wt;
	}

	// increase iteration number
	m_niter++;

	// do minor iterations callbacks
	fem.DoCallback(CB_MINOR_ITERS);

	return true;
}

//-----------------------------------------------------------------------------
//! calculates the residual vector
//! Note that the concentrated nodal forces are not calculated here.
//! This is because they do not depend on the geometry 
//! so we only calculate them once (in Quasin) and then add them here.

bool FEExplicitSolidSolver::Residual(vector<double>& R)
{
	TRACK_TIME(Timer_Residual);

	// get the time information
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	const FETimeInfo& tp = fem.GetTime();

	// initialize residual with concentrated nodal loads
	zero(R);

	// zero nodal reaction forces
	zero(m_Fr);

	// setup the global vector
	FEResidualVector RHS(fem, R, m_Fr);

	// zero rigid body reaction forces
	int NRB = fem.RigidBodies();
	for (int i=0; i<NRB; ++i)
	{
		FERigidBody& RB = *fem.GetRigidBody(i);
		RB.m_Fr = RB.m_Mr = vec3d(0,0,0);
	}

	// get the mesh
	FEMesh& mesh = fem.GetMesh();

	// calculate the internal (stress) forces
	for (int i=0; i<mesh.Domains(); ++i)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
		dom.InternalForces(RHS);
	}

	// calculate forces due to model loads
	int nml = fem.ModelLoads();
	for (int i=0; i<nml; ++i)
	{
		FEModelLoad* pml = fem.ModelLoad(i);
		if (pml->IsActive()) pml->LoadVector(RHS);
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

	// set the nodal reaction forces
	// TODO: Is this a good place to do this?
#pragma omp parallel for
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		node.set_load(m_dofU[0], 0);
		node.set_load(m_dofU[1], 0);
		node.set_load(m_dofU[2], 0);
		node.set_load(m_dofSU[0], 0);
		node.set_load(m_dofSU[1], 0);
		node.set_load(m_dofSU[2], 0);

		int n;
		if ((n = -node.m_ID[m_dofU[0]] - 2) >= 0) node.set_load(m_dofU[0], -m_Fr[n]);
		if ((n = -node.m_ID[m_dofU[1]] - 2) >= 0) node.set_load(m_dofU[1], -m_Fr[n]);
		if ((n = -node.m_ID[m_dofU[2]] - 2) >= 0) node.set_load(m_dofU[2], -m_Fr[n]);

		if ((n = -node.m_ID[m_dofSU[0]] - 2) >= 0) node.set_load(m_dofSU[0], -m_Fr[n]);
		if ((n = -node.m_ID[m_dofSU[1]] - 2) >= 0) node.set_load(m_dofSU[1], -m_Fr[n]);
		if ((n = -node.m_ID[m_dofSU[2]] - 2) >= 0) node.set_load(m_dofSU[2], -m_Fr[n]);
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
		if (pci->IsActive()) pci->LoadVector(R, tp);
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
		if (plc->IsActive()) plc->LoadVector(R, tp);
	}
}
