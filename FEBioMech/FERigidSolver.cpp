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
#include "FERigidSolver.h"
#include "FESolidSolver2.h"
#include "FERigidMaterial.h"
#include "FERigidBody.h"
#include <FECore/FEModel.h>
#include "RigidBC.h"
#include <FECore/FEAnalysis.h>
#include <FECore/SparseMatrix.h>
#include <FECore/log.h>
#include <FECore/FEMaterial.h>
#include "FEMechModel.h"
#include <FECore/FELinearSystem.h>
#include "FESolidAnalysis.h"

FERigidSolver::FERigidSolver(FEModel* fem)
{
	m_fem = dynamic_cast<FEMechModel*>(fem);

	m_dofX = m_dofY = m_dofZ = -1;

	m_bAllowMixedBCs = false;
}

int FERigidSolver::InitEquations(int neq)
{
	// Next, we assign equation numbers to the rigid body degrees of freedom
	if (m_fem == nullptr) return neq;

	int nrb = m_fem->RigidBodies();
	for (int i = 0; i<nrb; ++i)
	{
		FERigidBody& RB = *m_fem->GetRigidBody(i);
		for (int j = 0; j<6; ++j)
		{
			int bcj = RB.m_BC[j];
			int lmj = RB.m_LM[j];
			if      (bcj == DOF_OPEN      ) { RB.m_LM[j] = neq; neq++; }
			else if (bcj == DOF_PRESCRIBED) { RB.m_LM[j] = -neq - 2; neq++; }
			else if (bcj == DOF_FIXED     ) { RB.m_LM[j] = -1; }
			else { assert(false); return -1; }
		}
	}

	// get the DOF indices
	m_dofX = m_fem->GetDOFIndex("x");
	m_dofY = m_fem->GetDOFIndex("y");
	m_dofZ = m_fem->GetDOFIndex("z");
	m_dofVX = m_fem->GetDOFIndex("vx");
	m_dofVY = m_fem->GetDOFIndex("vy");
	m_dofVZ = m_fem->GetDOFIndex("vz");
    m_dofU = m_fem->GetDOFIndex("u");
    m_dofV = m_fem->GetDOFIndex("v");
    m_dofW = m_fem->GetDOFIndex("w");
    m_dofSX = m_fem->GetDOFIndex("sx");
    m_dofSY = m_fem->GetDOFIndex("sy");
    m_dofSZ = m_fem->GetDOFIndex("sz");
    m_dofSVX = m_fem->GetDOFIndex("svx");
    m_dofSVY = m_fem->GetDOFIndex("svy");
    m_dofSVZ = m_fem->GetDOFIndex("svz");
	int dofRU = m_fem->GetDOFIndex("Ru");
	int dofRV = m_fem->GetDOFIndex("Rv");
	int dofRW = m_fem->GetDOFIndex("Rw");

	// we assign the rigid body equation number to
	// Also make sure that the nodes are NOT constrained!
	FEMesh& mesh = m_fem->GetMesh();
	for (int i = 0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		if (node.m_rid >= 0)
		{
			FERigidBody& RB = *m_fem->GetRigidBody(node.m_rid);
			node.m_ID[m_dofX] = (RB.m_LM[0] >= 0 ? -RB.m_LM[0] - 2 : RB.m_LM[0]);
			node.m_ID[m_dofY] = (RB.m_LM[1] >= 0 ? -RB.m_LM[1] - 2 : RB.m_LM[1]);
			node.m_ID[m_dofZ] = (RB.m_LM[2] >= 0 ? -RB.m_LM[2] - 2 : RB.m_LM[2]);
			node.m_ID[dofRU] = (RB.m_LM[3] >= 0 ? -RB.m_LM[3] - 2 : RB.m_LM[3]);
			node.m_ID[dofRV] = (RB.m_LM[4] >= 0 ? -RB.m_LM[4] - 2 : RB.m_LM[4]);
			node.m_ID[dofRW] = (RB.m_LM[5] >= 0 ? -RB.m_LM[5] - 2 : RB.m_LM[5]);
			if (node.HasFlags(FENode::SHELL) && node.HasFlags(FENode::RIGID_CLAMP)) {
                node.m_ID[m_dofU] = (RB.m_LM[0] >= 0 ? -RB.m_LM[0] - 2 : RB.m_LM[0]);
                node.m_ID[m_dofV] = (RB.m_LM[1] >= 0 ? -RB.m_LM[1] - 2 : RB.m_LM[1]);
                node.m_ID[m_dofW] = (RB.m_LM[2] >= 0 ? -RB.m_LM[2] - 2 : RB.m_LM[2]);
                node.m_ID[m_dofSX] = (RB.m_LM[0] >= 0 ? -RB.m_LM[0] - 2 : RB.m_LM[0]);
                node.m_ID[m_dofSY] = (RB.m_LM[1] >= 0 ? -RB.m_LM[1] - 2 : RB.m_LM[1]);
                node.m_ID[m_dofSZ] = (RB.m_LM[2] >= 0 ? -RB.m_LM[2] - 2 : RB.m_LM[2]);
            }
		}
	}

	return neq;
}

//-----------------------------------------------------------------------------
//! Serialization
void FERigidSolver::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;
	ar & m_dofX & m_dofY & m_dofZ;
	ar & m_dofVX & m_dofVY & m_dofVZ;
	ar & m_dofU & m_dofV & m_dofW;
    ar & m_dofSX & m_dofSY & m_dofSZ;
    ar & m_dofSVX & m_dofSVY & m_dofSVZ;
	ar & m_bAllowMixedBCs;
}

//-----------------------------------------------------------------------------
//! return the FEModel
FEModel* FERigidSolver::GetFEModel()
{
	return m_fem;
}

//-----------------------------------------------------------------------------
// \todo: eliminate need for ui parameter
void FERigidSolver::PrepStep(const FETimeInfo& timeInfo, vector<double>& ui)
{
	if (m_fem == nullptr) return;
	FEMechModel& fem = *m_fem;

	int NO = fem.RigidBodies();
	for (int i = 0; i<NO; ++i) fem.GetRigidBody(i)->Init();

	// calculate local rigid displacements
	for (int i = 0; i<fem.RigidPrescribedBCs(); ++i)
	{
		FERigidPrescribedBC& DC = *fem.GetRigidPrescribedBC(i);
		if (DC.IsActive()) DC.InitTimeStep();
	}

	// calculate global rigid displacements
	for (int i = 0; i<NO; ++i)
	{
		FERigidBody* prb = fem.GetRigidBody(i);
		if (prb)
		{
			FERigidBody& RB = *prb;
			if (RB.m_prb == 0)
			{
				// if all rotation dofs are fixed or prescribed, set the flag
				if (m_bAllowMixedBCs==false)
				{
					RB.m_bpofr = false;
					if (RB.m_pDC[3] || RB.m_pDC[4] || RB.m_pDC[5])
					{
						bool bpofr[3] = { false };
						for (int j = 3; j<6; ++j) if (RB.m_pDC[j] || (RB.m_LM[j] < 0)) bpofr[j - 3] = true;
						if (bpofr[0] && bpofr[1] && bpofr[2]) RB.m_bpofr = true;
						else
						{
							feLogError("Rigid body rotations cannot mix prescribed and free components.\nRigid body: %d, Material: %d\n", RB.m_nID, RB.GetMaterialID());
							throw "FATAL ERROR";
						}

						// see if all or none prescribed dofs are relative
						bool br[3] = { false, false, false };
						for (int j = 3; j < 6; ++j)
						{
							FERigidPrescribedBC* dc = RB.m_pDC[j];
							if (dc && dc->GetRelativeFlag()) br[j - 3] = true;
						}

						bool bf = ((br[0] == br[1]) && (br[1] == br[2]) && (br[0] == br[2]));
						if (bf==false)
						{
							feLogError("All or none rigid rotations must have relative flag set.\nRigid body: %d, Material: %d\n", RB.m_nID, RB.GetMaterialID());
							throw "FATAL ERROR";
						}
					}
				}

				for (int j = 0; j<6; ++j) RB.m_du[j] = RB.m_dul[j];
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
	for (int i = 0; i<NO; ++i)
	{
		FERigidBody& RB = *fem.GetRigidBody(i);
		for (int j = 0; j<6; ++j)
		{
			int I = -RB.m_LM[j] - 2;
			if (I >= 0) ui[I] = RB.m_du[j];
		}
	}

	FEAnalysis* pstep = m_fem->GetCurrentStep();
	if (pstep->m_nanalysis == FESolidAnalysis::DYNAMIC)
	{
		FEMesh& mesh = m_fem->GetMesh();

		// set the initial velocities of all rigid nodes
		for (int i = 0; i<mesh.Nodes(); ++i)
		{
			FENode& n = mesh.Node(i);
			if (n.m_rid >= 0)
			{
				FERigidBody& rb = *fem.GetRigidBody(n.m_rid);
				vec3d V = rb.m_vt;
				vec3d W = rb.m_wt;
				vec3d r = n.m_rt - rb.m_rt;

				vec3d v = V + (W ^ r);
				n.m_vp = v;
				n.set_vec3d(m_dofVX, m_dofVY, m_dofVZ, v);

				vec3d a = (W ^ V)*2.0 + (W ^ (W ^ r));
				n.m_ap = n.m_at = a;
			}
		}
	}

	// store the current rigid body reaction forces
	for (int i = 0; i<fem.RigidBodies(); ++i)
	{
		FERigidBody& RB = *fem.GetRigidBody(i);
		RB.m_Fp = RB.m_Fr;
		RB.m_Mp = RB.m_Mr;
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the rigid stiffness matrices
void FERigidSolver::RigidStiffness(SparseMatrix& K, vector<double>& ui, vector<double>& F, const FEElementMatrix& ke, double alpha)
{
	if (m_fem == nullptr) return;

	const vector<int>& en = ke.Nodes();
    int n = (int)en.size();
    FEMesh& mesh = m_fem->GetMesh();
    
    bool bclamped_shell = false;
    for (int j = 0; j<n; ++j) 
	{
		if (en[j] >= 0)
		{
			if (mesh.Node(en[j]).HasFlags(FENode::SHELL) && mesh.Node(en[j]).HasFlags(FENode::RIGID_CLAMP)) {
				bclamped_shell = true;
				break;
			}
		}
    }
    if (bclamped_shell)
        RigidStiffnessShell(K, ui, F, en, ke.RowIndices(), ke.ColumnsIndices(), ke, alpha);
    else
        RigidStiffnessSolid(K, ui, F, en, ke.RowIndices(), ke.ColumnsIndices(), ke, alpha);
    return;
}

//-----------------------------------------------------------------------------
//! This function calculates the rigid stiffness matrices
//! correct stiffness matrix for rigid-solid interfaces
void FERigidSolver::RigidStiffnessSolid(SparseMatrix& K, vector<double>& ui, vector<double>& F, const vector<int>& en, const vector<int>& elmi, const std::vector<int>& elmj, const matrix& ke, double alpha)
{
	if (m_fem == nullptr) return;
	FEMechModel& fem = *m_fem;

	if (fem.RigidBodies() == 0) return;
	if (en.empty()) return;
    
    int i, j, k, l, n = (int)en.size();
    
    // get nodal DOFS
    DOFS& fedofs = m_fem->GetDOFS();
    int MAX_NDOFS = fedofs.GetTotalDOFS();

	int ndof = ke.columns() / n;
    
    matrix kij(ndof, ndof);
    matrix KF(ndof, 6);

    double KR[6][6];
    
    int *lmi, *lmj;
    int I, J;
    
    vec3d zi, zj;
    mat3d Zi, Zj;
    
   
    FEMesh& mesh = m_fem->GetMesh();
    
    // loop over columns
    for (j = 0; j<n; ++j)
    {
		if (en[j] >= 0)
		{
			FENode& nodej = mesh.Node(en[j]);
			if (nodej.m_rid >= 0)
			{
				// this is a rigid interface node
				// get the rigid body this node is attached to
				FERigidBody& RBj = *fem.GetRigidBody(nodej.m_rid);

				// get the rigid body equation nrs.
				lmj = RBj.m_LM;

				// get the relative distance to the center of mass
				zj = nodej.m_rt - RBj.m_rt;
				Zj.skew(zj);

				// loop over rows
				for (i = 0; i < n; ++i)
				{
					// get the element sub-matrix
					for (k = 0; k < ndof; ++k)
						for (l = 0; l < ndof; ++l)
							kij[k][l] = ke[ndof*i + k][ndof*j + l];

					mat3d Kuu(kij[0][0], kij[0][1], kij[0][2],
						kij[1][0], kij[1][1], kij[1][2],
						kij[2][0], kij[2][1], kij[2][2]);

					if (en[i] >= 0)
					{
						FENode& nodei = mesh.Node(en[i]);

						if (nodei.m_rid >= 0)
						{
							// node i is also a rigid body node
							// get the rigid body this node is attached to
							FERigidBody& RBi = *fem.GetRigidBody(nodei.m_rid);

							lmi = RBi.m_LM;

							// get the relative distance (use alpha rule)
							zi = (nodei.m_rt - RBi.m_rt)*alpha + (nodei.m_rp - RBi.m_rp)*(1 - alpha);
							Zi.skew(zi);

							mat3d M;

							// Kuu transformation to Krr
							M = Kuu * alpha;
							KR[0][0] = M[0][0]; KR[0][1] = M[0][1]; KR[0][2] = M[0][2];
							KR[1][0] = M[1][0]; KR[1][1] = M[1][1]; KR[1][2] = M[1][2];
							KR[2][0] = M[2][0]; KR[2][1] = M[2][1]; KR[2][2] = M[2][2];


							// Kuu transformation to Krq
							M = Kuu * Zj*(-alpha);
							KR[0][3] = M[0][0]; KR[0][4] = M[0][1]; KR[0][5] = M[0][2];
							KR[1][3] = M[1][0]; KR[1][4] = M[1][1]; KR[1][5] = M[1][2];
							KR[2][3] = M[2][0]; KR[2][4] = M[2][1]; KR[2][5] = M[2][2];


							// Kuu transformation to Kqr
							M = Zi * Kuu*alpha;
							KR[3][0] = M[0][0]; KR[3][1] = M[0][1]; KR[3][2] = M[0][2];
							KR[4][0] = M[1][0]; KR[4][1] = M[1][1]; KR[4][2] = M[1][2];
							KR[5][0] = M[2][0]; KR[5][1] = M[2][1]; KR[5][2] = M[2][2];


							// Kuu transformation to Kqq
							M = Zi * Kuu*Zj*(-alpha);
							KR[3][3] = M[0][0]; KR[3][4] = M[0][1]; KR[3][5] = M[0][2];
							KR[4][3] = M[1][0]; KR[4][4] = M[1][1]; KR[4][5] = M[1][2];
							KR[5][3] = M[2][0]; KR[5][4] = M[2][1]; KR[5][5] = M[2][2];

							// add the stiffness components to the Krr matrix
							for (k = 0; k < 6; ++k)
								for (l = 0; l < 6; ++l)
								{
									J = lmj[k];
									I = lmi[l];

									if (I >= 0)
									{
										// multiply KR by alpha for alpha rule
										if (J < -1) F[I] -= KR[l][k] * ui[-J - 2];
										else if (J >= 0) K.add(I, J, KR[l][k]);
									}
								}

							// we still need to couple the non-rigid degrees of node i to the
							// rigid dofs of node j
							for (k = 3; k < ndof; ++k) {
								vec3d kpu(kij[k][0], kij[k][1], kij[k][2]);
								vec3d m = kpu * alpha;
								KF[k][0] = m.x; KF[k][1] = m.y; KF[k][2] = m.z;
								m = Zj * kpu*alpha;
								KF[k][3] = m.x; KF[k][4] = m.y; KF[k][5] = m.z;
							}

							for (k = 0; k < 6; ++k)
								for (l = 3; l < ndof; ++l)
								{
									J = lmj[k];
									I = elmi[ndof*i + l];

									if (I >= 0)
									{
										// multiply KF by alpha for alpha rule
										if (J < -1) F[I] -= KF[l][k] * ui[-J - 2];
										else if (J >= 0) K.add(I, J, KF[l][k]);
									}
								}

							// now the transpose location
							for (l = 3; l < ndof; ++l) {
								vec3d kup(kij[0][l], kij[1][l], kij[2][l]);
								vec3d m = Zi * kup;
								KF[l][0] = kup.x; KF[l][1] = kup.y; KF[l][2] = kup.z;
								KF[l][3] = m.x; KF[l][4] = m.y; KF[l][5] = m.z;
							}

							for (k = 0; k < 6; ++k)
								for (l = 3; l < ndof; ++l)
								{
									J = elmj[ndof*j + l];
									I = lmi[k];

									if (I >= 0)
									{
										if (J < -1) F[I] -= KF[l][k] * ui[-J - 2];
										else if (J >= 0) K.add(I, J, KF[l][k]);
									}
								}

						}
						else
						{
							// node i is not a rigid body node
							// add the stiffness components to the Kfr matrix

							// Kij
							for (k = 0; k < ndof; ++k) {
								vec3d kpu(kij[k][0], kij[k][1], kij[k][2]);
								vec3d m = kpu * alpha;
								KF[k][0] = m.x; KF[k][1] = m.y; KF[k][2] = m.z;
								m = Zj * kpu*alpha;
								KF[k][3] = m.x; KF[k][4] = m.y; KF[k][5] = m.z;
							}

							for (k = 0; k < 6; ++k)
								for (l = 0; l < ndof; ++l)
								{
									J = lmj[k];
									I = elmi[ndof*i + l];

									if (I >= 0)
									{
										// multiply KF by alpha for alpha rule
										if (J < -1) F[I] -= KF[l][k] * ui[-J - 2];
										else if (J >= 0) K.add(I, J, KF[l][k]);
									}
								}
						}
					}
				}
			}
			else
			{
				// loop over rows
				for (i = 0; i < n; ++i)
				{
					if (en[i] >= 0)
					{
						FENode& nodei = mesh.Node(en[i]);
						if (nodei.m_rid >= 0)
						{
							// node i is a rigid body
							// get the rigid body this node is attached to
							FERigidBody& RBi = *fem.GetRigidBody(nodei.m_rid);

							// get the rigid body equation nrs.
							lmi = RBi.m_LM;

							// get the relative distance (use alpha rule)
							zi = (nodei.m_rt - RBi.m_rt)*alpha + (nodei.m_rp - RBi.m_rp)*(1 - alpha);
							Zi.skew(zi);

							// get the element sub-matrix
							for (k = 0; k < ndof; ++k)
								for (l = 0; l < ndof; ++l)
									kij[k][l] = ke[ndof*i + k][ndof*j + l];

							// add the stiffness components to the Krf matrix

							// Kij
							for (k = 0; k < ndof; ++k) {
								vec3d kup(kij[0][k], kij[1][k], kij[2][k]);
								vec3d m = Zi * kup;
								KF[k][0] = kup.x; KF[k][1] = kup.y; KF[k][2] = kup.z;
								KF[k][3] = m.x; KF[k][4] = m.y; KF[k][5] = m.z;
							}

							for (k = 0; k < 6; ++k)
								for (l = 0; l < ndof; ++l)
								{
									I = lmi[k];
									J = elmj[ndof*j + l];

									if (I >= 0)
									{
										if (J < -1) F[I] -= KF[l][k] * ui[-J - 2];
										else if (J >= 0) K.add(I, J, KF[l][k]);
									}
								}
						}
					}
				}
			}
		}
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the rigid stiffness matrices
//! correct stiffness matrix for rigid bodies accounting for rigid-body-deformable-shell interfaces
void FERigidSolver::RigidStiffnessShell(SparseMatrix& K, vector<double>& ui, vector<double>& F, const vector<int>& en, const vector<int>& elmi, const vector<int>& elmj, const matrix& ke, double alpha)
{
	if (m_fem == nullptr) return;
	FEMechModel& fem = *m_fem;

	if (fem.RigidBodies() == 0) return;
    
    int i, j, k, l, n = (int)en.size();
    
    // get nodal DOFS
    DOFS& fedofs = m_fem->GetDOFS();
    int MAX_NDOFS = fedofs.GetTotalDOFS();
    
    vector< vector<double> > kij; kij.assign(MAX_NDOFS, vector<double>(MAX_NDOFS));
    
    vector< vector<double> > KF; KF.assign(MAX_NDOFS, vector<double>(6));
    double KR[6][6];
    
    int *lmi, *lmj;
    int I, J;
    
    vec3d ai, aj, bi, bj;
    mat3d Ai, Aj, Bi, Bj;
    
    int ndof = ke.columns() / n;
    
    FEMesh& mesh = m_fem->GetMesh();
    
    // loop over columns
    for (j = 0; j<n; ++j)
    {
        FENode& nodej = mesh.Node(en[j]);
        if (nodej.m_rid >= 0)
        {
            // this is a rigid interface node
            // get the rigid body this node is attached to
            FERigidBody& RBj = *fem.GetRigidBody(nodej.m_rid);
            
            // get the rigid body equation nrs.
            lmj = RBj.m_LM;
            
            // get the relative distance to the center of mass
            aj = nodej.m_rt - RBj.m_rt;
            Aj.skew(aj);
            
            // get the shell director
            bj = aj - nodej.m_dt;
            Bj.skew(bj);
            
            // loop over rows
            for (i = 0; i<n; ++i)
            {
                // get the element sub-matrix
                for (k = 0; k<ndof; ++k)
                    for (l = 0; l<ndof; ++l)
                        kij[k][l] = ke[ndof*i + k][ndof*j + l];
                
                mat3d Kuu(kij[0][0], kij[0][1], kij[0][2],
                          kij[1][0], kij[1][1], kij[1][2],
                          kij[2][0], kij[2][1], kij[2][2]);
                
                mat3d Kuw(kij[0][3], kij[0][4], kij[0][5],
                          kij[1][3], kij[1][4], kij[1][5],
                          kij[2][3], kij[2][4], kij[2][5]);
                
                mat3d Kwu(kij[3][0], kij[3][1], kij[3][2],
                          kij[4][0], kij[4][1], kij[4][2],
                          kij[5][0], kij[5][1], kij[5][2]);
                
                mat3d Kww(kij[3][3], kij[3][4], kij[3][5],
                          kij[4][3], kij[4][4], kij[4][5],
                          kij[5][3], kij[5][4], kij[5][5]);
                
                FENode& nodei = mesh.Node(en[i]);
                
                if (nodei.m_rid >= 0)
                {
                    // node i is also a rigid body node
                    // get the rigid body this node is attached to
                    FERigidBody& RBi = *fem.GetRigidBody(nodei.m_rid);
                    
                    lmi = RBi.m_LM;
                    
                    // get the relative distance (use alpha rule)
                    ai = (nodei.m_rt - RBi.m_rt)*alpha + (nodei.m_rp - RBi.m_rp)*(1 - alpha);
                    Ai.skew(ai);
                    
                    // get the shell director
                    bi = ai - nodei.m_dt;
                    Bi.skew(bi);
                    
                    mat3d M;
                    
                    // Kuu transformation
                    M = (Kuu + Kwu + Kuw + Kww)*alpha;
                    KR[0][0] = M[0][0]; KR[0][1] = M[0][1]; KR[0][2] = M[0][2];
                    KR[1][0] = M[1][0]; KR[1][1] = M[1][1]; KR[1][2] = M[1][2];
                    KR[2][0] = M[2][0]; KR[2][1] = M[2][1]; KR[2][2] = M[2][2];
                    
                    
                    // Kuw transformation
                    M = ((Kuu + Kwu)*Aj + (Kuw + Kww)*Bj)*(-alpha);
                    KR[0][3] = M[0][0]; KR[0][4] = M[0][1]; KR[0][5] = M[0][2];
                    KR[1][3] = M[1][0]; KR[1][4] = M[1][1]; KR[1][5] = M[1][2];
                    KR[2][3] = M[2][0]; KR[2][4] = M[2][1]; KR[2][5] = M[2][2];
                    
                    
                    // Kwu transformation
                    M = (Ai*(Kuu + Kuw) + Bi*(Kwu + Kww))*alpha;
                    KR[3][0] = M[0][0]; KR[3][1] = M[0][1]; KR[3][2] = M[0][2];
                    KR[4][0] = M[1][0]; KR[4][1] = M[1][1]; KR[4][2] = M[1][2];
                    KR[5][0] = M[2][0]; KR[5][1] = M[2][1]; KR[5][2] = M[2][2];
                    
                    
                    // Kww transformation
                    M = ((Ai*Kuu + Bi*Kwu)*Aj + (Ai*Kuw + Bi*Kww)*Bj)*(-alpha);
                    KR[3][3] = M[0][0]; KR[3][4] = M[0][1]; KR[3][5] = M[0][2];
                    KR[4][3] = M[1][0]; KR[4][4] = M[1][1]; KR[4][5] = M[1][2];
                    KR[5][3] = M[2][0]; KR[5][4] = M[2][1]; KR[5][5] = M[2][2];
                    
                    // add the stiffness components to the Krr matrix
                    for (k = 0; k<6; ++k)
                        for (l = 0; l<6; ++l)
                        {
                            J = lmj[k];
                            I = lmi[l];
                            
                            if (I >= 0)
                            {
                                // multiply KR by alpha for alpha rule
                                if (J < -1) F[I] -= KR[l][k]*ui[-J - 2];
                                else if (J >= 0) K.add(I, J, KR[l][k]);
                            }
                        }
                    
                    // we still need to couple the non-rigid degrees of node i to the
                    // rigid dofs of node j
                    for (k = 6; k<ndof; ++k) {
                        vec3d kpu(kij[k][0], kij[k][1], kij[k][2]);
                        vec3d kpw(kij[k][3], kij[k][4], kij[k][5]);
                        vec3d m = (kpu + kpw)*alpha;
                        KF[k][0] = m.x; KF[k][1] = m.y; KF[k][2] = m.z;
                        m = (Aj*kpu + Bj*kpw)*alpha;
                        KF[k][3] = m.x; KF[k][4] = m.y; KF[k][5] = m.z;
                    }
                    
                    for (k = 0; k<6; ++k)
                        for (l = 6; l<ndof; ++l)
                        {
                            J = lmj[k];
                            I = elmi[ndof*i + l];
                            
                            if (I >= 0)
                            {
                                // multiply KF by alpha for alpha rule
                                if (J < -1) F[I] -= KF[l][k] * ui[-J - 2];
                                else if (J >= 0) K.add(I, J, KF[l][k]);
                            }
                        }
                    
                    // now the transpose location
                    for (l = 6; l<ndof; ++l) {
                        vec3d kup(kij[0][l], kij[1][l], kij[2][l]);
                        vec3d kwp(kij[3][l], kij[4][l], kij[5][l]);
                        vec3d m = kup + kwp;
                        KF[l][0] = m.x; KF[l][1] = m.y; KF[l][2] = m.z;
                        m = Ai*kup + Bi*kwp;
                        KF[l][3] = m.x; KF[l][4] = m.y; KF[l][5] = m.z;
                    }
                    
                    for (k = 0; k<6; ++k)
                        for (l = 6; l<ndof; ++l)
                        {
                            J = elmj[ndof*j + l];
                            I = lmi[k];
                            
                            if (I >= 0)
                            {
                                if (J < -1) F[I] -= KF[l][k] * ui[-J - 2];
                                else if (J >= 0) K.add(I, J, KF[l][k]);
                            }
                        }
                    
                }
                else
                {
                    // node i is not a rigid body node
                    // add the stiffness components to the Kfr matrix
                    
                    // Kij
                    for (k = 0; k<ndof; ++k) {
                        vec3d kpu(kij[k][0], kij[k][1], kij[k][2]);
                        vec3d kpw(kij[k][3], kij[k][4], kij[k][5]);
                        vec3d m = (kpu + kpw)*alpha;
                        KF[k][0] = m.x; KF[k][1] = m.y; KF[k][2] = m.z;
                        m = (Aj*kpu + Bj*kpw)*alpha;
                        KF[k][3] = m.x; KF[k][4] = m.y; KF[k][5] = m.z;
                    }
                    
                    for (k = 0; k<6; ++k)
                        for (l = 0; l<ndof; ++l)
                        {
                            J = lmj[k];
                            I = elmi[ndof*i + l];
                            
                            if (I >= 0)
                            {
                                // multiply KF by alpha for alpha rule
                                if (J < -1) F[I] -= KF[l][k] * ui[-J - 2];
                                else if (J >= 0) K.add(I, J, KF[l][k]);
                            }
                        }
                }
            }
        }
        else
        {
            // loop over rows
            for (i = 0; i<n; ++i)
            {
                FENode& nodei = mesh.Node(en[i]);
                if (nodei.m_rid >= 0)
                {
                    // node i is a rigid body
                    // get the rigid body this node is attached to
                    FERigidBody& RBi = *fem.GetRigidBody(nodei.m_rid);
                    
                    // get the rigid body equation nrs.
                    lmi = RBi.m_LM;
                    
                    // get the relative distance (use alpha rule)
                    ai = (nodei.m_rt - RBi.m_rt)*alpha + (nodei.m_rp - RBi.m_rp)*(1 - alpha);
                    Ai.skew(ai);
                    
                    // get the shell director
                    bi = ai - nodei.m_dt;
                    Bi.skew(bi);
                    
                    // get the element sub-matrix
                    for (k = 0; k<ndof; ++k)
                        for (l = 0; l<ndof; ++l)
                            kij[k][l] = ke[ndof*i + k][ndof*j + l];
                    
                    // add the stiffness components to the Krf matrix
                    
                    // Kij
                    for (k = 0; k<ndof; ++k) {
                        vec3d kup(kij[0][k], kij[1][k], kij[2][k]);
                        vec3d kwp(kij[3][k], kij[4][k], kij[5][k]);
                        vec3d m = kup + kwp;
                        KF[k][0] = m.x; KF[k][1] = m.y; KF[k][2] = m.z;
                        m = Ai*kup + Bi*kwp;
                        KF[k][3] = m.x; KF[k][4] = m.y; KF[k][5] = m.z;
                    }
                    
                    for (k = 0; k<6; ++k)
                        for (l = 0; l<ndof; ++l)
                        {
                            I = lmi[k];
                            J = elmj[ndof*j + l];
                            
                            if (I >= 0)
                            {
                                if (J < -1) F[I] -= KF[l][k] * ui[-J - 2];
                                else if (J >= 0) K.add(I, J, KF[l][k]);
                            }
                        }
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FERigidSolver::AssembleResidual(int node_id, int dof, double f, vector<double>& R)
{
	if (m_fem == nullptr) return;
	FEMechModel& fem = *m_fem;

	FEMesh& mesh = m_fem->GetMesh();
	
    // get the equation number
    FENode& node = mesh.Node(node_id);
    int n = node.m_ID[dof];
    
    // assemble into global vector
    if (n >= 0) R[n] += f;
    else if (node.m_rid >= 0)
    {
        // this is a rigid body node
        FERigidBody& RB = *fem.GetRigidBody(node.m_rid);
        
        // get the relative position
        vec3d a = node.m_rt - RB.m_rt;
        
        int* lm = RB.m_LM;
        if (dof == m_dofX)
        {
            if (lm[0] >= 0) R[lm[0]] += f;
            if (lm[4] >= 0) R[lm[4]] += a.z*f;
            if (lm[5] >= 0) R[lm[5]] += -a.y*f;
        }
        else if (dof == m_dofY)
        {
            if (lm[1] >= 0) R[lm[1]] += f;
            if (lm[3] >= 0) R[lm[3]] += -a.z*f;
            if (lm[5] >= 0) R[lm[5]] += a.x*f;
        }
        else if (dof == m_dofZ)
        {
            if (lm[2] >= 0) R[lm[2]] += f;
            if (lm[3] >= 0) R[lm[3]] += a.y*f;
            if (lm[4] >= 0) R[lm[4]] += -a.x*f;
        }
		if (node.HasFlags(FENode::SHELL) && node.HasFlags(FENode::RIGID_CLAMP)) {
            // get the shell director
            vec3d d = node.m_dt;
            vec3d b = a - d;
            if (dof == m_dofSX)
            {
                if (lm[0] >= 0) R[lm[0]] +=  f;
                if (lm[4] >= 0) R[lm[4]] += b.z*f;
                if (lm[5] >= 0) R[lm[5]] += -b.y*f;
            }
            else if (dof == m_dofSY)
            {
                if (lm[1] >= 0) R[lm[1]] +=  f;
                if (lm[3] >= 0) R[lm[3]] += -b.z*f;
                if (lm[5] >= 0) R[lm[5]] += b.x*f;
            }
            else if (dof == m_dofSZ)
            {
                if (lm[2] >= 0) R[lm[2]] +=  f;
                if (lm[3] >= 0) R[lm[3]] += b.y*f;
                if (lm[4] >= 0) R[lm[4]] += -b.x*f;
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FERigidSolver::Residual()
{
	if (m_fem == nullptr) return;
	FEMechModel& fem = *m_fem;

	int NRB = fem.RigidBodies();
	for (int i = 0; i<NRB; ++i)
	{
		FERigidBody& RB = *fem.GetRigidBody(i);
		RB.m_Fr = RB.m_Mr = vec3d(0, 0, 0);
	}
}

//-----------------------------------------------------------------------------
void FERigidSolver::StiffnessMatrix(SparseMatrix& K, const FETimeInfo& tp)
{
	if (m_fem == nullptr) return;
	FEMechModel& fem = *m_fem;

	// we still need to set the diagonal elements to 1
	// for the prescribed rigid body dofs.
	int NRB = fem.RigidBodies();
	for (int i = 0; i<NRB; ++i)
	{
		FERigidBody& rb = *fem.GetRigidBody(i);
		for (int j = 0; j<6; ++j)
		if (rb.m_LM[j] < -1)
		{
			int I = -rb.m_LM[j] - 2;
			K.set(I, I, 1);
		}
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the contribution to the mass matrix from the rigid bodies
void FERigidSolver::RigidMassMatrix(FELinearSystem& LS, const FETimeInfo& timeInfo)
{
	if (m_fem == nullptr) return;
	FEMechModel& fem = *m_fem;

	// element stiffness matrix
	vector<int> lm;
	FEElementMatrix ke;

	// 6 dofs per rigid body
	ke.resize(6, 6);

	// Newmark integration rule
	double dt = timeInfo.timeIncrement;
    double alpham = timeInfo.alpham;
	double beta = timeInfo.beta;
	double gamma = timeInfo.gamma;
	double a = 1. / (beta*dt*dt);

	for (int i=0; i<fem.RigidBodies(); ++i)
	{
		FERigidBody& RB = *fem.GetRigidBody(i);

		// mass matrix
		double M = RB.m_mass*a*alpham;

		ke.zero();
		ke[0][0] = M;
		ke[1][1] = M;
		ke[2][2] = M;

		// evaluate mass moment of inertia at t
		mat3d Rt = RB.GetRotation().RotationMatrix();
		mat3ds Jt = (Rt*RB.m_moi*Rt.transpose()).sym();

		// incremental rotation in spatial frame
		quatd q = RB.GetRotation()*RB.m_qp.Inverse();
		q.MakeUnit();                           // clean-up roundoff errors
		double theta = 2 * tan(q.GetAngle() / 2);   // get theta from Cayley transform
		vec3d e = q.GetVector();

		// skew-symmetric tensor whose axial vector is the incremental rotation
		mat3d qhat;
		qhat.skew(e*theta);

		// generate tensor T(theta)
		mat3d T = mat3dd(1) + qhat / 2 + dyad(e*theta) / 4;

		// skew-symmetric of angular momentum
		mat3d hhat;
		hhat.skew(Jt*RB.m_wt);

		// rotational inertia stiffness
		mat3d K = ((Jt*T)/(beta*dt) - hhat/gamma)*(alpham/dt);

		ke[3][3] = K(0, 0); ke[3][4] = K(0, 1); ke[3][5] = K(0, 2);
		ke[4][3] = K(1, 0); ke[4][4] = K(1, 1); ke[4][5] = K(1, 2);
		ke[5][3] = K(2, 0); ke[5][4] = K(2, 1); ke[5][5] = K(2, 2);

		lm.assign(RB.m_LM, RB.m_LM + 6);

		ke.SetIndices(lm);
		LS.Assemble(ke);
	}
}

//=================================================================================================
// FERigidSolverOld
//=================================================================================================

//-----------------------------------------------------------------------------
//! This function updates the rigid body linear and angular velocity by solving
//! an overdetermined system of linear equations using the least-square method. 
void FERigidSolverOld::UpdateRigidKinematics()
{
	// get the model and mesh
	if (m_fem == nullptr) return;
	FEMechModel& fem = *m_fem;

	FEMesh& mesh = fem.GetMesh();

	// loop over all rigid bodies
	int NRB = fem.RigidBodies();
	for (int j = 0; j<NRB; ++j)
	{
		// get the rigid body
		FERigidBody& rb = *fem.GetRigidBody(j);

		// right-hand side and least-square matrix
		vector<double> r; r.assign(6, 0.0);
		matrix m(6, 6); m.zero();

		// we need to loop over all domains that define this rigid body
		int ncnt = 0;
		int NDOM = mesh.Domains();
		for (int n = 0; n<NDOM; ++n)
		{
			FEDomain& dom = mesh.Domain(n);
			FERigidMaterial* prm = dynamic_cast<FERigidMaterial*>(dom.GetMaterial());
			if (prm)
			{
				if (prm->GetRigidBodyID() == j)
				{
					// now loop over all the nodes
					int NN = dom.Nodes();
					for (int i = 0; i<NN; ++i, ncnt++)
					{
						vec3d ri = dom.Node(i).m_rt - rb.m_rt;
						vec3d vi = dom.Node(i).get_vec3d(m_dofVX, m_dofVY, m_dofVZ);

						vec3d wi = ri ^ vi;

						// right-hand side
						r[0] += vi.x;
						r[1] += vi.y;
						r[2] += vi.z;
						r[3] += wi.x;
						r[4] += wi.y;
						r[5] += wi.z;

						// least-squares matrix
						m[0][0] += 1.0;
						m[1][1] += 1.0;
						m[2][2] += 1.0;

						m[0][4] += ri.z; m[0][5] += -ri.y;
						m[1][3] += -ri.z; m[1][5] += ri.x;
						m[2][3] += ri.y; m[2][4] += -ri.x;

						m[3][4] += -ri.z; m[3][5] += ri.y;
						m[4][3] += ri.z; m[4][5] += -ri.x;
						m[5][3] += -ri.y; m[5][4] += ri.x;

						m[3][3] += ri.y*ri.y + ri.z*ri.z; m[3][4] += -ri.x*ri.y; m[3][5] += -ri.x*ri.z;
						m[4][4] += ri.x*ri.x + ri.z*ri.z; m[4][3] += -ri.x*ri.y; m[4][5] += -ri.y*ri.z;
						m[5][5] += ri.x*ri.x + ri.y*ri.y; m[5][3] += -ri.x*ri.z; m[5][4] += -ri.y*ri.z;
					}
				}
			}
		}

		// solve for the rigid body velocity (if we have enough nodes)
		if (ncnt > 2)
		{
			vector<double> VR = r / m;
			rb.m_vt = vec3d(VR[0], VR[1], VR[2]);
			rb.m_wt = vec3d(VR[3], VR[4], VR[5]);
		}
	}
}

//-----------------------------------------------------------------------------
//! Updates the rigid body data
void FERigidSolverOld::UpdateRigidBodies(vector<double>& Ui, vector<double>& ui, bool bnewUpdate)
{
	if (m_fem == nullptr) return;
	FEMechModel& fem = *m_fem;

	const int NRB = fem.RigidBodies();

	// first calculate the rigid body displacement increments
	for (int i = 0; i<NRB; ++i)
	{
		// get the rigid body
		FERigidBody& RB = *fem.GetRigidBody(i);
		int *lm = RB.m_LM;
		double* du = RB.m_du;

		if (RB.m_prb == 0)
		{
			for (int j = 0; j<6; ++j)
			{
				du[j] = (lm[j] >= 0 ? Ui[lm[j]] + ui[lm[j]] : 0);
			}
		}
	}

	// for prescribed displacements, the displacement increments are evaluated differently
	// TODO: Is this really necessary? Why can't the ui vector contain the correct values?
	const int NRD = fem.RigidPrescribedBCs();
	for (int i = 0; i<NRD; ++i)
	{
		FERigidPrescribedBC& dc = *fem.GetRigidPrescribedBC(i);
		if (dc.IsActive())
		{
			FERigidBody& RB = *fem.GetRigidBody(dc.GetID());
			if (RB.m_prb == 0)
			{
				int I = dc.GetBC();
				RB.m_du[I] = dc.Value() - RB.m_Up[I];
			}
		}
	}

	// update the rigid bodies
	for (int i = 0; i<NRB; ++i)
	{
		// get the rigid body
		FERigidBody& RB = *fem.GetRigidBody(i);
		double* du = RB.m_du;

		if (bnewUpdate)
		{
			// This is the "new" update algorithm which addressesses a couple issues
			// with the old method, namely that prescribed rotational dofs aren't update correctly.
			// Unfortunately, it seems to produce worse convergence in some cases, especially with line search
			// and it doesn't work when rigid bodies are used in a hierarchy
			if (RB.m_prb) du = RB.m_dul;
			RB.m_Ut[0] = RB.m_Up[0] + du[0];
			RB.m_Ut[1] = RB.m_Up[1] + du[1];
			RB.m_Ut[2] = RB.m_Up[2] + du[2];
			RB.m_Ut[3] = RB.m_Up[3] + du[3];
			RB.m_Ut[4] = RB.m_Up[4] + du[4];
			RB.m_Ut[5] = RB.m_Up[5] + du[5];

			RB.m_rt = RB.m_r0 + vec3d(RB.m_Ut[0], RB.m_Ut[1], RB.m_Ut[2]);

			vec3d Rt(RB.m_Ut[3], RB.m_Ut[4], RB.m_Ut[5]);
			RB.SetRotation(quatd(Rt));
		}
		else
		{
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
	}

	// we need to update the position of rigid nodes
	fem.UpdateRigidMesh();

	// Since the rigid nodes are repositioned we need to update the displacement DOFS
	FEMesh& mesh = m_fem->GetMesh();
	int N = mesh.Nodes();
	for (int i = 0; i<N; ++i)
	{
		FENode& node = mesh.Node(i);
		if (node.m_rid >= 0)
		{
			vec3d ut = node.m_rt - node.m_r0;
			node.set_vec3d(m_dofX, m_dofY, m_dofZ, ut);
		}
	}
}

//=================================================================================================
// FERigidSolverNew
//=================================================================================================

//-----------------------------------------------------------------------------
void FERigidSolverNew::UpdateIncrements(vector<double>& Ui, vector<double>& ui, bool emap)
{
	if (m_fem == nullptr) return;
	FEMechModel& fem = *m_fem;

	int nrb = fem.RigidBodies();
	for (int i = 0; i<nrb; ++i)
	{
		// get the rigid body
		FERigidBody& RB = *fem.GetRigidBody(i);
		if (RB.m_prb == 0)
		{
			int *lm = RB.m_LM;
			vec3d v;
			quatd qUi;

			// first do the displacements
			for (int j = 0; j<3; ++j)
			if (lm[j] >= 0) Ui[lm[j]] += ui[lm[j]];

			// next, we do the rotations. We do this separately since
			// they need to be interpreted differently than displacements
			vec3d vUi(0, 0, 0);
			vec3d vui(0, 0, 0);
			if (lm[3] >= 0) { vUi.x = Ui[lm[3]]; vui.x = ui[lm[3]]; }
			if (lm[4] >= 0) { vUi.y = Ui[lm[4]]; vui.y = ui[lm[4]]; }
			if (lm[5] >= 0) { vUi.z = Ui[lm[5]]; vui.z = ui[lm[5]]; }
			if (emap) qUi = quatd(vUi);
			else qUi = quatd(2 * atan(vUi.norm() / 2), vUi);     // Cayley transform
			quatd qui(2 * atan(vui.norm() / 2), vui);            // Cayley transform
			quatd q = qui*qUi;
			q.MakeUnit();
			if (emap) v = q.GetVector()*q.GetAngle();
			else v = q.GetVector()*(2 * tan(q.GetAngle() / 2)); // Cayley transform
			if (lm[3] >= 0) Ui[lm[3]] = v.x;
			if (lm[4] >= 0) Ui[lm[4]] = v.y;
			if (lm[5] >= 0) Ui[lm[5]] = v.z;
		}
	}
}


//-----------------------------------------------------------------------------
//! Updates the rigid body data
void FERigidSolverNew::UpdateRigidBodies(vector<double>& Ui, vector<double>& ui)
{
	// update rigid bodies
	if (m_fem == nullptr) return;
	FEMechModel& fem = *m_fem;

	int nrb = fem.RigidBodies();
	for (int i = 0; i<nrb; ++i)
	{
		// get the rigid body
		FERigidBody& RB = *fem.GetRigidBody(i);
		int *lm = RB.m_LM;
		double* du = RB.m_du;

		// first do the displacements
		if (RB.m_prb == 0)
		{
			for (int j = 0; j<3; ++j)
			{
				FERigidPrescribedBC* pdc = RB.m_pDC[j];
				if (pdc)
				{
					// TODO: do I need to take the line search step into account here?
					du[j] = pdc->Value() - RB.m_Up[j];
				}
				else
				{
					du[j] = (lm[j] >= 0 ? Ui[lm[j]] + ui[lm[j]] : 0);
				}
			}
		}

		RB.m_rt.x = RB.m_rp.x + du[0];
		RB.m_rt.y = RB.m_rp.y + du[1];
		RB.m_rt.z = RB.m_rp.z + du[2];

		// update RB center of mass translations
		if (RB.m_prb) du = RB.m_dul;
		RB.m_Ut[0] = RB.m_Up[0] + du[0];
		RB.m_Ut[1] = RB.m_Up[1] + du[1];
		RB.m_Ut[2] = RB.m_Up[2] + du[2];

		// next, we do the rotations. We do this seperatly since
		// they need to be interpreted differently than displacements
		if (RB.m_prb == 0)
		{
			if (RB.m_bpofr) {
				// if all rotation components are known (prescribed or fixed)
				// evaluate net increment from load curve
				double Ut[3] = { 0 };
				for (int j = 3; j<6; ++j) {
					if (RB.m_pDC[j]) {
						Ut[j - 3] = RB.m_pDC[j]->Value();
					}
				}
				quatd qUt(vec3d(Ut[0], Ut[1], Ut[2]));
				RB.SetRotation(qUt);
				RB.m_Ut[3] = Ut[0];
				RB.m_Ut[4] = Ut[1];
				RB.m_Ut[5] = Ut[2];
			}
			else
			{
				// rotation components are either free or fixed
				vec3d vUi(0, 0, 0);   // initialize total increment so far
				vec3d vui(0, 0, 0);   // initialize current increment
				if (lm[3] >= 0) { vUi.x = Ui[lm[3]]; vui.x = ui[lm[3]]; }
				if (lm[4] >= 0) { vUi.y = Ui[lm[4]]; vui.y = ui[lm[4]]; }
				if (lm[5] >= 0) { vUi.z = Ui[lm[5]]; vui.z = ui[lm[5]]; }
				quatd qUi(2 * atan(vUi.norm() / 2), vUi);                   // Cayley transform
				quatd qui(2 * atan(vui.norm() / 2), vui);                   // Cayley transform
				quatd qdu = qui*qUi;                                        // quaternion of net increment

				qdu.MakeUnit();                                         // clean-up roundoff errors
				quatd Q = qdu*RB.m_qp;
				Q.MakeUnit();

				// update at the current time step
				RB.SetRotation(Q);

				// update RB rotations
				vec3d vUt = Q.GetRotationVector();
				RB.m_Ut[3] = vUt.x;
				RB.m_Ut[4] = vUt.y;
				RB.m_Ut[5] = vUt.z;
			}
		}
	}

	// update the mesh' nodes
	fem.UpdateRigidMesh();

	// Since the rigid nodes are repositioned we need to update the displacement DOFS
	FEMesh& mesh = m_fem->GetMesh();
	int N = mesh.Nodes();
	for (int i = 0; i<N; ++i)
	{
		FENode& node = mesh.Node(i);
		if (node.m_rid >= 0)
		{
			vec3d ut = node.m_rt - node.m_r0;
			node.set_vec3d(m_dofX, m_dofY, m_dofZ, ut);
            
			if (node.HasFlags(FENode::SHELL) && node.HasFlags(FENode::RIGID_CLAMP)) {
                // get the rigid body
                FERigidBody& RB = *fem.GetRigidBody(node.m_rid);
                // evaluate the director in the current configuration
				vec3d d = RB.GetRotation()*node.m_d0 - node.m_d0;
                // evaluate the back face displacement increments
                node.set_vec3d(m_dofSX, m_dofSY, m_dofSZ, ut - d);
            }
		}
	}
}

//-----------------------------------------------------------------------------
//! evaluate inertia forces
void FERigidSolverNew::InertialForces(FEGlobalVector& R, const FETimeInfo& timeInfo)
{
	// Newmark rule
    double alpham = timeInfo.alpham;
    double beta = timeInfo.beta;
    double gamma = timeInfo.gamma;
	double dt = timeInfo.timeIncrement;
	double a = 1.0 / (beta*dt);
	double b = a / dt;
	double c = 1.0 - 0.5 / beta;

	if (m_fem == nullptr) return;
	FEMechModel& fem = *m_fem;

	int nrb = fem.RigidBodies();
	for (int i = 0; i<nrb; ++i)
	{
		// get the rigid body
		FERigidBody& RB = *fem.GetRigidBody(i);

		// acceleration and velocity of center of mass
		RB.m_at = (RB.m_rt - RB.m_rp)*b - RB.m_vp*a + RB.m_ap*c;
		RB.m_vt = RB.m_vp + (RB.m_ap*(1.0 - gamma) + RB.m_at*gamma)*dt;
		// angular acceleration and velocity of rigid body
		quatd q = RB.GetRotation()*RB.m_qp.Inverse();
		q.MakeUnit();
		vec3d vq = q.GetVector()*(2 * tan(q.GetAngle() / 2));  // Cayley transform
		RB.m_wt = vq*(a*gamma) - RB.m_wp + (RB.m_wp + RB.m_alp*dt / 2.)*(2 - gamma / beta);
		q.RotateVector(RB.m_wt);
		RB.m_alt = vq*b - RB.m_wp*a + RB.m_alp*c;
		q.RotateVector(RB.m_alt);
	}

	// calculate rigid body inertial forces
	for (int i = 0; i<nrb; ++i)
	{
		FERigidBody& RB = *fem.GetRigidBody(i);

		// 6 dofs per rigid body
		vector<double> fe(6);
		vector<int>	LM(6);

		// rate of change of linear momentum = mass*acceleration
		vec3d F = (RB.m_at*alpham + RB.m_ap*(1-alpham))*RB.m_mass;

		fe[0] = -F.x;
		fe[1] = -F.y;
		fe[2] = -F.z;

		// evaluate mass moment of inertia at t and tp
		mat3d Rt = RB.GetRotation().RotationMatrix();
		mat3ds Jt = (Rt*RB.m_moi*Rt.transpose()).sym();
        RB.m_ht = Jt*RB.m_wt;

		// evaluate rate of change of angular momentum
		RB.m_dht = (RB.m_ht - RB.m_hp) / (gamma*dt) + RB.m_dhp*(1-1.0/gamma);
//        RB.m_dht = (RB.m_wt ^ RB.m_ht) + Jt*RB.m_alt;

        vec3d M = (RB.m_dht*alpham + RB.m_dhp*(1-alpham));

		fe[3] = -M.x;
		fe[4] = -M.y;
		fe[5] = -M.z;

		// pack the equation numbers
		LM[0] = RB.m_LM[0];
		LM[1] = RB.m_LM[1];
		LM[2] = RB.m_LM[2];
		LM[3] = RB.m_LM[3];
		LM[4] = RB.m_LM[4];
		LM[5] = RB.m_LM[5];
		R.Assemble(LM, fe);
	}
}
