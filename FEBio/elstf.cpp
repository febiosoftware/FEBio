#include "stdafx.h"
#include <math.h>
#include "FESolver.h"
#include "FEPoroElastic.h"
#include "tens4d.h"

//-----------------------------------------------------------------------------
//! Calculates global stiffness matrix.

bool FESolver::StiffnessMatrix()
{
	int iel;

	// get the stiffness matrix
	SparseMatrix& K = *m_pK;

	// zero stiffness matrix
	K.zero();

	// zero the residual adjustment vector
	m_Fd.zero();

	// element stiffness matrix
	matrix ke;

	// nodal degrees of freedom
	int i, j, I, ndof;

	// get the mesh
	FEMesh& mesh = m_fem.m_mesh;

	// get the number of solid elements
	int NE = mesh.SolidElements();

	// get the number of shell elements
	int NS = mesh.ShellElements();

	// repeat over all solid elements
	for (iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = mesh.SolidElement(iel);
		if (!el.IsRigid())
		{
			mesh.UnpackElement(el);

			// get the elements material
			FEMaterial* pmat = m_fem.GetMaterial(el.GetMatID());

			// skip rigid elements and poro-elastic elements
			if (dynamic_cast<FEPoroElastic*>(pmat) == 0)
			{
				// create the element's stiffness matrix
				ndof = 3*el.Nodes();
				ke.Create(ndof, ndof);
				ke.zero();

				// calculate the element stiffness matrix
				ElementStiffness(el, ke);

				// add the inertial stiffness for dynamics
				if (m_fem.m_pStep->m_itype == FE_DYNAMIC) ElementInertialStiffness(el, ke);

				// assemble element matrix in global stiffness matrix
				AssembleStiffness(el.m_node, el.LM(), ke);
			}
			else if (dynamic_cast<FEPoroElastic*>(pmat))
			{
				// allocate stiffness matrix
				ndof = el.Nodes()*4;
				ke.Create(ndof, ndof);
		
				// calculate the element stiffness matrix
				ElementPoroStiffness(el, ke);
	
				// assemble element matrix in global stiffness matrix
				AssembleStiffness(el.m_node, el.LM(), ke);
			}
		}
		else
		{
			// for dynamic analyses we do need to add the inertial stiffness of the rigid body
			if (m_fem.m_pStep->m_itype == FE_DYNAMIC)
			{
				mesh.UnpackElement(el);

				ndof = 3*el.Nodes();
				ke.Create(ndof, ndof);
				ke.zero();

				// add the inertial stiffness for dynamics
				ElementInertialStiffness(el, ke);

				// assemble element matrix in global stiffness matrix
				AssembleStiffness(el.m_node, el.LM(), ke);
			}
		}

		if (m_fem.m_pStep->GetPrintLevel() == FE_PRINT_MINOR_ITRS_EXP)
		{
			fprintf(stderr, "Calculating stiffness matrix: %.1lf %% \r", 100.0*iel/(NE + NS));
		}
	}

	// repeat over all shell elements
	for (iel=0; iel<NS; ++iel)
	{
		FEShellElement& el = mesh.ShellElement(iel);
		if (!el.IsRigid())
		{
			mesh.UnpackElement(el);

			// get the elements material
			FEMaterial* pmat = m_fem.GetMaterial(el.GetMatID());

			// skip rigid elements and poro-elastic elements
			if (dynamic_cast<FEPoroElastic*>(pmat) == 0)
			{
				// create the element's stiffness matrix
				ndof = 6*el.Nodes();
				ke.Create(ndof, ndof);

				// calculate the element stiffness matrix
				ElementStiffness(el, ke);

				// assemble element matrix in global stiffness matrix
				AssembleStiffness(el.m_node, el.LM(), ke);
			}
			else if (dynamic_cast<FEPoroElastic*>(pmat))
			{
				// TODO: implement poro-elasticity for shells
			}		
		}

		if (m_fem.m_pStep->GetPrintLevel() == FE_PRINT_MINOR_ITRS_EXP)
		{
			fprintf(stderr, "Calculating stiffness matrix: %.1lf %% \r", 100.0*(NE + iel)/(NE + NS));
		}
	}


	// calculate contact stiffness
	if (m_fem.m_bcontact) 
	{
		ContactStiffness();
	}

	// calculate joint stiffness 
	if (m_fem.m_nrj)
	{
		for (int i=0; i<m_fem.m_nrj; ++i) m_fem.m_RJ[i].JointStiffness();
	}

	// calculate pressure stiffness term
	if (m_fem.m_pStep->m_istiffpr != 0) 
	{
		int npr = m_fem.m_PC.size();
		for (int m=0; m<npr; ++m)
		{
			// get the surface element
			FESurfaceElement& el = m_fem.m_psurf->Element(m);

			// skip rigid surface elements
			// TODO: do we really need to skip rigid elements?
			if (!el.IsRigid())
			{
				mesh.UnpackElement(el);

				// calculate nodal pressures
				double* pt = el.pt();
				FEPressureLoad& pc = m_fem.m_PC[m];

				if (!pc.blinear)
				{
					double g = m_fem.GetLoadCurve(pc.lc)->Value();

					for (j=0; j<el.Nodes(); ++j) pt[j] = -g*pc.s[j];

					// get the element stiffness matrix
					ndof = 3*el.Nodes();
					ke.Create(ndof, ndof);

					// calculate pressure stiffness
					PressureStiffness(el, ke);

					// assemble element matrix in global stiffness matrix
					AssembleStiffness(el.m_node, el.LM(), ke);
				}
			}
		}
	}

	// discrete element stiffness
	if (m_fem.m_DE.size())
	{
		DiscreteElementStiffness();
	}

	// calculate linear constraint stiffness
	// note that this is the contribution of the 
	// constrainst enforced with augmented lagrangian
	LinearConstraintStiffness();

	// we still need to set the diagonal elements to 1
	// for the prescribed rigid body dofs.
	for (i=0; i<m_fem.m_nrb; ++i)
	{
		FERigidBody& rb = m_fem.m_RB[i];
		for (j=0; j<6; ++j)
			if (rb.m_LM[j] < -1)
			{
				I = -rb.m_LM[j]-2;
				K.set(I,I, 1);
			}
	}

	// let's check the stiffness matrix for zero diagonal elements
	if (m_fem.GetDebugFlag())
	{
		for (i=0; i<m_fem.m_neq; ++i)
		{
			if (K.diag(i) == 0) throw ZeroDiagonal(i, m_fem);
		}
	}

	return true;
}

//-----------------------------------------------------------------------------

void FESolver::LinearConstraintStiffness()
{
	int N = m_fem.m_LCSet.size();
	if (N > 0)
	{
		list<FELinearConstraintSet*>::iterator im = m_fem.m_LCSet.begin();
		for (int i=0; i<N; ++i, ++im) (*im)->Stiffness();
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the contact stiffness matrix

void FESolver::ContactStiffness()
{
	for (int i=0; i<m_fem.ContactInterfaces(); ++i) m_fem.m_CI[i].ContactStiffness();
}

//-----------------------------------------------------------------------------
//! This function calculates the rigid stiffness matrices

void FESolver::RigidStiffness(vector<int>& en, vector<int>& elm, matrix& ke)
{
	int i, j, k, l, n = en.size();
	double kij[MAX_NDOFS][MAX_NDOFS], Ri[3][3] = {0}, Rj[3][3] = {0};

	double KF[MAX_NDOFS][6];
	double KR[6][6];

	int *lmi, *lmj;
	int I, J;

	SparseMatrix& K = *m_pK;

	vec3d ai, aj;

	int ndof = ke.columns() / n;

	// loop over columns
	for (j=0; j<n; ++j)
	{
		FENode& nodej = m_fem.m_mesh.Node(en[j]);
		if (nodej.m_rid >= 0)
		{
			// this is a rigid interface node
			// get the rigid body this node is attached to
			FERigidBody& RBj = m_fem.m_RB[nodej.m_rid];

			// get the rigid body equation nrs.
			lmj = RBj.m_LM;

			// get the relative distance to the center of mass
			aj = nodej.m_rt - RBj.m_rt;
	
			Rj[0][1] = aj.z; Rj[0][2] =-aj.y;
			Rj[1][0] =-aj.z; Rj[1][2] = aj.x;
			Rj[2][0] = aj.y; Rj[2][1] =-aj.x;

			// loop over rows
			for (i=0; i<n; ++i)
			{
				// get the element sub-matrix
				for (k=0; k<ndof; ++k)
					for (l=0; l<ndof; ++l)
						kij[k][l] = ke[ndof*i+k][ndof*j+l];

				FENode& nodei = m_fem.m_mesh.Node(en[i]);

				if (nodei.m_rid>=0)
				{
					// node i is also a rigid body node
					// get the rigid body this node is attached to
					FERigidBody& RBi = m_fem.m_RB[nodei.m_rid];

					lmi = m_fem.m_RB[nodei.m_rid].m_LM;
					
					// get the relative distance
					ai = nodei.m_rt - RBi.m_rt;
	
					Ri[0][1] = ai.z; Ri[0][2] =-ai.y;
					Ri[1][0] =-ai.z; Ri[1][2] = ai.x;
					Ri[2][0] = ai.y; Ri[2][1] =-ai.x;

					// Kij
					KR[0][0] = kij[0][0]; KR[0][1] = kij[0][1]; KR[0][2] = kij[0][2];
					KR[1][0] = kij[1][0]; KR[1][1] = kij[1][1]; KR[1][2] = kij[1][2];
					KR[2][0] = kij[2][0]; KR[2][1] = kij[2][1]; KR[2][2] = kij[2][2];


					//Kij*Rj
					KR[0][3] = kij[0][0]*Rj[0][0]+kij[0][1]*Rj[1][0]+kij[0][2]*Rj[2][0];
					KR[0][4] = kij[0][0]*Rj[0][1]+kij[0][1]*Rj[1][1]+kij[0][2]*Rj[2][1];
					KR[0][5] = kij[0][0]*Rj[0][2]+kij[0][1]*Rj[1][2]+kij[0][2]*Rj[2][2];

					KR[1][3] = kij[1][0]*Rj[0][0]+kij[1][1]*Rj[1][0]+kij[1][2]*Rj[2][0];
					KR[1][4] = kij[1][0]*Rj[0][1]+kij[1][1]*Rj[1][1]+kij[1][2]*Rj[2][1];
					KR[1][5] = kij[1][0]*Rj[0][2]+kij[1][1]*Rj[1][2]+kij[1][2]*Rj[2][2];

					KR[2][3] = kij[2][0]*Rj[0][0]+kij[2][1]*Rj[1][0]+kij[2][2]*Rj[2][0];
					KR[2][4] = kij[2][0]*Rj[0][1]+kij[2][1]*Rj[1][1]+kij[2][2]*Rj[2][1];
					KR[2][5] = kij[2][0]*Rj[0][2]+kij[2][1]*Rj[1][2]+kij[2][2]*Rj[2][2];


					// Ri^T*Kij
					KR[3][0] = Ri[0][0]*kij[0][0]+Ri[1][0]*kij[1][0]+Ri[2][0]*kij[2][0];
					KR[3][1] = Ri[0][0]*kij[0][1]+Ri[1][0]*kij[1][1]+Ri[2][0]*kij[2][1];
					KR[3][2] = Ri[0][0]*kij[0][2]+Ri[1][0]*kij[1][2]+Ri[2][0]*kij[2][2];

					KR[4][0] = Ri[0][1]*kij[0][0]+Ri[1][1]*kij[1][0]+Ri[2][1]*kij[2][0];
					KR[4][1] = Ri[0][1]*kij[0][1]+Ri[1][1]*kij[1][1]+Ri[2][1]*kij[2][1];
					KR[4][2] = Ri[0][1]*kij[0][2]+Ri[1][1]*kij[1][2]+Ri[2][1]*kij[2][2];

					KR[5][0] = Ri[0][2]*kij[0][0]+Ri[1][2]*kij[1][0]+Ri[2][2]*kij[2][0];
					KR[5][1] = Ri[0][2]*kij[0][1]+Ri[1][2]*kij[1][1]+Ri[2][2]*kij[2][1];
					KR[5][2] = Ri[0][2]*kij[0][2]+Ri[1][2]*kij[1][2]+Ri[2][2]*kij[2][2];



					// Ri^T*Kij*Rj
					KR[3][3] = Ri[0][0]*KR[0][3]+Ri[1][0]*KR[1][3]+Ri[2][0]*KR[2][3];
					KR[3][4] = Ri[0][0]*KR[0][4]+Ri[1][0]*KR[1][4]+Ri[2][0]*KR[2][4];
					KR[3][5] = Ri[0][0]*KR[0][5]+Ri[1][0]*KR[1][5]+Ri[2][0]*KR[2][5];

					KR[4][3] = Ri[0][1]*KR[0][3]+Ri[1][1]*KR[1][3]+Ri[2][1]*KR[2][3];
					KR[4][4] = Ri[0][1]*KR[0][4]+Ri[1][1]*KR[1][4]+Ri[2][1]*KR[2][4];
					KR[4][5] = Ri[0][1]*KR[0][5]+Ri[1][1]*KR[1][5]+Ri[2][1]*KR[2][5];

					KR[5][3] = Ri[0][2]*KR[0][3]+Ri[1][2]*KR[1][3]+Ri[2][2]*KR[2][3];
					KR[5][4] = Ri[0][2]*KR[0][4]+Ri[1][2]*KR[1][4]+Ri[2][2]*KR[2][4];
					KR[5][5] = Ri[0][2]*KR[0][5]+Ri[1][2]*KR[1][5]+Ri[2][2]*KR[2][5];

					// add the stiffness components to the Krr matrix
					for (k=0; k<6; ++k)
						for (l=0; l<6; ++l)
						{
							J = lmj[k];
							I = lmi[l];

							if (I >= 0)
							{
								if (J < -1) m_Fd[I] -= KR[l][k]*m_ui[-J-2];
								else if (J >= 0) K.add(I,J, KR[l][k]);
							}
						}

					// we still need to couple the non-rigid degrees of node i to the
					// rigid dofs of node j
					for (k=3; k<ndof; ++k)
						for (l=0; l<3; ++l)
						{
							KF[k][l] = kij[k][l];
							KF[k][3+l] = kij[k][0]*Rj[0][l] + kij[k][1]*Rj[1][l] + kij[k][2]*Rj[2][l];
						}

					for (k=0; k<6; ++k)
						for (l=3; l<ndof; ++l)
						{
							J = lmj[k];
							I = elm[ndof*i+l];

							if (I >= 0)
							{
								if (J < -1) m_Fd[I] -= KF[l][k]*m_ui[-J-2];
								else if (J >= 0) K.add(I,J, KF[l][k]);
							}
						}

				}
				else
				{
					// node i is not a rigid body node
					// add the stiffness components to the Kfr matrix

					// Kij
					for (k=0; k<ndof; ++k)
						for (l=0; l<3; ++l)
						{
							KF[k][l] = kij[k][l];
							KF[k][3+l] = kij[k][0]*Rj[0][l] + kij[k][1]*Rj[1][l] + kij[k][2]*Rj[2][l];
						}

					for (k=0; k<6; ++k)
						for (l=0; l<ndof; ++l)
						{
							J = lmj[k];
							I = elm[ndof*i+l];

							if (I >= 0)
							{
								if (J < -1) m_Fd[I] -= KF[l][k]*m_ui[-J-2];
								else if (J >= 0) K.add(I,J, KF[l][k]);
							}
						}
				}
			}
		}
		else
		{
			// loop over rows
			for (i=0; i<n; ++i)
			{
				FENode& nodei = m_fem.m_mesh.Node(en[i]);
				if (nodei.m_rid>=0)
				{
					// node i is a rigid body
					// get the rigid body this node is attached to
					FERigidBody& RBi = m_fem.m_RB[nodei.m_rid];

					// get the rigid body equation nrs.
					lmi = RBi.m_LM;

					// get the relative distance to the center of mass
					ai = nodei.m_rt - RBi.m_rt;

					Ri[0][1] = ai.z; Ri[0][2] =-ai.y;
					Ri[1][0] =-ai.z; Ri[1][2] = ai.x;
					Ri[2][0] = ai.y; Ri[2][1] =-ai.x;

					// get the element sub-matrix
					for (k=0; k<ndof; ++k)
						for (l=0; l<ndof; ++l)
							kij[k][l] = ke[ndof*i+k][ndof*j+l];

					// add the stiffness components to the Krf matrix

					// Kij
					for (k=0; k<ndof; ++k)
						for (l=0; l<3; ++l)
						{
							KF[k][l] = kij[l][k];
							KF[k][3+l] = Ri[0][l]*kij[0][k] + Ri[1][l]*kij[1][k] + Ri[2][l]*kij[2][k];
						}

					for (k=0; k<6; ++k)
						for (l=0; l<ndof; ++l)
						{
							I = lmi[k];
							J = elm[ndof*j+l];

							if (I >= 0)
							{
								if (J < -1) m_Fd[I] -= KF[l][k]*m_ui[-J-2];
								else if (J >= 0) K.add(I,J, KF[l][k]);
							}
						}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//!  Assembles the element stiffness matrix into the global stiffness matrix.
//!  Also adjusts the global stiffness matrix and residual to take the 
//!  prescribed displacements into account.

// TODO: In stead of changing the global stiffness matrix to accomodate for 
// the rigid bodies and linear constraints, can I modify the element stiffness
// matrix prior to assembly? I might have to change the elm vector as well as 
// the element matrix size.

void FESolver::AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke)
{
	// assemble into global stiffness matrix
	m_pK->Assemble(ke, elm);

	// adjust for linear constraints
	if (m_fem.m_LinC.size() > 0)
	{
		int i, j, l;
		int nlin = m_fem.m_LinC.size();

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
							if ((J>=0) && (J<m_fem.m_nreq) && (I>=0)) m_Fd[I] -= kij*m_ui[J];
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
							if ((J>=0) && (J<m_fem.m_nreq) && (I>=0)) m_Fd[I] -= kij*m_ui[J];
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
								if ((J>=0) && (J<m_fem.m_nreq) && (I>=0)) m_Fd[I] -= kij*m_ui[J];
							}
						}
					}
				}
			}
		}
	}

	// adjust stiffness matrix for prescribed degrees of freedom
	if (m_fem.m_DC.size() > 0)
	{
		int i, j;
		int I, J;

		SparseMatrix& K = *m_pK;

		int N = ke.rows();

		int neq = m_fem.m_neq;

		// loop over columns
		for (j=0; j<N; ++j)
		{
			J = -elm[j]-2;
			if ((J >= 0) && (J<m_fem.m_nreq))
			{
				// dof j is a prescribed degree of freedom

				// loop over rows
				for (i=0; i<N; ++i)
				{
					I = elm[i];
					if (I >= 0)
					{
						// dof i is not a prescribed degree of freedom
						m_Fd[I] -= ke[i][j]*m_ui[J];
					}
				}

				// set the diagonal element of K to 1
				K.set(J,J, 1);			
			}
		}
	}

	// see if there are any rigid body dofs here
	if (m_fem.m_nrb > 0) RigidStiffness(en, elm, ke);
}

//-----------------------------------------------------------------------------
//! Calculates the contact forces
void FESolver::ContactForces(vector<double>& R)
{
	for (int i=0; i<m_fem.ContactInterfaces(); ++i) m_fem.m_CI[i].ContactForces(R);
}

//-----------------------------------------------------------------------------
//! calculates the residual vector
//! Note that the concentrated nodal forces are not calculated here.
//! This is because they do not depend on the geometry 
//! so we only calculate them once (in Quasin) and then add them here.

bool FESolver::Residual(vector<double>& R)
{
	int i, j, J, neln;
	int neq = m_fem.m_neq;
	int ndof;

	// initialize residual with concentrated nodal loads
	R = m_Fn;

	// zero nodal reaction forces
	m_Fr.zero();

	// element force vector
	vector<double> fe;

	// zero rigid body reaction forces
	for (i=0; i<m_fem.m_nrb; ++i)
	{
		FERigidBody& RB = m_fem.m_RB[i];
		RB.m_Fr = RB.m_Mr = vec3d(0,0,0);
	}

	// get the mesh
	FEMesh& mesh = m_fem.m_mesh;

	// loop over solid elements
	int NE = mesh.SolidElements();
	for (i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = mesh.SolidElement(i);

		// unpack the element
		if (!el.IsRigid()) mesh.UnpackElement(el);

		FEMaterial* pm = m_fem.GetMaterial(el.GetMatID());

		// get the element force vector and initialize it to zero
		ndof = 3*el.Nodes();
		fe.create(ndof);
		fe.zero();

		// skip rigid elements for internal force calculations
		if (!el.IsRigid())
		{
			// calculate internal force vector
			if (el.Type() != FE_UDGHEX) InternalForces(el, fe);
			else UDGInternalForces(el, fe);

			// apply body forces
			if (m_fem.UseBodyForces())
			{
				BodyForces(el, fe);
			}

			// assemble element 'fe'-vector into global R vector
			AssembleResidual(el.m_node, el.LM(), fe, R);
		}
		else if (m_fem.UseBodyForces())
		{
			// unpack the element
			mesh.UnpackElement(el);

			// apply body force to rigid elements
			BodyForces(el, fe);

			// assemble element 'fe'-vector into global R vector
			AssembleResidual(el.m_node, el.LM(), fe, R);
		}

		// do poro-elastic forces
		if ((m_fem.m_pStep->m_itype == FE_STATIC_PORO)&&(dynamic_cast<FEPoroElastic*>(pm)))
		{
			// calculate fluid internal work
			InternalFluidWork(el, fe);

			// add fluid work to global residual
			neln = el.Nodes();
			int* lm = el.LM();
			for (j=0; j<neln; ++j)
			{
				J = lm[3*neln+j];
				if (J >= 0) R[J] += fe[j];
			}
		}
	}

	// loop over shell elements
	int NS = mesh.ShellElements();
	for (i=0; i<NS; ++i)
	{
		// get the element
		FEShellElement& el = mesh.ShellElement(i);

		// create the element force vector and initialize to zero
		ndof = 6*el.Nodes();
		fe.create(ndof);
		fe.zero();

		if (!el.IsRigid())
		{
			mesh.UnpackElement(el);

			// skip rigid elements for internal force calculation
			InternalForces(el, fe);

			// apply body forces to shells
			if (m_fem.UseBodyForces())
			{
				BodyForces(el, fe);
			}

			// assemble the residual
			AssembleResidual(el.m_node, el.LM(), fe, R);
		}
		else if (m_fem.UseBodyForces())
		{
			mesh.UnpackElement(el);

			// apply body forces to shells
			BodyForces(el, fe);

			// assemble the residual
			AssembleResidual(el.m_node, el.LM(), fe, R);
		}

		// TODO: Do poro-elasticity for shells
	}

	// calculate inertial forces for dynamic problems
	if (m_fem.m_pStep->m_itype == FE_DYNAMIC) InertialForces(R);

	// calculate forces due to pressure and add them to the residual
	// loop over surface elements
	int npr = m_fem.m_PC.size();
	for (i=0; i<npr; ++i)
	{
		FESurfaceElement& el = m_fem.m_psurf->Element(i);
		mesh.UnpackElement(el);

		// calculate nodal pressures
		double* pt = el.pt();
		FEPressureLoad& pc = m_fem.m_PC[i];

		double g = m_fem.GetLoadCurve(pc.lc)->Value();

		for (j=0; j<el.Nodes(); ++j) pt[j] = -g*pc.s[j];

		ndof = 3*el.Nodes();
		fe.create(ndof);

		if (pc.blinear) LinearPressureForce(el, fe); else PressureForce(el, fe);

		// add element force vector to global force vector
		AssembleResidual(el.m_node, el.LM(), fe, R);
	}

	// rigid joint forces
	if (m_fem.m_nrj)
	{
		for (int i=0; i<m_fem.m_nrj; ++i) m_fem.m_RJ[i].JointForces(R);
	}

	// calculate contact forces
	if (m_fem.m_bcontact)
	{
		ContactForces(R);
	}

	// calculate linear constraint forces
	// note that these are the linear constraints
	// enforced using the augmented lagrangian
	LinearConstraintForces(R);

	// add discrete element forces
	if (m_fem.m_DE.size())
	{
		DiscreteElementForces(R);
	}

	// increase RHS counter
	m_nrhs++;

	return true;
}

//-----------------------------------------------------------------------------
//! calculate the linear constraint forces 

void FESolver::LinearConstraintForces(vector<double> &R)
{
	int N = m_fem.m_LCSet.size();
	if (N>0)
	{
		list<FELinearConstraintSet*>::iterator im = m_fem.m_LCSet.begin();
		for (int i=0; i<N; ++i, ++im) (*im)->Residual(R);
	}
}

//-----------------------------------------------------------------------------
//!  Assembles the element into the global residual. This function
//!  also checks for rigid dofs and assembles the residual using a condensing
//!  procedure in the case of rigid dofs.

void FESolver::AssembleResidual(vector<int>& en, vector<int>& elm, vector<double>& fe, vector<double>& R)
{
	int i, j, I, n, l;

	int neq = m_fem.m_neq;

	vec3d a, d;

	// assemble the element residual into the global residual
	int ndof = fe.size();
	for (i=0; i<ndof; ++i)
	{
		I = elm[i];
		if ( I >= 0) R[I] += fe[i];
		else if (-I-2 >= 0) m_Fr[-I-2] -= fe[i];
	}

	int ndn = ndof / en.size();

	// if there are linear constraints we need to apply them
	if (m_fem.m_LinC.size() > 0)
	{
		// loop over all degrees of freedom of this element
		for (i=0; i<ndof; ++i)
		{
			// see if this dof belongs to a linear constraint
			n = MAX_NDOFS*(en[i/ndn]) + i%ndn;
			l = m_fem.m_LCT[n];
			if (l >= 0)
			{
				// if so, get the linear constraint
				FELinearConstraint& lc = *m_fem.m_LCA[l];
				assert(elm[i] == -1);
	
				// now loop over all "slave" nodes and
				// add the contribution to the residual
				int ns = lc.slave.size();
				list<FELinearConstraint::SlaveDOF>::iterator is = lc.slave.begin();
				for (j=0; j<ns; ++j, ++is)
				{
					I = is->neq;
					if (I >= 0)
					{
						double A = is->val;
						R[I] += A*fe[i];
					}
				}
			}
		}
	}

	// If there are rigid bodies we need to look for rigid dofs
	if (m_fem.m_nrb > 0)
	{
		int *lm;

		for (i=0; i<ndof; i+=ndn)
		{
			FENode& node = m_fem.m_mesh.Node(en[i/ndn]);
			if (node.m_rid >= 0)
			{
				vec3d F(fe[i], fe[i+1], fe[i+2]);

				// this is an interface dof
				// get the rigid body this node is connected to
				FERigidBody& RB = m_fem.m_RB[node.m_rid];
				lm = RB.m_LM;

				// add to total torque of this body
				a = node.m_rt - RB.m_rt;

				n = lm[3]; if (n >= 0) R[n] += a.y*F.z-a.z*F.y; RB.m_Mr.x -= a.y*F.z-a.z*F.y;
				n = lm[4]; if (n >= 0) R[n] += a.z*F.x-a.x*F.z; RB.m_Mr.y -= a.z*F.x-a.x*F.z;
				n = lm[5]; if (n >= 0) R[n] += a.x*F.y-a.y*F.x; RB.m_Mr.z -= a.x*F.y-a.y*F.x;
/*
				// if the rotational degrees of freedom are constrained for a rigid node
				// then we need to add an additional component to the residual
				if (node.m_ID[7] == lm[3])
				{
					d = node.m_Dt;
					n = lm[3]; if (n >= 0) R[n] += d.y*F.z-d.z*F.y; RB.m_Mr.x -= d.y*F.z-d.z*F.y;
					n = lm[4]; if (n >= 0) R[n] += d.z*F.x-d.x*F.z; RB.m_Mr.y -= d.z*F.x-d.x*F.z;
					n = lm[5]; if (n >= 0) R[n] += d.x*F.y-d.y*F.x; RB.m_Mr.z -= d.x*F.y-d.y*F.x;
				}
*/
				// add to global force vector
				n = lm[0]; if (n >= 0) R[n] += F.x; RB.m_Fr.x -= F.x;
				n = lm[1]; if (n >= 0) R[n] += F.y; RB.m_Fr.y -= F.y;
				n = lm[2]; if (n >= 0) R[n] += F.z; RB.m_Fr.z -= F.z;
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates element inertial stiffness matrix

void FESolver::ElementInertialStiffness(FESolidElement& el, matrix& ke)
{
	int i, j, n;

	// shape functions
	double* H;

	// jacobian
	double detJ0;

	// get the material
	FEMaterial* pm = m_fem.GetMaterial(el.GetMatID());

	double a = 4.0 / (m_fem.m_pStep->m_dt*m_fem.m_pStep->m_dt);
	double d = pm->Density();
	double kab;

	// Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	const int ndof = 3*neln;

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// calculate element stiffness matrix
	for (n=0; n<nint; ++n)
	{
		H = el.H(n);
		detJ0 = el.detJ0(n)*gw[n];
		for (i=0; i<neln; ++i)
			for (j=i; j<neln; ++j)
			{
				kab = a*H[i]*H[j]*detJ0*d;
				ke[3*i  ][3*j  ] += kab;
				ke[3*i+1][3*j+1] += kab;
				ke[3*i+2][3*j+2] += kab;
			}	
	}

	// assign symmetic parts
	// TODO: Can this be omitted by changing the Assemble routine so that it only
	// grabs elements from the upper diagonal matrix?
	for (i=0; i<ndof; ++i)
		for (j=i+1; j<ndof; ++j)
			ke[j][i] = ke[i][j];

}

//-----------------------------------------------------------------------------
//! This function calculates the element stiffness matrix. It calls the material
//! stiffness function, the geometrical stiffness function and, if necessary, the
//! dilatational stiffness function. Note that these three functions only calculate
//! the upper diagonal matrix due to the symmetry of the element stiffness matrix
//! The last section of this function fills the rest of the element stiffness matrix.

void FESolver::ElementStiffness(FESolidElement& el, matrix& ke)
{
	// see if the material is incompressible
	FEElasticMaterial* pme = m_fem.GetElasticMaterial(el.GetMatID());
	bool bdilst = false;
	if (dynamic_cast<FEIncompressibleMaterial*>(pme)) bdilst = true;

	if (el.Type() == FE_UDGHEX)
	{
		// calculate material stiffness
		UDGMaterialStiffness(el, ke);

		// calculate geometrical stiffness
		UDGGeometricalStiffness(el, ke);

		// Calculate dilatational stiffness, if necessary
		if (bdilst) UDGDilatationalStiffness(el, ke);

		// add hourglass stiffness
		UDGHourglassStiffness(el, ke);
	}
	else
	{
		// calculate material stiffness (i.e. constitutive component)
		MaterialStiffness(el, ke);

		// calculate geometrical stiffness
		GeometricalStiffness(el, ke);

		// Calculate dilatational stiffness, if necessary
		if (bdilst) DilatationalStiffness(el, ke);
	}

	// assign symmetic parts
	// TODO: Can this be omitted by changing the Assemble routine so that it only
	// grabs elements from the upper diagonal matrix?
	int ndof = 3*el.Nodes();
	int i, j;
	for (i=0; i<ndof; ++i)
		for (j=i+1; j<ndof; ++j)
			ke[j][i] = ke[i][j];
}

//-----------------------------------------------------------------------------
//! Calculates element material stiffness element matrix

void FESolver::MaterialStiffness(FESolidElement &el, matrix &ke)
{
	int i, i3, j, j3, n;

	// Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	const int ndof = 3*neln;

	// global derivatives of shape functions
	// NOTE: hard-coding of hex elements!
	// Gx = dH/dx
	double Gx[8], Gy[8], Gz[8];

	double Gxi, Gyi, Gzi;
	double Gxj, Gyj, Gzj;

	// The 'D' matrix
	double D[6][6] = {0};	// The 'D' matrix

	// The 'D*BL' matrix
	double DBL[6][3];

	double *Grn, *Gsn, *Gtn;
	double Gr, Gs, Gt;

	// jacobian
	double Ji[3][3], detJt;
	
	// weights at gauss points
	const double *gw = el.GaussWeights();

	// see if this is a poroelastic material
	FEMaterial* pmat = m_fem.GetMaterial(el.GetMatID());
	bool bporo = false;
	if ((m_fem.m_pStep->m_itype == FE_STATIC_PORO) && (dynamic_cast<FEPoroElastic*>(pmat))) bporo = true;

	// calculate element stiffness matrix
	for (n=0; n<nint; ++n)
	{
		// calculate jacobian
		el.invjact(Ji, n);
		detJt = el.detJt(n)*gw[n];

		Grn = el.Gr(n);
		Gsn = el.Gs(n);
		Gtn = el.Gt(n);

		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// setup the material point
		// NOTE: deformation gradient and determinant have already been evaluated in the stress routine
//		el.defgrad(pt.F, n);
//		pt.J = el.detF(n);
		pt.avgJ = el.m_eJ;
		pt.avgp = el.m_ep;

		if (bporo)
		{
			FEPoroElasticMaterialPoint& pt = *(mp.ExtractData<FEPoroElasticMaterialPoint>());

			// evaluate fluid pressure at gauss-point
			pt.m_p = el.Evaluate(el.pt(), n);
		}

		// get the 'D' matrix
		tens4ds C = pmat->Tangent(mp);
		C.extract(D);

/*		if (m_fem.GetDebugFlag())
		{
			tens4ds t(D);
			if (IsPositiveDefinite(t) == false)
			{
				m_fem.m_log.printbox("WARNING", "Elasticity tensor is not positive-definite for\nelement %d at integration point %d.", el.m_nID, n+1);
			}
		}
*/
		for (i=0; i<neln; ++i)
		{
			Gr = Grn[i];
			Gs = Gsn[i];
			Gt = Gtn[i];

			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx[i] = Ji[0][0]*Gr+Ji[1][0]*Gs+Ji[2][0]*Gt;
			Gy[i] = Ji[0][1]*Gr+Ji[1][1]*Gs+Ji[2][1]*Gt;
			Gz[i] = Ji[0][2]*Gr+Ji[1][2]*Gs+Ji[2][2]*Gt;
		}

		// we only calculate the upper triangular part
		// since ke is symmetric. The other part is
		// determined below using this symmetry.
		for (i=0, i3=0; i<neln; ++i, i3 += 3)
		{
			Gxi = Gx[i];
			Gyi = Gy[i];
			Gzi = Gz[i];

			for (j=i, j3 = i3; j<neln; ++j, j3 += 3)
			{
				Gxj = Gx[j];
				Gyj = Gy[j];
				Gzj = Gz[j];

				// calculate D*BL matrices
				DBL[0][0] = (D[0][0]*Gxj+D[0][3]*Gyj+D[0][5]*Gzj);
				DBL[0][1] = (D[0][1]*Gyj+D[0][3]*Gxj+D[0][4]*Gzj);
				DBL[0][2] = (D[0][2]*Gzj+D[0][4]*Gyj+D[0][5]*Gxj);

				DBL[1][0] = (D[1][0]*Gxj+D[1][3]*Gyj+D[1][5]*Gzj);
				DBL[1][1] = (D[1][1]*Gyj+D[1][3]*Gxj+D[1][4]*Gzj);
				DBL[1][2] = (D[1][2]*Gzj+D[1][4]*Gyj+D[1][5]*Gxj);

				DBL[2][0] = (D[2][0]*Gxj+D[2][3]*Gyj+D[2][5]*Gzj);
				DBL[2][1] = (D[2][1]*Gyj+D[2][3]*Gxj+D[2][4]*Gzj);
				DBL[2][2] = (D[2][2]*Gzj+D[2][4]*Gyj+D[2][5]*Gxj);

				DBL[3][0] = (D[3][0]*Gxj+D[3][3]*Gyj+D[3][5]*Gzj);
				DBL[3][1] = (D[3][1]*Gyj+D[3][3]*Gxj+D[3][4]*Gzj);
				DBL[3][2] = (D[3][2]*Gzj+D[3][4]*Gyj+D[3][5]*Gxj);

				DBL[4][0] = (D[4][0]*Gxj+D[4][3]*Gyj+D[4][5]*Gzj);
				DBL[4][1] = (D[4][1]*Gyj+D[4][3]*Gxj+D[4][4]*Gzj);
				DBL[4][2] = (D[4][2]*Gzj+D[4][4]*Gyj+D[4][5]*Gxj);

				DBL[5][0] = (D[5][0]*Gxj+D[5][3]*Gyj+D[5][5]*Gzj);
				DBL[5][1] = (D[5][1]*Gyj+D[5][3]*Gxj+D[5][4]*Gzj);
				DBL[5][2] = (D[5][2]*Gzj+D[5][4]*Gyj+D[5][5]*Gxj);

				ke[i3  ][j3  ] += (Gxi*DBL[0][0] + Gyi*DBL[3][0] + Gzi*DBL[5][0] )*detJt;
				ke[i3  ][j3+1] += (Gxi*DBL[0][1] + Gyi*DBL[3][1] + Gzi*DBL[5][1] )*detJt;
				ke[i3  ][j3+2] += (Gxi*DBL[0][2] + Gyi*DBL[3][2] + Gzi*DBL[5][2] )*detJt;

				ke[i3+1][j3  ] += (Gyi*DBL[1][0] + Gxi*DBL[3][0] + Gzi*DBL[4][0] )*detJt;
				ke[i3+1][j3+1] += (Gyi*DBL[1][1] + Gxi*DBL[3][1] + Gzi*DBL[4][1] )*detJt;
				ke[i3+1][j3+2] += (Gyi*DBL[1][2] + Gxi*DBL[3][2] + Gzi*DBL[4][2] )*detJt;

				ke[i3+2][j3  ] += (Gzi*DBL[2][0] + Gyi*DBL[4][0] + Gxi*DBL[5][0] )*detJt;
				ke[i3+2][j3+1] += (Gzi*DBL[2][1] + Gyi*DBL[4][1] + Gxi*DBL[5][1] )*detJt;
				ke[i3+2][j3+2] += (Gzi*DBL[2][2] + Gyi*DBL[4][2] + Gxi*DBL[5][2] )*detJt;
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates element's geometrical stiffness component for integration point n

void FESolver::GeometricalStiffness(FESolidElement &el, matrix &ke)
{
	int n, i, j;

	double Gx[8], Gy[8], Gz[8];
	double *Grn, *Gsn, *Gtn;
	double Gr, Gs, Gt;

	// nr of nodes
	int neln = el.Nodes();

	// nr of integration points
	int nint = el.GaussPoints();

	// jacobian
	double Ji[3][3], detJt;

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// stiffness component for the initial stress component of stiffness matrix
	double kab;

	// calculate geometrical element stiffness matrix
	for (n=0; n<nint; ++n)
	{
		// calculate jacobian
		el.invjact(Ji, n);
		detJt = el.detJt(n)*gw[n];

		Grn = el.Gr(n);
		Gsn = el.Gs(n);
		Gtn = el.Gt(n);

		for (i=0; i<neln; ++i)
		{
			Gr = Grn[i];
			Gs = Gsn[i];
			Gt = Gtn[i];

			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx[i] = Ji[0][0]*Gr+Ji[1][0]*Gs+Ji[2][0]*Gt;
			Gy[i] = Ji[0][1]*Gr+Ji[1][1]*Gs+Ji[2][1]*Gt;
			Gz[i] = Ji[0][2]*Gr+Ji[1][2]*Gs+Ji[2][2]*Gt;
		}

		// get the material point data
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// element's Cauchy-stress tensor at gauss point n
		// s is the voight vector
		mat3ds& s = pt.s;

		for (i=0; i<neln; ++i)
			for (j=i; j<neln; ++j)
			{
				kab = (Gx[i]*(s.xx()*Gx[j]+s.xy()*Gy[j]+s.xz()*Gz[j]) +
					   Gy[i]*(s.xy()*Gx[j]+s.yy()*Gy[j]+s.yz()*Gz[j]) + 
					   Gz[i]*(s.xz()*Gx[j]+s.yz()*Gy[j]+s.zz()*Gz[j]))*detJt;

				ke[3*i  ][3*j  ] += kab;
				ke[3*i+1][3*j+1] += kab;
				ke[3*i+2][3*j+2] += kab;
			}
	}
}

//-----------------------------------------------------------------------------
//! calculates dilatational element stiffness component for element iel

void FESolver::DilatationalStiffness(FESolidElement& elem, matrix& ke)
{
	int i, j, n;

	const int nint = elem.GaussPoints();
	const int neln = elem.Nodes();
	const int ndof = 3*neln;

	// get the elements material
	FEElasticMaterial* pm = m_fem.GetElasticMaterial(elem.GetMatID());

	FEIncompressibleMaterial* pmi = dynamic_cast<FEIncompressibleMaterial*>(pm);
	assert(pmi);

	// average global derivatives
	vector<double> gradN(3*neln);
	gradN.zero();

	// initial element volume
	double Ve = 0;

	// global derivatives of shape functions
	double Gx, Gy, Gz;
	const double *gw = elem.GaussWeights();

	// jacobian
	double Ji[3][3], detJt, detJ0;

	double *Gr, *Gs, *Gt;

	// repeat over gauss-points
	for (n=0; n<nint; ++n)
	{
		// calculate jacobian
		detJ0 = elem.detJ0(n);
		detJt = elem.detJt(n);
		elem.invjact(Ji, n);

		detJt *= gw[n];

		Ve += detJ0*gw[n];

		Gr = elem.Gr(n);
		Gs = elem.Gs(n);
		Gt = elem.Gt(n);

		// calculate global gradient of shape functions
		// note that we need the transposed of Ji, not Ji itself !
		for (i=0; i<neln; ++i)
		{
			Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

			gradN[3*i  ] += Gx*detJt;
			gradN[3*i+1] += Gy*detJt;
			gradN[3*i+2] += Gz*detJt;
		}
	}

	// get effective modulus
	double k = pmi->Upp(elem.m_eJ);

	// next, we add the Lagrangian contribution
	// note that this term will always be zero if the material does not
	// use the augmented lagrangian
	k += elem.m_Lk*pmi->hpp(elem.m_eJ);

	// divide by initial volume
	k /= Ve;

	// calculate dilatational stiffness component
	// we only calculate the upper triangular part
	// since ke is symmetric.
	for (i=0; i<ndof; ++i)
		for (j=i; j<ndof; ++j)
			ke[i][j] += k*gradN[i]*gradN[j];
}

//-----------------------------------------------------------------------------

void FESolver::UDGMaterialStiffness(FESolidElement &el, matrix &ke)
{
	// make sure we have the right element type
	assert(el.Type() == FE_UDGHEX);

	int i, i3, j, j3;

	// Get the current element's data
	const int neln = el.Nodes();

	double Gxi, Gyi, Gzi;
	double Gxj, Gyj, Gzj;

	// The 'D' matrix
	double D[6][6] = {0};	// The 'D' matrix

	// The 'D*BL' matrix
	double DBL[6][3];

	// see if this is a poroelastic material
	FEMaterial* pmat = m_fem.GetMaterial(el.GetMatID());
	bool bporo = false;
	if ((m_fem.m_pStep->m_itype == FE_STATIC_PORO) && (dynamic_cast<FEPoroElastic*>(pmat))) bporo = true;
	
	// for now we do not allow this element to be used in a poroelastic simulation
	assert(bporo == false);

	FEMaterialPoint& mp = *el.m_State[0];
	FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

	// calculate the average cartesian derivatives
	double GX[8], GY[8], GZ[8];
	AvgCartDerivs(el, GX, GY, GZ);

	// calculate average deformation gradient Fbar
	mat3d Fb;
	AvgDefGrad(el, Fb, GX, GY, GZ);

	// calculate the transposed inverse of Fbar
	mat3d Fti = Fb.transinv();

	// calculate current element volume
	double ve = HexVolume(el, 1);

	// current averaged shape derivatives
	double Gx[8], Gy[8], Gz[8];
	for (i=0; i<8; ++i)
	{
		Gx[i] = Fti(0,0)*GX[i]+Fti(0,1)*GY[i]+Fti(0,2)*GZ[i];
		Gy[i] = Fti(1,0)*GX[i]+Fti(1,1)*GY[i]+Fti(1,2)*GZ[i];
		Gz[i] = Fti(2,0)*GX[i]+Fti(2,1)*GY[i]+Fti(2,2)*GZ[i];
	}

	// setup the material point
	// NOTE: deformation gradient and determinant have already been evaluated in the stress routine
//	el.defgrad(pt.F, n);
//	pt.J = el.detF(n);
	pt.avgJ = el.m_eJ;
	pt.avgp = el.m_ep;

	// get the 'D' matrix
	tens4ds C = pmat->Tangent(mp);
	C.extract(D);

	// we only calculate the upper triangular part
	// since ke is symmetric. The other part is
	// determined below using this symmetry.
	for (i=0, i3=0; i<neln; ++i, i3 += 3)
	{
		Gxi = Gx[i];
		Gyi = Gy[i];
		Gzi = Gz[i];

		for (j=i, j3 = i3; j<neln; ++j, j3 += 3)
		{
			Gxj = Gx[j];
			Gyj = Gy[j];
			Gzj = Gz[j];

			// calculate D*BL matrices
			DBL[0][0] = (D[0][0]*Gxj+D[0][3]*Gyj+D[0][5]*Gzj);
			DBL[0][1] = (D[0][1]*Gyj+D[0][3]*Gxj+D[0][4]*Gzj);
			DBL[0][2] = (D[0][2]*Gzj+D[0][4]*Gyj+D[0][5]*Gxj);

			DBL[1][0] = (D[1][0]*Gxj+D[1][3]*Gyj+D[1][5]*Gzj);
			DBL[1][1] = (D[1][1]*Gyj+D[1][3]*Gxj+D[1][4]*Gzj);
			DBL[1][2] = (D[1][2]*Gzj+D[1][4]*Gyj+D[1][5]*Gxj);

			DBL[2][0] = (D[2][0]*Gxj+D[2][3]*Gyj+D[2][5]*Gzj);
			DBL[2][1] = (D[2][1]*Gyj+D[2][3]*Gxj+D[2][4]*Gzj);
			DBL[2][2] = (D[2][2]*Gzj+D[2][4]*Gyj+D[2][5]*Gxj);

			DBL[3][0] = (D[3][0]*Gxj+D[3][3]*Gyj+D[3][5]*Gzj);
			DBL[3][1] = (D[3][1]*Gyj+D[3][3]*Gxj+D[3][4]*Gzj);
			DBL[3][2] = (D[3][2]*Gzj+D[3][4]*Gyj+D[3][5]*Gxj);

			DBL[4][0] = (D[4][0]*Gxj+D[4][3]*Gyj+D[4][5]*Gzj);
			DBL[4][1] = (D[4][1]*Gyj+D[4][3]*Gxj+D[4][4]*Gzj);
			DBL[4][2] = (D[4][2]*Gzj+D[4][4]*Gyj+D[4][5]*Gxj);

			DBL[5][0] = (D[5][0]*Gxj+D[5][3]*Gyj+D[5][5]*Gzj);
			DBL[5][1] = (D[5][1]*Gyj+D[5][3]*Gxj+D[5][4]*Gzj);
			DBL[5][2] = (D[5][2]*Gzj+D[5][4]*Gyj+D[5][5]*Gxj);

			ke[i3  ][j3  ] += (Gxi*DBL[0][0] + Gyi*DBL[3][0] + Gzi*DBL[5][0] )*ve;
			ke[i3  ][j3+1] += (Gxi*DBL[0][1] + Gyi*DBL[3][1] + Gzi*DBL[5][1] )*ve;
			ke[i3  ][j3+2] += (Gxi*DBL[0][2] + Gyi*DBL[3][2] + Gzi*DBL[5][2] )*ve;

			ke[i3+1][j3  ] += (Gyi*DBL[1][0] + Gxi*DBL[3][0] + Gzi*DBL[4][0] )*ve;
			ke[i3+1][j3+1] += (Gyi*DBL[1][1] + Gxi*DBL[3][1] + Gzi*DBL[4][1] )*ve;
			ke[i3+1][j3+2] += (Gyi*DBL[1][2] + Gxi*DBL[3][2] + Gzi*DBL[4][2] )*ve;

			ke[i3+2][j3  ] += (Gzi*DBL[2][0] + Gyi*DBL[4][0] + Gxi*DBL[5][0] )*ve;
			ke[i3+2][j3+1] += (Gzi*DBL[2][1] + Gyi*DBL[4][1] + Gxi*DBL[5][1] )*ve;
			ke[i3+2][j3+2] += (Gzi*DBL[2][2] + Gyi*DBL[4][2] + Gxi*DBL[5][2] )*ve;
		}
	}
}

//-----------------------------------------------------------------------------

void FESolver::UDGGeometricalStiffness(FESolidElement& el, matrix& ke)
{
	int i, j;

	// stiffness component for the initial stress component of stiffness matrix
	double kab;

	// get the material point data
	FEMaterialPoint& mp = *el.m_State[0];
	FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

	// element's Cauchy-stress tensor at gauss point n
	// s is the voight vector
	mat3ds& s = pt.s;

	// calculate the average cartesian derivatives
	double GX[8], GY[8], GZ[8];
	AvgCartDerivs(el, GX, GY, GZ);

	// calculate average deformation gradient Fbar
	mat3d Fb;
	AvgDefGrad(el, Fb, GX, GY, GZ);

	// calculate the transposed inverse of Fbar
	mat3d Fti = Fb.transinv();

	// calculate current element volume
	double ve = HexVolume(el, 1);

	// current averaged shape derivatives
	double Gx[8], Gy[8], Gz[8];
	for (i=0; i<8; ++i)
	{
		Gx[i] = Fti(0,0)*GX[i]+Fti(0,1)*GY[i]+Fti(0,2)*GZ[i];
		Gy[i] = Fti(1,0)*GX[i]+Fti(1,1)*GY[i]+Fti(1,2)*GZ[i];
		Gz[i] = Fti(2,0)*GX[i]+Fti(2,1)*GY[i]+Fti(2,2)*GZ[i];
	}

	for (i=0; i<8; ++i)
		for (j=i; j<8; ++j)
		{
			kab = (Gx[i]*(s.xx()*Gx[j]+s.xy()*Gy[j]+s.xz()*Gz[j]) +
				   Gy[i]*(s.xy()*Gx[j]+s.yy()*Gy[j]+s.yz()*Gz[j]) + 
				   Gz[i]*(s.xz()*Gx[j]+s.yz()*Gy[j]+s.zz()*Gz[j]))*ve;

			ke[3*i  ][3*j  ] += kab;
			ke[3*i+1][3*j+1] += kab;
			ke[3*i+2][3*j+2] += kab;
		}
}

//-----------------------------------------------------------------------------

void FESolver::UDGDilatationalStiffness(FESolidElement& el, matrix& ke)
{
	int i, j;

	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	const int ndof = 3*neln;

	// get the elements material
	FEElasticMaterial* pm = m_fem.GetElasticMaterial(el.GetMatID());

	FEIncompressibleMaterial* pmi = dynamic_cast<FEIncompressibleMaterial*>(pm);
	assert(pmi);

	// calculate the average cartesian derivatives
	double Gx[8], Gy[8], Gz[8];
	AvgCartDerivs(el, Gx, Gy, Gz, 1);

	// calculate element volume
	double ve = HexVolume(el, 1);
	double Ve = HexVolume(el, 0);

	// get effective modulus
	double k = pmi->Upp(el.m_eJ);

	// next, we add the Lagrangian contribution
	// note that this term will always be zero if the material does not
	// use the augmented lagrangian
	k += el.m_Lk*pmi->hpp(el.m_eJ);

	// multiply with volume
	k *= (ve*ve)/Ve;

	// calculate dilatational stiffness component
	// we only calculate the upper triangular part
	// since ke is symmetric.
	for (i=0; i<8; ++i)
		for (j=i; j<8; ++j)
		{
			ke[3*i  ][3*j  ] += k*Gx[i]*Gx[j];
			ke[3*i  ][3*j+1] += k*Gx[i]*Gy[j];
			ke[3*i  ][3*j+2] += k*Gx[i]*Gz[j];

			ke[3*i+1][3*j  ] += k*Gy[i]*Gx[j];
			ke[3*i+1][3*j+1] += k*Gy[i]*Gy[j];
			ke[3*i+1][3*j+2] += k*Gy[i]*Gz[j];

			ke[3*i+2][3*j  ] += k*Gz[i]*Gx[j];
			ke[3*i+2][3*j+1] += k*Gz[i]*Gy[j];
			ke[3*i+2][3*j+2] += k*Gz[i]*Gz[j];
		}
}

//-----------------------------------------------------------------------------
//! calculates the hourglass stiffness for UDG hex elements

void FESolver::UDGHourglassStiffness(FESolidElement& el, matrix& ke)
{
	int i, j;

	const double h4[8] = { 1,-1, 1,-1, 1,-1, 1,-1 };
	const double h5[8] = { 1,-1,-1, 1,-1, 1, 1,-1 };
	const double h6[8] = { 1, 1,-1,-1,-1,-1, 1, 1 };
	const double h7[8] = {-1, 1,-1, 1, 1,-1, 1,-1 };

	vec3d* r0 = el.r0();
	vec3d* rt = el.rt();

	double x4 = 0, x5 = 0, x6 = 0, x7 = 0;
	double y4 = 0, y5 = 0, y6 = 0, y7 = 0;
	double z4 = 0, z5 = 0, z6 = 0, z7 = 0;

	double X4 = 0, X5 = 0, X6 = 0, X7 = 0;
	double Y4 = 0, Y5 = 0, Y6 = 0, Y7 = 0;
	double Z4 = 0, Z5 = 0, Z6 = 0, Z7 = 0;

	for (i=0; i<8; ++i)
	{
		X4 += h4[i]*r0[i].x; Y4 += h4[i]*r0[i].y; Z4 += h4[i]*r0[i].z;
		X5 += h5[i]*r0[i].x; Y5 += h5[i]*r0[i].y; Z5 += h5[i]*r0[i].z;
		X6 += h6[i]*r0[i].x; Y6 += h6[i]*r0[i].y; Z6 += h6[i]*r0[i].z;
		X7 += h7[i]*r0[i].x; Y7 += h7[i]*r0[i].y; Z7 += h7[i]*r0[i].z;

		x4 += h4[i]*rt[i].x; y4 += h4[i]*rt[i].y; z4 += h4[i]*rt[i].z;
		x5 += h5[i]*rt[i].x; y5 += h5[i]*rt[i].y; z5 += h5[i]*rt[i].z;
		x6 += h6[i]*rt[i].x; y6 += h6[i]*rt[i].y; z6 += h6[i]*rt[i].z;
		x7 += h7[i]*rt[i].x; y7 += h7[i]*rt[i].y; z7 += h7[i]*rt[i].z;
	}

	double GX[8], GY[8], GZ[8];
	AvgCartDerivs(el, GX, GY, GZ);

	double g4[8] = {0}, g5[8] = {0}, g6[8] = {0}, g7[8] = {0};

	for (i=0; i<8; ++i)
	{
		g4[i] = h4[i] - (GX[i]*X4 + GY[i]*Y4 + GZ[i]*Z4);
		g5[i] = h5[i] - (GX[i]*X5 + GY[i]*Y5 + GZ[i]*Z5);
		g6[i] = h6[i] - (GX[i]*X6 + GY[i]*Y6 + GZ[i]*Z6);
		g7[i] = h7[i] - (GX[i]*X7 + GY[i]*Y7 + GZ[i]*Z7);
	}

	// calculate hourglass stiffness
	double hg = m_fem.m_pStep->m_hg;

	double kab;

	for (i=0; i<8; ++i)
	{
		for (j=i; j<8; ++j)
		{
			kab = hg*(g4[i]*g4[j] + g5[i]*g5[j] + g6[i]*g6[j] + g7[i]*g7[j]);

			ke[3*i  ][3*j  ] += kab;
			ke[3*i+1][3*j+1] += kab;
			ke[3*i+2][3*j+2] += kab;
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates dilatational element stiffness component for element iel

void FESolver::DilatationalStiffness(FEShellElement& elem, matrix& ke)
{
	int i, j, n;

	const int nint = elem.GaussPoints();
	const int neln = elem.Nodes();
	const int ndof = 6*neln;

	// get the elements material
	FEIncompressibleMaterial* pmi = dynamic_cast<FEIncompressibleMaterial*>(m_fem.GetMaterial(elem.GetMatID()));

	// average global derivatives
	vector<double> gradN(6*neln);
	gradN.zero();

	// initial element volume
	double Ve = 0;

	// global derivatives of shape functions
	double Nx, Ny, Nz, Mx, My, Mz;
	const double *gw = elem.GaussWeights();

	// jacobian
	double Ji[3][3], detJt, detJ0;

	double *Gr, *Gs, *H;

	// calculate the average thickness
	double* h0 = elem.m_h0, gt, za;

	// repeat over gauss-points
	for (n=0; n<nint; ++n)
	{
		// calculate jacobian
		detJ0 = elem.detJ0(n);
		detJt = elem.detJt(n);
		elem.invjact(Ji, n);

		detJt *= gw[n];

		Ve += detJ0*gw[n];

		Gr = elem.Hr(n);
		Gs = elem.Hs(n);
		H = elem.H(n);
		gt = elem.gt(n);

		// calculate global gradient of shape functions
		// note that we need the transposed of Ji, not Ji itself !
		for (i=0; i<neln; ++i)
		{
			za = 0.5*gt*h0[i];

			Nx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i];
			Ny = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i];
			Nz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i];

			Mx = za*Ji[0][0]*Gr[i] + za*Ji[1][0]*Gs[i] + Ji[2][0]*0.5*h0[i]*H[i];
			My = za*Ji[0][1]*Gr[i] + za*Ji[1][1]*Gs[i] + Ji[2][1]*0.5*h0[i]*H[i];
			Mz = za*Ji[0][2]*Gr[i] + za*Ji[1][2]*Gs[i] + Ji[2][2]*0.5*h0[i]*H[i];

			gradN[6*i  ] += Nx*detJt;
			gradN[6*i+1] += Ny*detJt;
			gradN[6*i+2] += Nz*detJt;

			gradN[6*i+3] += Mx*detJt;
			gradN[6*i+4] += My*detJt;
			gradN[6*i+5] += Mz*detJt;
		}
	}

	// get bulk modulus
	double k = pmi->m_K;

	// next, we add the Lagrangian contribution
	// note that this term will always be zero if the material does not
	// use the augmented lagrangian
	k += elem.m_Lk*pmi->hpp(elem.m_eJ);

	// divide by initial volume
	k /= Ve;

	for (i=0; i<ndof; ++i)
		for (j=i; j<ndof; ++j)
			ke[i][j] += k*gradN[i]*gradN[j];
}


//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FESolver::InternalForces(FESolidElement& el, vector<double>& fe)
{
	int i, n;

	// jacobian matrix, inverse jacobian matrix and determinants
	double Ji[3][3], detJt;

	double Gx, Gy, Gz;

	const double* Gr, *Gs, *Gt;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	double*	gw = el.GaussWeights();

	// repeat for all integration points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// calculate the jacobian
		el.invjact(Ji, n);
		detJt = el.detJt(n);

		detJt *= gw[n];

		// get the stress vector for this integration point
		mat3ds& s = pt.s;

		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);

		for (i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

			// calculate internal force
			// the '-' sign is so that the internal forces get subtracted
			// from the global residual vector
			fe[3*i  ] -= ( Gx*s.xx() +
				           Gy*s.xy() +
					       Gz*s.xz() )*detJt;

			fe[3*i+1] -= ( Gy*s.yy() +
				           Gx*s.xy() +
					       Gz*s.yz() )*detJt;

			fe[3*i+2] -= ( Gz*s.zz() +
				           Gy*s.yz() +
					       Gx*s.xz() )*detJt;
		}
	}
}


//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for enhanced strain
//! solid elements.

void FESolver::UDGInternalForces(FESolidElement& el, vector<double>& fe)
{
	// make sure this element is of the correct type
	assert(el.Type() == FE_UDGHEX);

	// get the stress data
	FEMaterialPoint& mp = *el.m_State[0];
	FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
	mat3ds& s = pt.s;

	// calculate the average cartesian derivatives
	double GX[8], GY[8], GZ[8];
	AvgCartDerivs(el, GX, GY, GZ);

	// calculate average deformation gradient Fbar
	mat3d Fb;
	AvgDefGrad(el, Fb, GX, GY, GZ);

	// calculate the transposed inverse of Fbar
	mat3d Fti = Fb.transinv();

	// calculate current element volume
	double ve = HexVolume(el, 1);

	// current averaged shape derivatives
	double Gx, Gy, Gz;

	// calculate the internal force
	for (int i=0; i<8; ++i)
	{
		Gx = Fti(0,0)*GX[i]+Fti(0,1)*GY[i]+Fti(0,2)*GZ[i];
		Gy = Fti(1,0)*GX[i]+Fti(1,1)*GY[i]+Fti(1,2)*GZ[i];
		Gz = Fti(2,0)*GX[i]+Fti(2,1)*GY[i]+Fti(2,2)*GZ[i];

		fe[3*i  ] -= ve*(Gx*s.xx() + Gy*s.xy() + Gz*s.xz());
		fe[3*i+1] -= ve*(Gx*s.xy() + Gy*s.yy() + Gz*s.yz());
		fe[3*i+2] -= ve*(Gx*s.xz() + Gy*s.yz() + Gz*s.zz());
	}

	// add hourglass forces
	UDGHourglassForces(el, fe);
}

//-----------------------------------------------------------------------------
//! calculates the hourglass forces

void FESolver::UDGHourglassForces(FESolidElement &el, vector<double> &fe)
{
	int i;

	const double h4[8] = { 1,-1, 1,-1, 1,-1, 1,-1 };
	const double h5[8] = { 1,-1,-1, 1,-1, 1, 1,-1 };
	const double h6[8] = { 1, 1,-1,-1,-1,-1, 1, 1 };
	const double h7[8] = {-1, 1,-1, 1, 1,-1, 1,-1 };

	vec3d* r0 = el.r0();
	vec3d* rt = el.rt();

	double x4 = 0, x5 = 0, x6 = 0, x7 = 0;
	double y4 = 0, y5 = 0, y6 = 0, y7 = 0;
	double z4 = 0, z5 = 0, z6 = 0, z7 = 0;

	double X4 = 0, X5 = 0, X6 = 0, X7 = 0;
	double Y4 = 0, Y5 = 0, Y6 = 0, Y7 = 0;
	double Z4 = 0, Z5 = 0, Z6 = 0, Z7 = 0;

	for (i=0; i<8; ++i)
	{
		X4 += h4[i]*r0[i].x; Y4 += h4[i]*r0[i].y; Z4 += h4[i]*r0[i].z;
		X5 += h5[i]*r0[i].x; Y5 += h5[i]*r0[i].y; Z5 += h5[i]*r0[i].z;
		X6 += h6[i]*r0[i].x; Y6 += h6[i]*r0[i].y; Z6 += h6[i]*r0[i].z;
		X7 += h7[i]*r0[i].x; Y7 += h7[i]*r0[i].y; Z7 += h7[i]*r0[i].z;

		x4 += h4[i]*rt[i].x; y4 += h4[i]*rt[i].y; z4 += h4[i]*rt[i].z;
		x5 += h5[i]*rt[i].x; y5 += h5[i]*rt[i].y; z5 += h5[i]*rt[i].z;
		x6 += h6[i]*rt[i].x; y6 += h6[i]*rt[i].y; z6 += h6[i]*rt[i].z;
		x7 += h7[i]*rt[i].x; y7 += h7[i]*rt[i].y; z7 += h7[i]*rt[i].z;
	}

	double GX[8], GY[8], GZ[8];
	AvgCartDerivs(el, GX, GY, GZ);

	mat3d F;
	AvgDefGrad(el, F, GX, GY, GZ);
/*
	double Ji0[3][3], Jt[3][3];
	el.invjac0(Ji0, 0);
	el.jact(Jt, 0);

	F[0][0] = Ji0[0][0]*Jt[0][0]+Ji0[0][1]*Jt[1][0]+Ji0[0][2]*Jt[2][0];
	F[0][1] = Ji0[0][0]*Jt[0][1]+Ji0[0][1]*Jt[1][1]+Ji0[0][2]*Jt[2][1];
	F[0][2] = Ji0[0][0]*Jt[0][2]+Ji0[0][1]*Jt[1][2]+Ji0[0][2]*Jt[2][2];

	F[1][0] = Ji0[1][0]*Jt[0][0]+Ji0[1][1]*Jt[1][0]+Ji0[1][2]*Jt[2][0];
	F[1][1] = Ji0[1][0]*Jt[0][1]+Ji0[1][1]*Jt[1][1]+Ji0[1][2]*Jt[2][1];
	F[1][2] = Ji0[1][0]*Jt[0][2]+Ji0[1][1]*Jt[1][2]+Ji0[1][2]*Jt[2][2];

	F[2][0] = Ji0[2][0]*Jt[0][0]+Ji0[2][1]*Jt[1][0]+Ji0[2][2]*Jt[2][0];
	F[2][1] = Ji0[2][0]*Jt[0][1]+Ji0[2][1]*Jt[1][1]+Ji0[2][2]*Jt[2][1];
	F[2][2] = Ji0[2][0]*Jt[0][2]+Ji0[2][1]*Jt[1][2]+Ji0[2][2]*Jt[2][2];
*/
	double u4 = 0, u5 = 0, u6 = 0, u7 = 0;
	double v4 = 0, v5 = 0, v6 = 0, v7 = 0;
	double w4 = 0, w5 = 0, w6 = 0, w7 = 0;

	u4 = x4 - (F[0][0]*X4 + F[0][1]*Y4 + F[0][2]*Z4);
	v4 = y4 - (F[1][0]*X4 + F[1][1]*Y4 + F[1][2]*Z4);
	w4 = z4 - (F[2][0]*X4 + F[2][1]*Y4 + F[2][2]*Z4);

	u5 = x5 - (F[0][0]*X5 + F[0][1]*Y5 + F[0][2]*Z5);
	v5 = y5 - (F[1][0]*X5 + F[1][1]*Y5 + F[1][2]*Z5);
	w5 = z5 - (F[2][0]*X5 + F[2][1]*Y5 + F[2][2]*Z5);

	u6 = x6 - (F[0][0]*X6 + F[0][1]*Y6 + F[0][2]*Z6);
	v6 = y6 - (F[1][0]*X6 + F[1][1]*Y6 + F[1][2]*Z6);
	w6 = z6 - (F[2][0]*X6 + F[2][1]*Y6 + F[2][2]*Z6);

	u7 = x7 - (F[0][0]*X7 + F[0][1]*Y7 + F[0][2]*Z7);
	v7 = y7 - (F[1][0]*X7 + F[1][1]*Y7 + F[1][2]*Z7);
	w7 = z7 - (F[2][0]*X7 + F[2][1]*Y7 + F[2][2]*Z7);

	double g4[8] = {0}, g5[8] = {0}, g6[8] = {0}, g7[8] = {0};

	for (i=0; i<8; ++i)
	{
		g4[i] = h4[i] - (GX[i]*X4 + GY[i]*Y4 + GZ[i]*Z4);
		g5[i] = h5[i] - (GX[i]*X5 + GY[i]*Y5 + GZ[i]*Z5);
		g6[i] = h6[i] - (GX[i]*X6 + GY[i]*Y6 + GZ[i]*Z6);
		g7[i] = h7[i] - (GX[i]*X7 + GY[i]*Y7 + GZ[i]*Z7);
	}

	// calculate hourglass forces
	double hg = m_fem.m_pStep->m_hg;

	for (i=0; i<8; ++i)
	{
		fe[3*i  ] -= hg*(g4[i]*u4 + g5[i]*u5 + g6[i]*u6 + g7[i]*u7);
		fe[3*i+1] -= hg*(g4[i]*v4 + g5[i]*v5 + g6[i]*v6 + g7[i]*v7);
		fe[3*i+2] -= hg*(g4[i]*w4 + g5[i]*w5 + g6[i]*w6 + g7[i]*w7);
	}
}

//-----------------------------------------------------------------------------
//! Calculates the average Cartesian derivatives
//! Note that we assume that the GX, GY and GX contain the averaged 
//! Cartesian derivatives

void FESolver::AvgDefGrad(FESolidElement& el, mat3d& F, double GX[8], double GY[8], double GZ[8])
{
	vec3d* rt = el.rt();

	F.zero();
	for (int i=0; i<8; ++i)
	{
		F[0][0] += rt[i].x*GX[i];
		F[0][1] += rt[i].x*GY[i];
		F[0][2] += rt[i].x*GZ[i];

		F[1][0] += rt[i].y*GX[i];
		F[1][1] += rt[i].y*GY[i];
		F[1][2] += rt[i].y*GZ[i];

		F[2][0] += rt[i].z*GX[i];
		F[2][1] += rt[i].z*GY[i];
		F[2][2] += rt[i].z*GZ[i];
	}
}

//-----------------------------------------------------------------------------
//! Calculates the average Cartesian derivatives

void FESolver::AvgCartDerivs(FESolidElement& el, double GX[8], double GY[8], double GZ[8], int nstate)
{
	// get the nodal coordinates
	vec3d* r = (nstate == 0? el.r0() : el.rt());
	double x1 = r[0].x, y1 = r[0].y, z1 = r[0].z;
	double x2 = r[1].x, y2 = r[1].y, z2 = r[1].z;
	double x3 = r[2].x, y3 = r[2].y, z3 = r[2].z;
	double x4 = r[3].x, y4 = r[3].y, z4 = r[3].z;
	double x5 = r[4].x, y5 = r[4].y, z5 = r[4].z;
	double x6 = r[5].x, y6 = r[5].y, z6 = r[5].z;
	double x7 = r[6].x, y7 = r[6].y, z7 = r[6].z;
	double x8 = r[7].x, y8 = r[7].y, z8 = r[7].z;

	const double f12 = 1.0/12.0;

	// set up the B-matrix
	// we use the G arrays to store the B-matrix
	GX[0] = f12*(y2*((z6-z3)-(z4-z5))+y3*(z2-z4)+y4*((z3-z8)-(z5-z2))+y5*((z8-z6)-(z2-z4))+y6*(z5-z2)+y8*(z4-z5));
	GY[0] = f12*(z2*((x6-x3)-(x4-x5))+z3*(x2-x4)+z4*((x3-x8)-(x5-x2))+z5*((x8-x6)-(x2-x4))+z6*(x5-x2)+z8*(x4-x5));
	GZ[0] = f12*(x2*((y6-y3)-(y4-y5))+x3*(y2-y4)+x4*((y3-y8)-(y5-y2))+x5*((y8-y6)-(y2-y4))+x6*(y5-y2)+x8*(y4-y5));

	GX[1] = f12*(y3*((z7-z4)-(z1-z6))+y4*(z3-z1)+y1*((z4-z5)-(z6-z3))+y6*((z5-z7)-(z3-z1))+y7*(z6-z3)+y5*(z1-z6));
	GY[1] = f12*(z3*((x7-x4)-(x1-x6))+z4*(x3-x1)+z1*((x4-x5)-(x6-x3))+z6*((x5-x7)-(x3-x1))+z7*(x6-x3)+z5*(x1-x6));
	GZ[1] = f12*(x3*((y7-y4)-(y1-y6))+x4*(y3-y1)+x1*((y4-y5)-(y6-y3))+x6*((y5-y7)-(y3-y1))+x7*(y6-y3)+x5*(y1-y6));

	GX[2] = f12*(y4*((z8-z1)-(z2-z7))+y1*(z4-z2)+y2*((z1-z6)-(z7-z4))+y7*((z6-z8)-(z4-z2))+y8*(z7-z4)+y6*(z2-z7));
	GY[2] = f12*(z4*((x8-x1)-(x2-x7))+z1*(x4-x2)+z2*((x1-x6)-(x7-x4))+z7*((x6-x8)-(x4-x2))+z8*(x7-x4)+z6*(x2-x7));
	GZ[2] = f12*(x4*((y8-y1)-(y2-y7))+x1*(y4-y2)+x2*((y1-y6)-(y7-y4))+x7*((y6-y8)-(y4-y2))+x8*(y7-y4)+x6*(y2-y7));

	GX[3] = f12*(y1*((z5-z2)-(z3-z8))+y2*(z1-z3)+y3*((z2-z7)-(z8-z1))+y8*((z7-z5)-(z1-z3))+y5*(z8-z1)+y7*(z3-z8));
	GY[3] = f12*(z1*((x5-x2)-(x3-x8))+z2*(x1-x3)+z3*((x2-x7)-(x8-x1))+z8*((x7-x5)-(x1-x3))+z5*(x8-x1)+z7*(x3-x8));
	GZ[3] = f12*(x1*((y5-y2)-(y3-y8))+x2*(y1-y3)+x3*((y2-y7)-(y8-y1))+x8*((y7-y5)-(y1-y3))+x5*(y8-y1)+x7*(y3-y8));

	GX[4] = f12*(y8*((z4-z7)-(z6-z1))+y7*(z8-z6)+y6*((z7-z2)-(z1-z8))+y1*((z2-z4)-(z8-z6))+y4*(z1-z8)+y2*(z6-z1));
	GY[4] = f12*(z8*((x4-x7)-(x6-x1))+z7*(x8-x6)+z6*((x7-x2)-(x1-x8))+z1*((x2-x4)-(x8-x6))+z4*(x1-x8)+z2*(x6-x1));
	GZ[4] = f12*(x8*((y4-y7)-(y6-y1))+x7*(y8-y6)+x6*((y7-y2)-(y1-y8))+x1*((y2-y4)-(y8-y6))+x4*(y1-y8)+x2*(y6-y1));

	GX[5] = f12*(y5*((z1-z8)-(z7-z2))+y8*(z5-z7)+y7*((z8-z3)-(z2-z5))+y2*((z3-z1)-(z5-z7))+y1*(z2-z5)+y3*(z7-z2));
	GY[5] = f12*(z5*((x1-x8)-(x7-x2))+z8*(x5-x7)+z7*((x8-x3)-(x2-x5))+z2*((x3-x1)-(x5-x7))+z1*(x2-x5)+z3*(x7-x2));
	GZ[5] = f12*(x5*((y1-y8)-(y7-y2))+x8*(y5-y7)+x7*((y8-y3)-(y2-y5))+x2*((y3-y1)-(y5-y7))+x1*(y2-y5)+x3*(y7-y2));

	GX[6] = f12*(y6*((z2-z5)-(z8-z3))+y5*(z6-z8)+y8*((z5-z4)-(z3-z6))+y3*((z4-z2)-(z6-z8))+y2*(z3-z6)+y4*(z8-z3));
	GY[6] = f12*(z6*((x2-x5)-(x8-x3))+z5*(x6-x8)+z8*((x5-x4)-(x3-x6))+z3*((x4-x2)-(x6-x8))+z2*(x3-x6)+z4*(x8-x3));
	GZ[6] = f12*(x6*((y2-y5)-(y8-y3))+x5*(y6-y8)+x8*((y5-y4)-(y3-y6))+x3*((y4-y2)-(y6-y8))+x2*(y3-y6)+x4*(y8-y3));

	GX[7] = f12*(y7*((z3-z6)-(z5-z4))+y6*(z7-z5)+y5*((z6-z1)-(z4-z7))+y4*((z1-z3)-(z7-z5))+y3*(z4-z7)+y1*(z5-z4));
	GY[7] = f12*(z7*((x3-x6)-(x5-x4))+z6*(x7-x5)+z5*((x6-x1)-(x4-x7))+z4*((x1-x3)-(x7-x5))+z3*(x4-x7)+z1*(x5-x4));
	GZ[7] = f12*(x7*((y3-y6)-(y5-y4))+x6*(y7-y5)+x5*((y6-y1)-(y4-y7))+x4*((y1-y3)-(y7-y5))+x3*(y4-y7)+x1*(y5-y4));

	// calculate the volume
	double Vi = 1./(x1*GX[0]+x2*GX[1]+x3*GX[2]+x4*GX[3]+x5*GX[4]+x6*GX[5]+x7*GX[6]+x8*GX[7]);

	// divide the B-matrix by the volume
	GX[0] *= Vi; GY[0] *= Vi; GZ[0] *= Vi;
	GX[1] *= Vi; GY[1] *= Vi; GZ[1] *= Vi;
	GX[2] *= Vi; GY[2] *= Vi; GZ[2] *= Vi;
	GX[3] *= Vi; GY[3] *= Vi; GZ[3] *= Vi;
	GX[4] *= Vi; GY[4] *= Vi; GZ[4] *= Vi;
	GX[5] *= Vi; GY[5] *= Vi; GZ[5] *= Vi;
	GX[6] *= Vi; GY[6] *= Vi; GZ[6] *= Vi;
	GX[7] *= Vi; GY[7] *= Vi; GZ[7] *= Vi;
}

//-----------------------------------------------------------------------------
//! Calculates the exact volume of a hexahedral element

double FESolver::HexVolume(FESolidElement& el, int state)
{
	// let's make sure this is indeed a hex element
//	assert(el.Type() == FE_HEX);

	// get the nodal coordinates
	vec3d* r = (state == 0? el.r0() : el.rt());
	double x1 = r[0].x, y1 = r[0].y, z1 = r[0].z;
	double x2 = r[1].x, y2 = r[1].y, z2 = r[1].z;
	double x3 = r[2].x, y3 = r[2].y, z3 = r[2].z;
	double x4 = r[3].x, y4 = r[3].y, z4 = r[3].z;
	double x5 = r[4].x, y5 = r[4].y, z5 = r[4].z;
	double x6 = r[5].x, y6 = r[5].y, z6 = r[5].z;
	double x7 = r[6].x, y7 = r[6].y, z7 = r[6].z;
	double x8 = r[7].x, y8 = r[7].y, z8 = r[7].z;

	// set up the B-matrix
	double B1[8];
//	double B2[8];
//	double B3[8];

	const double f12 = 1.0/12.0;

	B1[0] = f12*(y2*((z6-z3)-(z4-z5))+y3*(z2-z4)+y4*((z3-z8)-(z5-z2))+y5*((z8-z6)-(z2-z4))+y6*(z5-z2)+y8*(z4-z5));
//	B2[0] = f12*(z2*((x6-x3)-(x4-x5))+z3*(x2-x4)+z4*((x3-x8)-(x5-x2))+z5*((x8-x6)-(x2-x4))+z6*(x5-x2)+z8*(x4-x5));
//	B3[0] = f12*(x2*((y6-y3)-(y4-y5))+x3*(y2-y4)+x4*((y3-y8)-(y5-y2))+x5*((y8-y6)-(y2-y4))+x6*(y5-y2)+x8*(y4-y5));

	B1[1] = f12*(y3*((z7-z4)-(z1-z6))+y4*(z3-z1)+y1*((z4-z5)-(z6-z3))+y6*((z5-z7)-(z3-z1))+y7*(z6-z3)+y5*(z1-z6));
//	B2[1] = f12*(z3*((x7-x4)-(x1-x6))+z4*(x3-x1)+z1*((x4-x5)-(x6-x3))+z6*((x5-x7)-(x3-x1))+z7*(x6-x3)+z5*(x1-x6));
//	B3[1] = f12*(x3*((y7-y4)-(y1-y6))+x4*(y3-y1)+x1*((y4-y5)-(y6-y3))+x6*((y5-y7)-(y3-y1))+x7*(y6-y3)+x5*(y1-y6));

	B1[2] = f12*(y4*((z8-z1)-(z2-z7))+y1*(z4-z2)+y2*((z1-z6)-(z7-z4))+y7*((z6-z8)-(z4-z2))+y8*(z7-z4)+y6*(z2-z7));
//	B2[2] = f12*(z4*((x8-x1)-(x2-x7))+z1*(x4-x2)+z2*((x1-x6)-(x7-x4))+z7*((x6-x8)-(x4-x2))+z8*(x7-x4)+z6*(x2-x7));
//	B3[2] = f12*(x4*((y8-y1)-(y2-y7))+x1*(y4-y2)+x2*((y1-y6)-(y7-y4))+x7*((y6-y8)-(y4-y2))+x8*(y7-y4)+x6*(y2-y7));

	B1[3] = f12*(y1*((z5-z2)-(z3-z8))+y2*(z1-z3)+y3*((z2-z7)-(z8-z1))+y8*((z7-z5)-(z1-z3))+y5*(z8-z1)+y7*(z3-z8));
//	B2[3] = f12*(z1*((x5-x2)-(x3-x8))+z2*(x1-x3)+z3*((x2-x7)-(x8-x1))+z8*((x7-x5)-(x1-x3))+z5*(x8-x1)+z7*(x3-x8));
//	B3[3] = f12*(x1*((y5-y2)-(y3-y8))+x2*(y1-y3)+x3*((y2-y7)-(y8-y1))+x8*((y7-y5)-(y1-y3))+x5*(y8-y1)+x7*(y3-y8));

	B1[4] = f12*(y8*((z4-z7)-(z6-z1))+y7*(z8-z6)+y6*((z7-z2)-(z1-z8))+y1*((z2-z4)-(z8-z6))+y4*(z1-z8)+y2*(z6-z1));
//	B2[4] = f12*(z8*((x4-x7)-(x6-x1))+z7*(x8-x6)+z6*((x7-x2)-(x1-x8))+z1*((x2-x4)-(x8-x6))+z4*(x1-x8)+z2*(x6-x1));
//	B3[4] = f12*(x8*((y4-y7)-(y6-y1))+x7*(y8-y6)+x6*((y7-y2)-(y1-y8))+x1*((y2-y4)-(y8-y6))+x4*(y1-y8)+x2*(y6-y1));

	B1[5] = f12*(y5*((z1-z8)-(z7-z2))+y8*(z5-z7)+y7*((z8-z3)-(z2-z5))+y2*((z3-z1)-(z5-z7))+y1*(z2-z5)+y3*(z7-z2));
//	B2[5] = f12*(z5*((x1-x8)-(x7-x2))+z8*(x5-x7)+z7*((x8-x3)-(x2-x5))+z2*((x3-x1)-(x5-x7))+z1*(x2-x5)+z3*(x7-x2));
//	B3[5] = f12*(x5*((y1-y8)-(y7-y2))+x8*(y5-y7)+x7*((y8-y3)-(y2-y5))+x2*((y3-y1)-(y5-y7))+x1*(y2-y5)+x3*(y7-y2));

	B1[6] = f12*(y6*((z2-z5)-(z8-z3))+y5*(z6-z8)+y8*((z5-z4)-(z3-z6))+y3*((z4-z2)-(z6-z8))+y2*(z3-z6)+y4*(z8-z3));
//	B2[6] = f12*(z6*((x2-x5)-(x8-x3))+z5*(x6-x8)+z8*((x5-x4)-(x3-x6))+z3*((x4-x2)-(x6-x8))+z2*(x3-x6)+z4*(x8-x3));
//	B3[6] = f12*(x6*((y2-y5)-(y8-y3))+x5*(y6-y8)+x8*((y5-y4)-(y3-y6))+x3*((y4-y2)-(y6-y8))+x2*(y3-y6)+x4*(y8-y3));

	B1[7] = f12*(y7*((z3-z6)-(z5-z4))+y6*(z7-z5)+y5*((z6-z1)-(z4-z7))+y4*((z1-z3)-(z7-z5))+y3*(z4-z7)+y1*(z5-z4));
//	B2[7] = f12*(z7*((x3-x6)-(x5-x4))+z6*(x7-x5)+z5*((x6-x1)-(x4-x7))+z4*((x1-x3)-(x7-x5))+z3*(x4-x7)+z1*(x5-x4));
//	B3[7] = f12*(x7*((y3-y6)-(y5-y4))+x6*(y7-y5)+x5*((y6-y1)-(y4-y7))+x4*((y1-y3)-(y7-y5))+x3*(y4-y7)+x1*(y5-y4));

	// calculate the volume V= xi*B1[i] = yi*B2[i] = zi*B3[i] (sum over i)
	return (x1*B1[0]+x2*B1[1]+x3*B1[2]+x4*B1[3]+x5*B1[4]+x6*B1[5]+x7*B1[6]+x8*B1[7]);
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal forces due to hydrostatic pressure

bool FESolver::LinearPressureForce(FESurfaceElement& el, vector<double>& fe)
{
	int i, n;

	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// pressure at nodes
	double *pn = el.pt();

	// nodal coordinates
	vec3d *r0 = el.r0();

	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();

	// pressure at integration points
	double pr;

	vec3d dxr, dxs;

	// force vector
	vec3d f;

	// repeat over integration points
	fe.zero();
	for (n=0; n<nint; ++n)
	{
		N  = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		pr = 0;
		dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			pr += N[i]*pn[i];

			dxr.x += Gr[i]*r0[i].x;
			dxr.y += Gr[i]*r0[i].y;
			dxr.z += Gr[i]*r0[i].z;

			dxs.x += Gs[i]*r0[i].x;
			dxs.y += Gs[i]*r0[i].y;
			dxs.z += Gs[i]*r0[i].z;
		}

		f = (dxr ^ dxs)*pr*w[n];

		for (i=0; i<neln; ++i)
		{
			fe[3*i  ] += N[i]*f.x;
			fe[3*i+1] += N[i]*f.y;
			fe[3*i+2] += N[i]*f.z;
		}
	}

	return true;
}


//-----------------------------------------------------------------------------
//! calculates the equivalent nodal forces due to hydrostatic pressure

bool FESolver::PressureForce(FESurfaceElement& el, vector<double>& fe)
{
	int i, n;

	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// pressure at nodes
	double *pn = el.pt();

	// nodal coordinates
	vec3d *rt = el.rt();

	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();

	// pressure at integration points
	double pr;

	vec3d dxr, dxs;

	// force vector
	vec3d f;

	// repeat over integration points
	fe.zero();
	for (n=0; n<nint; ++n)
	{
		N  = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		pr = 0;
		dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			pr += N[i]*pn[i];

			dxr.x += Gr[i]*rt[i].x;
			dxr.y += Gr[i]*rt[i].y;
			dxr.z += Gr[i]*rt[i].z;

			dxs.x += Gs[i]*rt[i].x;
			dxs.y += Gs[i]*rt[i].y;
			dxs.z += Gs[i]*rt[i].z;
		}

		f = (dxr ^ dxs)*pr*w[n];

		for (i=0; i<neln; ++i)
		{
			fe[3*i  ] += N[i]*f.x;
			fe[3*i+1] += N[i]*f.y;
			fe[3*i+2] += N[i]*f.z;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure

bool FESolver::PressureStiffness(FESurfaceElement& el, matrix& ke)
{
	int i, j, n;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// pressure at nodes
	double *pn = el.pt();

	// pressure at integration point
	double p;

	// gauss weights
	double* w = el.GaussWeights();

	// nodal coordinates
	vec3d* rt = el.rt();

	// jacobian
	double J[3][2];

	double t1, t2;
	double kab[3];

	ke.zero();

	double* N, *Gr, *Gs;

	// repeat over integration points
	for (n=0; n<nint; ++n)
	{
		N = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		// calculate pressure at integration point
		p = 0;
		for (i=0; i<neln; ++i) p += N[i]*pn[i];

		// calculate jacobian
		J[0][0] = J[0][1] = 0;
		J[1][0] = J[1][1] = 0;
		J[2][0] = J[2][1] = 0;
		for (i=0; i<neln; ++i)
		{
			J[0][0] += Gr[i]*rt[i].x; J[0][1] += Gs[i]*rt[i].x;
			J[1][0] += Gr[i]*rt[i].y; J[1][1] += Gs[i]*rt[i].y;
			J[2][0] += Gr[i]*rt[i].z; J[2][1] += Gs[i]*rt[i].z;
		}

		// calculate stiffness component
		for (i=0; i<neln; ++i)
			for (j=0; j<neln; ++j)
			{
				t1 = 0.5*(Gs[i]*N[j] - Gs[j]*N[i]);
				t2 = 0.5*(Gr[i]*N[j] - Gr[j]*N[i]);

				kab[0] = p*(J[0][0]*t1 - J[0][1]*t2)*w[n];
				kab[1] = p*(J[1][0]*t1 - J[1][1]*t2)*w[n];
				kab[2] = p*(J[2][0]*t1 - J[2][1]*t2)*w[n];

				ke[3*i  ][3*j  ] +=       0; //(0,0,0)*kab[0]+(0,0,1)*kab[1]+(0,0,2)*kab[2];
				ke[3*i  ][3*j+1] +=  kab[2]; //(0,1,0)*kab[0]+(0,1,1)*kab[1]+(0,1,2)*kab[2];
				ke[3*i  ][3*j+2] += -kab[1]; //(0,2,0)*kab[0]+(0,2,1)*kab[1]+(0,2,2)*kab[2];

				ke[3*i+1][3*j  ] += -kab[2]; //(1,0,0)*kab[0]+(1,0,1)*kab[1]+(1,0,2)*kab[2];
				ke[3*i+1][3*j+1] +=       0; //(1,1,0)*kab[0]+(1,1,1)*kab[1]+(1,1,2)*kab[2];
				ke[3*i+1][3*j+2] +=  kab[0]; //(1,2,0)*kab[0]+(1,2,1)*kab[1]+(1,2,2)*kab[2];

				ke[3*i+2][3*j  ] +=  kab[1]; //(2,0,0)*kab[0]+(2,0,1)*kab[1]+(2,0,2)*kab[2];
				ke[3*i+2][3*j+1] += -kab[0]; //(2,1,0)*kab[0]+(2,1,1)*kab[1]+(2,1,2)*kab[2];
				ke[3*i+2][3*j+2] +=       0; //(2,2,0)*kab[0]+(2,2,1)*kab[1]+(2,2,2)*kab[2];
			}
	}

	return true;
}


//-----------------------------------------------------------------------------
//! calculates the body forces

void FESolver::BodyForces(FESolidElement& el, vector<double>& fe)
{
	int i, n;
	double *H;

	// jacobian
	double detJ, dens;

	FEMaterial* pme = m_fem.GetMaterial(el.GetMatID());

	dens = pme->Density();

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	double* gw = el.GaussWeights();

	// loop over integration points
	vec3d g = m_fem.m_acc*dens;
	for (n=0; n<nint; ++n)
	{
		detJ = el.detJ0(n)*gw[n];

		H = el.H(n);

		for (i=0; i<neln; ++i)
		{
			fe[3*i  ] += H[i]*g.x*detJ;
			fe[3*i+1] += H[i]*g.y*detJ;
			fe[3*i+2] += H[i]*g.z*detJ;
		}						
	}
}

//-----------------------------------------------------------------------------
//! calculates the concentrated nodal forces

void FESolver::NodalForces(vector<double>& F)
{
	int i, id, bc, lc, n;
	double s, f;
	vec3d a;
	int* lm;

	// zero nodal force vector
	F.zero();

	FEMesh& mesh = m_fem.m_mesh;

	// loop over nodal force cards
	int ncnf = m_fem.m_FC.size();
	FENodalForce* FC = m_fem.m_FC;
	for (i=0; i<ncnf; ++i)
	{
		if (FC[i].IsActive())
		{
			id	 = FC[i].node;	// node ID
			bc   = FC[i].bc;	// direction of force
			lc   = FC[i].lc;	// loadcurve number
			s    = FC[i].s;		// force scale factor

			FENode& node = mesh.Node(id);

			n = node.m_ID[bc];
		
			f = s*m_fem.GetLoadCurve(lc)->Value();
			if (bc == 6) f *= m_fem.m_pStep->m_dt;	// GAA

			if (n >= 0) F[n] = f;
			else if (node.m_rid >=0)
			{
				// this is a rigid body node
				FERigidBody& RB = m_fem.m_RB[node.m_rid];

				// get the relative position
				a = node.m_rt - RB.m_rt;

				lm = RB.m_LM;
				switch (bc)
				{
				case 0:
					if (lm[0] >= 0) F[lm[0]] +=  f;
					if (lm[4] >= 0) F[lm[4]] +=  a.z*f;
					if (lm[5] >= 0) F[lm[5]] += -a.y*f;
					break;
				case 1:
					if (lm[1] >= 0) F[lm[1]] +=  f;
					if (lm[3] >= 0) F[lm[3]] += -a.z*f;
					if (lm[5] >= 0) F[lm[5]] +=  a.x*f;
					break;
				case 2:
					if (lm[2] >= 0) F[lm[2]] +=  f;
					if (lm[3] >= 0) F[lm[3]] +=  a.y*f;
					if (lm[4] >= 0) F[lm[4]] += -a.x*f;
					break;
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces due to the fluid work
//! Note that we only use the first n entries in fe, where n is the number
//! of nodes

bool FESolver::InternalFluidWork(FESolidElement& el, vector<double>& fe)
{
	int i, n;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// jacobian
	double Ji[3][3], detJ;

	double *Gr, *Gs, *Gt, *H;
	double Gx, Gy, Gz;

	vec3d* vt = el.vt();
	double* pn = el.pt();

	double trd;	// trace of velocity gradient

	// Bp-matrix
	vector<double[3]> B(neln);

	// gauss-weights
	double* wg = el.GaussWeights();

	// get the element's material
	FEPoroElastic* pm = dynamic_cast<FEPoroElastic*> (m_fem.GetMaterial(el.GetMatID()));
	if (pm == 0)
	{
		m_log.printbox("FATAL ERROR", "Incorrect material type\n");
		return false;
	}

	fe.zero();

	double dt = m_fem.m_pStep->m_dt;

	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEPoroElasticMaterialPoint& pt = *(el.m_State[n]->ExtractData<FEPoroElasticMaterialPoint>());

		// calculate jacobian
		el.invjact(Ji, n);
		detJ = el.detJt(n);

		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);

		H = el.H(n);

		trd = 0;

		for (i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

			// calculate trace of velocity gradient
			trd += Gx*vt[i].x + Gy*vt[i].y + Gz*vt[i].z;

			// calculate Bp matrix
			B[i][0] = Gx;
			B[i][1] = Gy;
			B[i][2] = Gz;
		}

		// get the flux
		vec3d& w = pt.m_w;

		// update force vector
		for (i=0; i<neln; ++i)
		{
			fe[i] -= dt*(B[i][0]*w.x+B[i][1]*w.y+B[i][2]*w.z - trd*H[i])*detJ*wg[n];
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! calculates element stiffness matrix for element iel

bool FESolver::ElementPoroStiffness(FESolidElement& el, matrix& ke)
{
	int i, j, l, n;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	double *Gr, *Gs, *Gt, *H;
	double Gx, Gy, Gz;

	// jacobian
	double Ji[3][3], detJ;

	// Bp-matrix
	vector<double[3]> B(neln);

	// permeability tensor
	double k[3][3];

	// zero stiffness matrix
	ke.zero();

	// calculate solid stiffness matrix
	int ndof = 3*el.Nodes();
	matrix ks(ndof, ndof); ks.zero();
	ElementStiffness(el, ks);

	// copy solid stiffness matrix into ke
	for (i=0; i<3*neln; ++i)
		for (j=0; j<3*neln; ++j)
		{
			ke[i][j] = ks[i][j];
		}

	// get the element's material
	FEPoroElastic* pm = dynamic_cast<FEPoroElastic*> (m_fem.GetMaterial(el.GetMatID()));
	if (pm == 0)
	{
		m_log.printbox("FATAL ERROR", "Incorrect material type\n");
		return false;
	}

	// gauss-weights
	double* gw = el.GaussWeights();

	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// calculate jacobian
		el.invjact(Ji, n);
		detJ = el.detJt(n);

		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);

		H = el.H(n);

		for (i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

			// calculate Bp matrix
			B[i][0] = Gx;
			B[i][1] = Gy;
			B[i][2] = Gz;
		}

		// get the permeability tensor
		// TODO: note that we pass a material point although it is not initialized
		// we can get away with this for now since the permeability tensor that
		// is returned is constant, but if this is no longer the case we should
		// initialize the material point data.
		pm->Permeability(k, mp);

		// calculate the Q = Bt*k*B matrix
		double dt = m_fem.m_pStep->m_dt;
		for (i=0; i<neln; ++i)
			for (j=0; j<neln; ++j)
			{
				for (l=0; l<3; ++l)
				{
					ke[3*neln+i][3*neln+j] -= dt*detJ*gw[n]*(B[i][0]*k[0][l]+B[i][1]*k[1][l]+B[i][2]*k[2][l])*B[j][l];
				}
			}

		// calculate the G-matrix
		for (i=0; i<neln; ++i)
			for (j=0; j<neln; ++j)
			{
				ke[3*i  ][3*neln+j] -= detJ*gw[n]*B[i][0]*H[j];
				ke[3*i+1][3*neln+j] -= detJ*gw[n]*B[i][1]*H[j];
				ke[3*i+2][3*neln+j] -= detJ*gw[n]*B[i][2]*H[j];

				ke[3*neln+j][3*i  ] -= detJ*gw[n]*B[i][0]*H[j];
				ke[3*neln+j][3*i+1] -= detJ*gw[n]*B[i][1]*H[j];
				ke[3*neln+j][3*i+2] -= detJ*gw[n]*B[i][2]*H[j];
			}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! This function calculates the inertial forces for dynamic problems

void FESolver::InertialForces(vector<double>& R)
{
	int i, j, iel, n;
	int nint, neln;
	double *H, kab;

	// get the mesh
	FEMesh& mesh = m_fem.m_mesh;

	// allocate F
	vector<double> F(3*mesh.Nodes());
	F.zero();

	vector<double> fe;

	// calculate F
	double a = 4.0 / m_fem.m_pStep->m_dt;
	double b = a / m_fem.m_pStep->m_dt;
	for (i=0; i<mesh.Nodes(); ++i)
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
	// first do the solid elements
	matrix ke;
	for (iel=0; iel<mesh.SolidElements(); ++iel)
	{
		FESolidElement& el = mesh.SolidElement(iel);
		mesh.UnpackElement(el);

		FEMaterial* pme = m_fem.GetMaterial(el.GetMatID());

		double d = pme->Density();

		nint = el.GaussPoints();
		neln = el.Nodes();

		ke.Create(3*neln, 3*neln);
		ke.zero();

		fe.create(3*neln);
		
		// create the element mass matrix
		for (n=0; n<nint; ++n)
		{
			double detJ0 = el.detJ0(n)*el.GaussWeights()[n];

			H = el.H(n);
			for (i=0; i<neln; ++i)
				for (j=0; j<neln; ++j)
				{
					kab = H[i]*H[j]*detJ0*d;
					ke[3*i  ][3*j  ] += kab;
					ke[3*i+1][3*j+1] += kab;
					ke[3*i+2][3*j+2] += kab;
				}	
		}

		// now, multiply M with F and add to R
		int* en = el.m_node;
		for (i=0; i<3*neln; ++i)
		{
			fe[i] = 0;
			for (j=0; j<3*neln; ++j)
			{
				fe[i] -= ke[i][j]*F[3*(en[j/3]) + j%3];
			}
		}

		// assemble fe into R
		AssembleResidual(el.m_node, el.LM(), fe, R);
	}

	// TODO: do dynamics for shell elements
}

//-----------------------------------------------------------------------------
//! Calculates the forces due to discrete elements (i.e. springs)

void FESolver::DiscreteElementForces(vector<double>& R)
{
	if (m_fem.m_DE.size() == 0) return;

	FEMesh& mesh = m_fem.m_mesh;

	vector<double> fe(6);
	vec3d u1, u2;
	double E;

	vector<int> en(2), lm(6);

	for (int i=0; i<m_fem.m_DE.size(); ++i)
	{
		FE_DISCRETE_ELEMENT& el = m_fem.m_DE[i];
		E = el.E;

		FENode& n1 = mesh.Node(el.n1);
		FENode& n2 = mesh.Node(el.n2);

		vec3d& r01 = n1.m_r0;
		vec3d& r02 = n2.m_r0;
		vec3d& rt1 = n1.m_rt;
		vec3d& rt2 = n2.m_rt;

		u1 = rt1 - r01;
		u2 = rt2 - r02;

		fe[0] = -E*(u1.x - u2.x);
		fe[1] = -E*(u1.y - u2.y);
		fe[2] = -E*(u1.z - u2.z);
		fe[3] = -E*(u2.x - u1.x);
		fe[4] = -E*(u2.y - u1.y);
		fe[5] = -E*(u2.z - u1.z);

		en[0] = el.n1;
		en[1] = el.n2;

		lm[0] = n1.m_ID[0];
		lm[1] = n1.m_ID[1];
		lm[2] = n1.m_ID[2];
		lm[3] = n2.m_ID[0];
		lm[4] = n2.m_ID[1];
		lm[5] = n2.m_ID[2];

		AssembleResidual(en, lm, fe, R);
	}
}

//-----------------------------------------------------------------------------
//! Calculates the discrete element stiffness

void FESolver::DiscreteElementStiffness()
{
	if (m_fem.m_DE.size() == 0) return;

	FEMesh& mesh = m_fem.m_mesh;

	double E;
	matrix ke(6,6);
	ke.zero();

	vector<int> en(2), lm(6);

	for (int i=0; i<m_fem.m_DE.size(); ++i)
	{
		FE_DISCRETE_ELEMENT& el = m_fem.m_DE[i];
		E = el.E;

		FENode& n1 = mesh.Node(el.n1);
		FENode& n2 = mesh.Node(el.n2);

		ke[0][0] = ke[1][1] = ke[2][2] =  E;
		ke[0][3] = ke[1][4] = ke[2][5] =  0;
		ke[3][3] = ke[4][4] = ke[5][5] =  E;
		ke[3][0] = ke[4][1] = ke[5][2] =  0;

		en[0] = el.n1;
		en[1] = el.n2;

		lm[0] = n1.m_ID[0];
		lm[1] = n1.m_ID[1];
		lm[2] = n1.m_ID[2];
		lm[3] = n2.m_ID[0];
		lm[4] = n2.m_ID[1];
		lm[5] = n2.m_ID[2];

		AssembleStiffness(en, lm, ke);
	}
}
