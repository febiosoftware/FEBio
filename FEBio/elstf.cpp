#include "stdafx.h"
#include <math.h>
#include "FESolidSolver.h"
#include "FEPoroElastic.h"
#include "FEMicroMaterial.h"
#include "tens4d.h"
#include "log.h"

//-----------------------------------------------------------------------------
//! Calculates global stiffness matrix.

bool FESolidSolver::StiffnessMatrix()
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

	// total nr of elements
	int TNE = mesh.Elements();

	// repeat over all solid elements
	FESolidDomain& bd = mesh.SolidDomain();
	int NE = bd.size();
	for (iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = mesh.SolidElement(iel);
		if (!el.IsRigid())
		{
			bd.UnpackElement(el);

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
				bd.ElementStiffness(m_fem, el, ke);

				// add the inertial stiffness for dynamics
				if (m_fem.m_pStep->m_nanalysis == FE_DYNAMIC) bd.ElementInertialStiffness(m_fem, el, ke);

				// assemble element matrix in global stiffness matrix
				AssembleStiffness(el.m_node, el.LM(), ke);
			}
			else if (dynamic_cast<FEPoroElastic*>(pmat))
			{
				// allocate stiffness matrix
				int neln = el.Nodes();
				ndof = neln*4;
				ke.Create(ndof, ndof);
		
				// calculate the element stiffness matrix
				bd.ElementPoroStiffness(m_fem, el, ke);

				// TODO: the problem here is that the LM array that is returned by the UnpackElement
				// function does not give the equation numbers in the right order. For this reason we
				// have to create a new lm array and place the equation numbers in the right order.
				// What we really ought to do is fix the UnpackElement function so that it returns
				// the LM vector in the right order for poroelastic elements.
				vector<int> lm(ndof);
				for (int i=0; i<neln; ++i)
				{
					lm[4*i  ] = el.LM()[3*i];
					lm[4*i+1] = el.LM()[3*i+1];
					lm[4*i+2] = el.LM()[3*i+2];
					lm[4*i+3] = el.LM()[3*neln+i];
				}
				
				// assemble element matrix in global stiffness matrix
				AssembleStiffness(el.m_node, lm, ke);
			}
		}
		else
		{
			// for dynamic analyses we do need to add the inertial stiffness of the rigid body
			if (m_fem.m_pStep->m_nanalysis == FE_DYNAMIC)
			{
				bd.UnpackElement(el);

				ndof = 3*el.Nodes();
				ke.Create(ndof, ndof);
				ke.zero();

				// add the inertial stiffness for dynamics
				bd.ElementInertialStiffness(m_fem, el, ke);

				// assemble element matrix in global stiffness matrix
				AssembleStiffness(el.m_node, el.LM(), ke);
			}
		}

		if (m_fem.m_pStep->GetPrintLevel() == FE_PRINT_MINOR_ITRS_EXP)
		{
			fprintf(stderr, "Calculating stiffness matrix: %.1lf %% \r", 100.0*iel/ TNE);
		}
	}

	// repeat over all shell elements
	FEShellDomain& sd = mesh.ShellDomain();
	int NS = sd.size();
	for (iel=0; iel<NS; ++iel)
	{
		FEShellElement& el = mesh.ShellElement(iel);
		if (!el.IsRigid())
		{
			sd.UnpackElement(el);

			// get the elements material
			FEMaterial* pmat = m_fem.GetMaterial(el.GetMatID());

			// skip rigid elements and poro-elastic elements
			if (dynamic_cast<FEPoroElastic*>(pmat) == 0)
			{
				// create the element's stiffness matrix
				ndof = 6*el.Nodes();
				ke.Create(ndof, ndof);

				// calculate the element stiffness matrix
				sd.ElementStiffness(m_fem, el, ke);

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
			fprintf(stderr, "Calculating stiffness matrix: %.1lf %% \r", 100.0*(NE + iel)/ TNE);
		}
	}

	// repeat over truss elements
	FETrussDomain& td = mesh.TrussDomain();
	int NT = td.size();
	for (iel =0; iel<NT; ++iel)
	{
		FETrussElement& el = td.Element(iel);
		td.UnpackElement(el);
		td.ElementStiffness(m_fem, el, ke);
		AssembleStiffness(el.m_node, el.LM(), ke);
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
			FEPressureLoad& pc = m_fem.m_PC[m];
			if (pc.bc == 0)
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
			if (K.diag(i) == 0) 
			{
				throw ZeroDiagonal(i, m_fem);
			}
		}
	}

	return true;
}

//-----------------------------------------------------------------------------

void FESolidSolver::LinearConstraintStiffness()
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

void FESolidSolver::ContactStiffness()
{
	for (int i=0; i<m_fem.ContactInterfaces(); ++i) m_fem.m_CI[i].ContactStiffness();
}

//-----------------------------------------------------------------------------
//! This function calculates the rigid stiffness matrices

void FESolidSolver::RigidStiffness(vector<int>& en, vector<int>& elm, matrix& ke)
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

void FESolidSolver::AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke)
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
void FESolidSolver::ContactForces(vector<double>& R)
{
	for (int i=0; i<m_fem.ContactInterfaces(); ++i) m_fem.m_CI[i].ContactForces(R);
}

//-----------------------------------------------------------------------------
//! calculates the residual vector
//! Note that the concentrated nodal forces are not calculated here.
//! This is because they do not depend on the geometry 
//! so we only calculate them once (in Quasin) and then add them here.

bool FESolidSolver::Residual(vector<double>& R)
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
	FESolidDomain& bd = mesh.SolidDomain();
	int NE = bd.size();
	for (i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = bd.Element(i);

		// unpack the element
		if (!el.IsRigid()) bd.UnpackElement(el);

		FEMaterial* pm = m_fem.GetMaterial(el.GetMatID());

		// get the element force vector and initialize it to zero
		ndof = 3*el.Nodes();
		fe.assign(ndof, 0);

		// skip rigid elements for internal force calculations
		if (!el.IsRigid())
		{
			// calculate internal force vector
			if (el.Type() != FE_UDGHEX) bd.InternalForces(el, fe);
			else bd.UDGInternalForces(m_fem, el, fe);

			// apply body forces
			if (m_fem.UseBodyForces())
			{
				bd.BodyForces(m_fem, el, fe);
			}

			// assemble element 'fe'-vector into global R vector
			AssembleResidual(el.m_node, el.LM(), fe, R);
		}
		else if (m_fem.UseBodyForces())
		{
			// unpack the element
			bd.UnpackElement(el);

			// apply body force to rigid elements
			bd.BodyForces(m_fem, el, fe);

			// assemble element 'fe'-vector into global R vector
			AssembleResidual(el.m_node, el.LM(), fe, R);
		}

		// do poro-elastic forces
		if ((m_fem.m_pStep->m_nModule == FE_POROELASTIC)&&(dynamic_cast<FEPoroElastic*>(pm)))
		{
			// calculate fluid internal work
			bd.InternalFluidWork(m_fem, el, fe);

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
	FEShellDomain& sd = mesh.ShellDomain();
	for (i=0; i<NS; ++i)
	{
		// get the element
		FEShellElement& el = mesh.ShellElement(i);

		// create the element force vector and initialize to zero
		ndof = 6*el.Nodes();
		fe.assign(ndof, 0);

		if (!el.IsRigid())
		{
			sd.UnpackElement(el);

			// skip rigid elements for internal force calculation
			sd.InternalForces(el, fe);

			// apply body forces to shells
			if (m_fem.UseBodyForces())
			{
				sd.BodyForces(m_fem, el, fe);
			}

			// assemble the residual
			AssembleResidual(el.m_node, el.LM(), fe, R);
		}
		else if (m_fem.UseBodyForces())
		{
			sd.UnpackElement(el);

			// apply body forces to shells
			sd.BodyForces(m_fem, el, fe);

			// assemble the residual
			AssembleResidual(el.m_node, el.LM(), fe, R);
		}

		// TODO: Do poro-elasticity for shells
	}

	// loop over truss elements
	FETrussDomain& td = mesh.TrussDomain();
	int NT = td.size();
	for (i=0; i<NT; ++i)
	{
		FETrussElement& el = td.Element(i);
		td.UnpackElement(el);
		td.InternalForces(el, fe);
		AssembleResidual(el.m_node, el.LM(), fe, R);
	}

	// calculate inertial forces for dynamic problems
	if (m_fem.m_pStep->m_nanalysis == FE_DYNAMIC) InertialForces(R);

	// calculate forces due to pressure and add them to the residual
	// loop over surface elements
	int npr = m_fem.m_PC.size();
	for (i=0; i<npr; ++i)
	{
		FEPressureLoad& pc = m_fem.m_PC[i];
		if (pc.bc == 0)
		{
			FESurfaceElement& el = m_fem.m_psurf->Element(i);
			mesh.UnpackElement(el);

			// calculate nodal pressures
			double* pt = el.pt();

			double g = m_fem.GetLoadCurve(pc.lc)->Value();

			for (j=0; j<el.Nodes(); ++j) pt[j] = -g*pc.s[j];

			ndof = 3*el.Nodes();
			fe.resize(ndof);

			if (pc.blinear) LinearPressureForce(el, fe); else PressureForce(el, fe);

			// add element force vector to global force vector
			AssembleResidual(el.m_node, el.LM(), fe, R);
		}
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

void FESolidSolver::LinearConstraintForces(vector<double> &R)
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

void FESolidSolver::AssembleResidual(vector<int>& en, vector<int>& elm, vector<double>& fe, vector<double>& R)
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
//! calculates the equivalent nodal forces due to hydrostatic pressure

bool FESolidSolver::LinearPressureForce(FESurfaceElement& el, vector<double>& fe)
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

bool FESolidSolver::PressureForce(FESurfaceElement& el, vector<double>& fe)
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

bool FESolidSolver::PressureStiffness(FESurfaceElement& el, matrix& ke)
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
//! calculates the concentrated nodal forces

void FESolidSolver::NodalForces(vector<double>& F)
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
			
			// TODO: this next line only needs to be done when we
			//		 are using the symmetric version of poroelasticity!
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
//! This function calculates the inertial forces for dynamic problems

void FESolidSolver::InertialForces(vector<double>& R)
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
	FESolidDomain& bd = mesh.SolidDomain();
	for (iel=0; iel<bd.size(); ++iel)
	{
		FESolidElement& el = bd.Element(iel);
		bd.UnpackElement(el);

		FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(m_fem.GetMaterial(el.GetMatID()));

		double d = pme->Density();

		nint = el.GaussPoints();
		neln = el.Nodes();

		ke.Create(3*neln, 3*neln);
		ke.zero();

		fe.resize(3*neln);
		
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

void FESolidSolver::DiscreteElementForces(vector<double>& R)
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

void FESolidSolver::DiscreteElementStiffness()
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
