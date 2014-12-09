// FEStiffnessMatrix.cpp: implementation of the FEStiffnessMatrix class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEStiffnessMatrix.h"
#include "FEUT4Domain.h"
#include "FEPointConstraint.h"
#include "FEAugLagLinearConstraint.h"
#include "FECore/FERigidBody.h"
#include "FECore/DOFS.h"
#include "FERigidJoint.h"
#include "FEBioMech/FERigidSphericalJoint.h"
#include "FEBioMech/FERigidPinJoint.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEStiffnessMatrix::FEStiffnessMatrix(SparseMatrix* pK) : FEGlobalMatrix(pK)
{
}

FEStiffnessMatrix::~FEStiffnessMatrix()
{
}

//-----------------------------------------------------------------------------
//! Constructs the stiffness matrix from a FEM object. 
//! First a MatrixProfile object is constructed. This is done in two steps. First
//! the "static" profile is constructed which contains the contribution of the
//! static elements to the stiffness profile. The static profile is constructed 
//! only once during the first call to Create(). For next calls it is simply copied.
//! After the static profile is created (or copied) the dynamic elements are added
//! to the profile. Dynamic elements can change connectivity in between calls to
//! Create() and therefore have to be added explicitly every time.

bool FEStiffnessMatrix::Create(FEModel* pfem, int neq, bool breset)
{
	int i, j, k, l, m, n;

    // get nodal DOFS
    DOFS& fedofs = *DOFS::GetInstance();
    int MAX_NDOFS = fedofs.GetNDOFS();

	// keep a pointer to the FEM object
	m_pfem = pfem;
	FEModel& fem = *m_pfem;
	FEAnalysis* pstep = fem.GetCurrentStep();
	FEMesh& mesh = fem.GetMesh();

	// The first time we come here we build the "static" profile.
	// This static profile stores the contribution to the matrix profile
	// of the "elements" that do not change. Most elements are static except
	// for instance contact elements which can change connectivity in between
	// calls to the Create() function. Storing the static profile instead of
	// reconstructing it every time we come here saves us a lot of time. The 
	// static profile is stored in the variable m_MP.

	// begin building the profile
	build_begin(neq);
	{
		// The first time we are here we construct the "static"
		// profile. This profile contains the contribution from
		// all static elements. A static element is defined as
		// an element that never changes its connectity. This 
		// static profile is stored in the MP object. Next time
		// we come here we simply copy the MP object in stead
		// of building it from scratch.
		if (breset)
		{
			m_MPs.clear();

			vector<int> elm;

			// Add all elements to the profile
			// Loop over all active domains
			for (int nd=0; nd<pstep->Domains(); ++nd)
			{
				FEDomain& d = *pstep->Domain(nd);

				if (dynamic_cast<FEUT4Domain*>(&d) == 0)
				{
					for (int j=0; j<d.Elements(); ++j)
					{
						FEElement& el = d.ElementRef(j);
						if (!el.IsRigid())
						{
							d.UnpackLM(el, elm);
							build_add(elm);
						}
					}
				}
				else
				{
					// The UT4 Domain requires a slightly different form
					FEUT4Domain& ut4 = dynamic_cast<FEUT4Domain&>(d);

					// we'll need the node-element list
					FENodeElemList& NEL = ut4.GetNodeElemList();
					assert(NEL.Size() > 0);

					vector<int> LM;
					for (int i=0; i<mesh.Nodes(); ++i)
					{
						int NE = NEL.Valence(i);
						if (NE > 0)
						{
							LM.assign(NE*4*MAX_NDOFS, -1);
							FEElement** ppe = NEL.ElementList(i);
							for (int n=0; n<NE; ++n)
							{
								ut4.UnpackLM(*ppe[n], elm);
								for (int j=0; j<elm.size(); ++j) LM[n*4*MAX_NDOFS + j] = elm[j];
							}
							build_add(LM);
						}
					}
				}
			}

			// Add rigid bodies to the profile
			if (fem.Objects())
			{
				vector<int> lm(6);
				int nrb = fem.Objects();
				for (int i=0; i<nrb; ++i)
				{
					FERigidBody& rb = static_cast<FERigidBody&>(*fem.Object(i));
					for (int j=0; j<6; ++j) lm[j] = rb.m_LM[j];
					build_add(lm);
				}
			}

			// Add linear constraints to the profile
			// TODO: we need to add a function build_add(lmi, lmj) for
			// this type of "elements". Now we allocate too much memory
			if (fem.m_LinC.size() > 0)
			{
				int nlin = fem.m_LinC.size();
				vector<int> lm, elm;
				
				// do the cross-term
				// TODO: I have to make this easier. For instance,
				// keep a list that stores for each node the list of
				// elements connected to that node.
				// loop over all solid elements
				for (int nd=0; nd<pstep->Domains(); ++nd)
				{
					FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(pstep->Domain(nd));
					if (pbd)
					{
						for (i=0; i<pbd->Elements(); ++i)
						{
							FESolidElement& el = pbd->Element(i);
							if (!el.IsRigid())
							{
								pbd->UnpackLM(el, elm);
								int ne = elm.size();

								// see if this element connects to the 
								// master node of a linear constraint ...
								m = el.Nodes();
								for (j=0; j<m; ++j)
								{
									for (k=0; k<MAX_NDOFS; ++k)
									{
										n = fem.m_LCT[el.m_node[j]*MAX_NDOFS + k];

										if (n >= 0)
										{
											// ... it does so we need to connect the 
											// element to the linear constraint
											FELinearConstraint* plc = fem.m_LCA[n];

											int ns = plc->slave.size();

											lm.resize(ne + ns);
											for (l=0; l<ne; ++l) lm[l] = elm[l];

											list<FELinearConstraint::SlaveDOF>::iterator is = plc->slave.begin();
											for (l=ne; l<ne+ns; ++l, ++is) lm[l] = is->neq;

											build_add(lm);
		
											break;
										}
									}
								}
							}
						}
					}
				}

				// TODO: do the same thing for shell elements

				// do the constraint term
				int ni;
				list<FELinearConstraint>::iterator ic = fem.m_LinC.begin();
				n = 0;
				for (i=0; i<nlin; ++i, ++ic) n += ic->slave.size();
				lm.resize(n);
				ic = fem.m_LinC.begin();
				n = 0;
				for (i=0; i<nlin; ++i, ++ic)
				{
					ni = ic->slave.size();
					list<FELinearConstraint::SlaveDOF>::iterator is = ic->slave.begin();
					for (j=0; j<ni; ++j, ++is) lm[n++] = is->neq;
				}
				build_add(lm);
			}

			// do the nonlinear constraints
			int M = fem.NonlinearConstraints();
			for (int m=0; m<M; ++m)
			{
				FENLConstraint* pnlc = fem.NonlinearConstraint(m);
				if (dynamic_cast<FEPointConstraint*>(pnlc))
				{
					FEPointConstraint& pc = dynamic_cast<FEPointConstraint&>(*pnlc);
					vector<int> lm(3*9);
					FENode& n0 = mesh.Node(pc.m_node);
					lm[0] = n0.m_ID[DOF_X];
					lm[1] = n0.m_ID[DOF_Y];
					lm[2] = n0.m_ID[DOF_Z];
					for (j=0; i<8; ++i)
					{
						FENode& nj = mesh.Node(pc.m_pel->m_node[j]);
						lm[3*(j+1)  ] = nj.m_ID[DOF_X];
						lm[3*(j+1)+1] = nj.m_ID[DOF_Y];
						lm[3*(j+1)+2] = nj.m_ID[DOF_Z];
					}
					build_add(lm);
				}
				else if (dynamic_cast<FELinearConstraintSet*>(pnlc))
				{
					FELinearConstraintSet& lcs = dynamic_cast<FELinearConstraintSet&>(*pnlc);
					list<FEAugLagLinearConstraint*>& LC = lcs.m_LC;
					vector<int> lm;
					int N = LC.size();
					list<FEAugLagLinearConstraint*>::iterator it = LC.begin();
					for (i=0; i<N; ++i, ++it)
					{
						int n = (*it)->m_dof.size();
						lm.resize(n);
						FEAugLagLinearConstraint::Iterator is = (*it)->m_dof.begin();
						for (j=0; j<n; ++j, ++is) lm[j] = is->neq;
		
						build_add(lm);
					}
				}
				else if (dynamic_cast<FERigidJoint*>(pnlc))
				{
					FERigidJoint& rj = dynamic_cast<FERigidJoint&>(*pnlc);
					vector<int> lm(12);
			
					int* lm1 = dynamic_cast<FERigidBody*>(fem.Object(rj.m_nRBa))->m_LM;
					int* lm2 = dynamic_cast<FERigidBody*>(fem.Object(rj.m_nRBb))->m_LM;

					for (j=0; j<6; ++j) lm[j  ] = lm1[j];
					for (j=0; j<6; ++j) lm[j+6] = lm2[j];
					build_add(lm);
				}
				else if (dynamic_cast<FERigidSphericalJoint*>(pnlc))
				{
					FERigidSphericalJoint& rj = dynamic_cast<FERigidSphericalJoint&>(*pnlc);
					vector<int> lm(12);
                    
					int* lm1 = dynamic_cast<FERigidBody*>(fem.Object(rj.m_nRBa))->m_LM;
					int* lm2 = dynamic_cast<FERigidBody*>(fem.Object(rj.m_nRBb))->m_LM;
                    
					for (j=0; j<6; ++j) lm[j  ] = lm1[j];
					for (j=0; j<6; ++j) lm[j+6] = lm2[j];
					build_add(lm);
				}
				else if (dynamic_cast<FERigidPinJoint*>(pnlc))
				{
					FERigidPinJoint& rj = dynamic_cast<FERigidPinJoint&>(*pnlc);
					vector<int> lm(12);
                    
					int* lm1 = dynamic_cast<FERigidBody*>(fem.Object(rj.m_nRBa))->m_LM;
					int* lm2 = dynamic_cast<FERigidBody*>(fem.Object(rj.m_nRBb))->m_LM;
                    
					for (j=0; j<6; ++j) lm[j  ] = lm1[j];
					for (j=0; j<6; ++j) lm[j+6] = lm2[j];
					build_add(lm);
				}
			}

			// copy the static profile to the MP object
			// Make sure the LM buffer is flushed first.
			build_flush();
			m_MPs = *m_pMP;
		}
		else
		{
			// copy the old static profile
			*m_pMP = m_MPs;
		}

		// All following "elements" are nonstatic. That is, they can change
		// connectivity between calls to this function. All of these elements
		// are related to contact analysis (at this point).
		if (fem.SurfacePairInteractions() > 0)
		{
			// Add all contact interface elements
			for (i=0; i<fem.SurfacePairInteractions(); ++i)
			{
				FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairInteraction(i));
				if (pci->IsActive()) pci->BuildMatrixProfile(*this);
			}
		}
	}
	// All done! We can now finish building the profile and create 
	// the actual sparse matrix. This is done in the following function
	build_end();

	return true;
}
