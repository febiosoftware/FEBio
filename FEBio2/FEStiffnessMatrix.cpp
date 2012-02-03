// FEStiffnessMatrix.cpp: implementation of the FEStiffnessMatrix class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEStiffnessMatrix.h"
#include "fem.h"
#include "FEBioLib/FESlidingInterface.h"
#include "FEBioLib/FETiedInterface.h"
#include "FEBioLib/FERigidWallInterface.h"
#include "FEBioLib/FEFacet2FacetSliding.h"
#include "FEBioLib/FESlidingInterface2.h"
#include "FEBioLib/FESlidingInterface3.h"
#include "FEBioLib/FEPeriodicBoundary.h"
#include "FEBioLib/FESurfaceConstraint.h"
#include "FESolver.h"
#include "FEBioLib/FEUT4Domain.h"
#include "FEBioLib/FEPointConstraint.h"
#include "FEBioLib/FEAugLagLinearConstraint.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEStiffnessMatrix::FEStiffnessMatrix(SparseMatrix* pK)
{
	m_pA = pK;
	m_LM.resize(MAX_LM_SIZE);
	m_pMP = 0;
	m_nlm = 0;
}

FEStiffnessMatrix::~FEStiffnessMatrix()
{
	delete m_pA;
	m_pA = 0;
	if (m_pMP) delete m_pMP;
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

bool FEStiffnessMatrix::Create(FESolver* psolver, int neq, bool breset)
{
	int i, j, k, l, m, n;

	// keep a pointer to the FEM object
	FESolver& solver = *psolver;
	FEM& fem = dynamic_cast<FEM&>(psolver->GetFEModel());
	m_pfem = &fem;

	// The first time we come here we build the "static" profile.
	// This static profile stores the contribution to the matrix profile
	// of the "elements" that do not change. Most elements are static except
	// for instance contact elements which can change connectivity in between
	// calls to the Create() function. Storing the static profile instead of
	// reconstructing it every time we come here saves us a lot of time.
	static SparseMatrixProfile MP;

	// begin building the profile
	build_begin(neq);
	{
		FEMesh& mesh = fem.m_mesh;

		// The first time we are here we construct the "static"
		// profile. This profile contains the contribution from
		// all static elements. A static element is defined as
		// an element that never changes its connectity. This 
		// static profile is stored in the MP object. Next time
		// we come here we simply copy the MP object in stead
		// of building it from scratch.
		if (breset)
		{
			MP.clear();

			vector<int> elm;

			// Add all elements to the profile
			// Loop over all active domains
			for (int nd=0; nd<solver.Domains(); ++nd)
			{
				FEDomain& d = *solver.Domain(nd);

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
								FESolidElement& el = dynamic_cast<FESolidElement&>(*ppe[n]);
								ut4.UnpackLM(el, elm);
								for (int j=0; j<4*MAX_NDOFS; ++j) LM[n*4*MAX_NDOFS + j] = elm[j];
							}
							build_add(LM);
						}
					}
				}
			}

			// Add rigid bodies to the profile
			if (fem.m_RB.empty() == false)
			{
				vector<int> lm(6);
				int nrb = fem.m_RB.size();
				for (int i=0; i<nrb; ++i)
				{
					FERigidBody& rb = fem.m_RB[i];
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
				for (int nd=0; nd<solver.Domains(); ++nd)
				{
					FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(solver.Domain(nd));
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
				switch (pnlc->Type())
				{
				case FE_POINT_CONSTRAINT: // add point constraints
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
					break;
				case FE_LINEAR_CONSTRAINT:
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
					break;
				case FE_RIGID_JOINT:
					{
						FERigidJoint& rj = dynamic_cast<FERigidJoint&>(*pnlc);
						vector<int> lm(12);
			
						int* lm1 = fem.m_RB[ rj.m_nRBa ].m_LM;
						int* lm2 = fem.m_RB[ rj.m_nRBb ].m_LM;

						for (j=0; j<6; ++j) lm[j  ] = lm1[j];
						for (j=0; j<6; ++j) lm[j+6] = lm2[j];
						build_add(lm);
					}
					break;
				default:
					assert(false);
				}
			}

			// copy the static profile to the MP object
			// Make sure the LM buffer is flushed first.
			build_flush();
			MP = *m_pMP;
		}
		else
		{
			// copy the old static profile
			*m_pMP = MP;
		}

		// All following "elements" are nonstatic. That is, they can change
		// connectivity between calls to this function. All of these elements
		// are related to contact analysis (at this point).
		if (fem.ContactInterfaces() > 0)
		{
			int *id;
			// Add all contact interface elements
			for (i=0; i<fem.ContactInterfaces(); ++i)
			{
				FEContactInterface* pci = fem.ContactInterface(i);

				// add sliding interface elements
				FESlidingInterface* psi = dynamic_cast<FESlidingInterface*>(pci);
				if (psi)
				{
					vector<int> lm(6*5);

					int npass = (psi->m_btwo_pass?2:1);
					for (int np=0; np<npass; ++np)
					{
						FESlidingSurface& ss = (np==0? psi->m_ss : psi->m_ms);
						FESlidingSurface& ms = (np==0? psi->m_ms : psi->m_ss);

						for (j=0; j<ss.Nodes(); ++j)
						{
							FESurfaceElement* pe = ss.m_pme[j];

							if (pe != 0)
							{
								FESurfaceElement& me = *pe;
								int* en = &me.m_node[0];

								// Note that we need to grab the rigid degrees of freedom as well
								// this is in case one of the nodes belongs to a rigid body.
								n = me.Nodes();
								if (n == 3)
								{
									lm[6*(3+1)  ] = -1;
									lm[6*(3+1)+1] = -1;
									lm[6*(3+1)+2] = -1;
									lm[6*(3+1)+3] = -1;
									lm[6*(3+1)+4] = -1;
									lm[6*(3+1)+5] = -1;
								}

								lm[0] = ss.Node(j).m_ID[DOF_X];
								lm[1] = ss.Node(j).m_ID[DOF_Y];
								lm[2] = ss.Node(j).m_ID[DOF_Z];
								lm[3] = ss.Node(j).m_ID[DOF_RU];
								lm[4] = ss.Node(j).m_ID[DOF_RV];
								lm[5] = ss.Node(j).m_ID[DOF_RW];

								for (k=0; k<n; ++k)
								{
									id = fem.m_mesh.Node(en[k]).m_ID;
									lm[6*(k+1)  ] = id[DOF_X];
									lm[6*(k+1)+1] = id[DOF_Y];
									lm[6*(k+1)+2] = id[DOF_Z];
									lm[6*(k+1)+3] = id[DOF_RU];
									lm[6*(k+1)+4] = id[DOF_RV];
									lm[6*(k+1)+5] = id[DOF_RW];
								}

								build_add(lm);
							}
						}
					}
				}

				// facet-to-facet sliding interfaces
				FEFacet2FacetSliding* pfi = dynamic_cast<FEFacet2FacetSliding*>(pci);
				if (pfi)
				{
					vector<int> lm(6*8);

					int npass = (pfi->m_btwo_pass?2:1);
					for (int np=0; np<npass; ++np)
					{
						FEFacetSlidingSurface& ss = (np == 0? pfi->m_ss : pfi->m_ms);
						FEFacetSlidingSurface& ms = (np == 0? pfi->m_ms : pfi->m_ss);

						int ni = 0, k, l;
						for (j=0; j<ss.Elements(); ++j)
						{
							FESurfaceElement& se = ss.Element(j);
							int nint = se.GaussPoints();
							int* sn = &se.m_node[0];
							for (k=0; k<nint; ++k, ++ni)
							{
								FESurfaceElement* pe = ss.m_pme[ni];
								if (pe != 0)
								{
									FESurfaceElement& me = dynamic_cast<FESurfaceElement&> (*pe);
									int* mn = &me.m_node[0];

									set(lm, -1);

									int nseln = se.Nodes();
									int nmeln = me.Nodes();

									for (l=0; l<nseln; ++l)
									{
										id = fem.m_mesh.Node(sn[l]).m_ID;
										lm[6*l  ] = id[DOF_X];
										lm[6*l+1] = id[DOF_Y];
										lm[6*l+2] = id[DOF_Z];
										lm[6*l+3] = id[DOF_RU];
										lm[6*l+4] = id[DOF_RV];
										lm[6*l+5] = id[DOF_RW];
									}

									for (l=0; l<nmeln; ++l)
									{
										id = fem.m_mesh.Node(mn[l]).m_ID;
										lm[6*(l+nseln)  ] = id[DOF_X];
										lm[6*(l+nseln)+1] = id[DOF_Y];
										lm[6*(l+nseln)+2] = id[DOF_Z];
										lm[6*(l+nseln)+3] = id[DOF_RU];
										lm[6*(l+nseln)+4] = id[DOF_RV];
										lm[6*(l+nseln)+5] = id[DOF_RW];
									}

									build_add(lm);
								}
							}
						}
					}
				}

				// sliding2 interfaces
				FESlidingInterface2* ps2 = dynamic_cast<FESlidingInterface2*>(pci);
				if (ps2)
				{
					vector<int> lm(7*8);

					int npass = (ps2->m_btwo_pass?2:1);
					for (int np=0; np<npass; ++np)
					{
						FESlidingSurface2& ss = (np == 0? ps2->m_ss : ps2->m_ms);
						FESlidingSurface2& ms = (np == 0? ps2->m_ms : ps2->m_ss);

						int ni = 0, k, l;
						for (j=0; j<ss.Elements(); ++j)
						{
							FESurfaceElement& se = ss.Element(j);
							int nint = se.GaussPoints();
							int* sn = &se.m_node[0];
							for (k=0; k<nint; ++k, ++ni)
							{
								FESurfaceElement* pe = ss.m_pme[ni];
								if (pe != 0)
								{
									FESurfaceElement& me = dynamic_cast<FESurfaceElement&> (*pe);
									int* mn = &me.m_node[0];

									set(lm, -1);

									int nseln = se.Nodes();
									int nmeln = me.Nodes();

									for (l=0; l<nseln; ++l)
									{
										id = fem.m_mesh.Node(sn[l]).m_ID;
										lm[7*l  ] = id[DOF_X];
										lm[7*l+1] = id[DOF_Y];
										lm[7*l+2] = id[DOF_Z];
										lm[7*l+3] = id[DOF_P];
										lm[7*l+4] = id[DOF_RU];
										lm[7*l+5] = id[DOF_RV];
										lm[7*l+6] = id[DOF_RW];
									}

									for (l=0; l<nmeln; ++l)
									{
										id = fem.m_mesh.Node(mn[l]).m_ID;
										lm[7*(l+nseln)  ] = id[DOF_X];
										lm[7*(l+nseln)+1] = id[DOF_Y];
										lm[7*(l+nseln)+2] = id[DOF_Z];
										lm[7*(l+nseln)+3] = id[DOF_P];
										lm[7*(l+nseln)+4] = id[DOF_RU];
										lm[7*(l+nseln)+5] = id[DOF_RV];
										lm[7*(l+nseln)+6] = id[DOF_RW];
									}

									build_add(lm);
								}
							}
						}
					}
				}

				// sliding3 interfaces
				FESlidingInterface3* ps3 = dynamic_cast<FESlidingInterface3*>(pci);
				if (ps3)
				{
					vector<int> lm(8*8);
					
					int npass = (ps3->m_btwo_pass?2:1);
					for (int np=0; np<npass; ++np)
					{
						FESlidingSurface3& ss = (np == 0? ps3->m_ss : ps3->m_ms);
						FESlidingSurface3& ms = (np == 0? ps3->m_ms : ps3->m_ss);
						
						int ni = 0, k, l;
						for (j=0; j<ss.Elements(); ++j)
						{
							FESurfaceElement& se = ss.Element(j);
							int nint = se.GaussPoints();
							int* sn = &se.m_node[0];
							for (k=0; k<nint; ++k, ++ni)
							{
								FESurfaceElement* pe = ss.m_pme[ni];
								if (pe != 0)
								{
									FESurfaceElement& me = dynamic_cast<FESurfaceElement&> (*pe);
									int* mn = &me.m_node[0];
									
									set(lm, -1);
									
									int nseln = se.Nodes();
									int nmeln = me.Nodes();
									
									for (l=0; l<nseln; ++l)
									{
										id = fem.m_mesh.Node(sn[l]).m_ID;
										lm[8*l  ] = id[DOF_X];
										lm[8*l+1] = id[DOF_Y];
										lm[8*l+2] = id[DOF_Z];
										lm[8*l+3] = id[DOF_P];
										lm[8*l+4] = id[DOF_RU];
										lm[8*l+5] = id[DOF_RV];
										lm[8*l+6] = id[DOF_RW];
										lm[8*l+7] = id[DOF_C];
									}
									
									for (l=0; l<nmeln; ++l)
									{
										id = fem.m_mesh.Node(mn[l]).m_ID;
										lm[8*(l+nseln)  ] = id[DOF_X];
										lm[8*(l+nseln)+1] = id[DOF_Y];
										lm[8*(l+nseln)+2] = id[DOF_Z];
										lm[8*(l+nseln)+3] = id[DOF_P];
										lm[8*(l+nseln)+4] = id[DOF_RU];
										lm[8*(l+nseln)+5] = id[DOF_RV];
										lm[8*(l+nseln)+6] = id[DOF_RW];
										lm[8*(l+nseln)+7] = id[DOF_C];
									}
									
									build_add(lm);
								}
							}
						}
					}
				}
				
				// add tied interface elements
				FETiedInterface* pti = dynamic_cast<FETiedInterface*>(pci);
				if (pti)
				{
					vector<int> lm(6*5);

					FETiedContactSurface& ss = pti->ss;
					FETiedContactSurface& ms = pti->ms;

					for (j=0; j<ss.Nodes(); ++j)
					{
						FEElement* pe = ss.m_pme[j];

						if (pe != 0)
						{
							FESurfaceElement& me = dynamic_cast<FESurfaceElement&> (*pe);
							int* en = &me.m_node[0];

							n = me.Nodes();
							if (n == 3)
							{
								lm[6*(3+1)  ] = -1;
								lm[6*(3+1)+1] = -1;
								lm[6*(3+1)+2] = -1;
								lm[6*(3+1)+3] = -1;
								lm[6*(3+1)+4] = -1;
								lm[6*(3+1)+5] = -1;
							}

							lm[0] = ss.Node(j).m_ID[DOF_X];
							lm[1] = ss.Node(j).m_ID[DOF_Y];
							lm[2] = ss.Node(j).m_ID[DOF_Z];
							lm[3] = ss.Node(j).m_ID[DOF_RU];
							lm[4] = ss.Node(j).m_ID[DOF_RV];
							lm[5] = ss.Node(j).m_ID[DOF_RW];

							for (k=0; k<n; ++k)
							{
								id = fem.m_mesh.Node(en[k]).m_ID;
								lm[6*(k+1)  ] = id[DOF_X];
								lm[6*(k+1)+1] = id[DOF_Y];
								lm[6*(k+1)+2] = id[DOF_Z];
								lm[6*(k+1)+3] = id[DOF_RU];
								lm[6*(k+1)+4] = id[DOF_RV];
								lm[6*(k+1)+5] = id[DOF_RW];
							}

							build_add(lm);
						}
					}
				}

				// add periodic boundary elements
				// TODO: what if two_pass ??
				FEPeriodicBoundary* pbi = dynamic_cast<FEPeriodicBoundary*>(pci);
				if (pbi)
				{
					vector<int> lm(6*5);

					FEPeriodicSurface& ss = pbi->m_ss;
					FEPeriodicSurface& ms = pbi->m_ms;

					for (j=0; j<ss.Nodes(); ++j)
					{
						FESurfaceElement& me = *ss.m_pme[j];
						int* en = &me.m_node[0];

						n = me.Nodes();
						if (n == 3)
						{
							lm[6*(3+1)  ] = -1;
							lm[6*(3+1)+1] = -1;
							lm[6*(3+1)+2] = -1;
							lm[6*(3+1)+3] = -1;
							lm[6*(3+1)+4] = -1;
							lm[6*(3+1)+5] = -1;
						}

						lm[0] = ss.Node(j).m_ID[DOF_X];
						lm[1] = ss.Node(j).m_ID[DOF_Y];
						lm[2] = ss.Node(j).m_ID[DOF_Z];
						lm[3] = ss.Node(j).m_ID[DOF_RU];
						lm[4] = ss.Node(j).m_ID[DOF_RV];
						lm[5] = ss.Node(j).m_ID[DOF_RW];

						for (k=0; k<n; ++k)
						{
							id = fem.m_mesh.Node(en[k]).m_ID;
							lm[6*(k+1)  ] = id[DOF_X];
							lm[6*(k+1)+1] = id[DOF_Y];
							lm[6*(k+1)+2] = id[DOF_Z];
							lm[6*(k+1)+3] = id[DOF_RU];
							lm[6*(k+1)+4] = id[DOF_RV];
							lm[6*(k+1)+5] = id[DOF_RW];
						}

						build_add(lm);
					}
				}

				// add surface constraints
				// TODO: what if two_pass ??
				FESurfaceConstraint* psc = dynamic_cast<FESurfaceConstraint*>(pci);
				if (psc)
				{
					vector<int> lm(6*5);

					FESurfaceConstraintSurface& ss = psc->m_ss;
					FESurfaceConstraintSurface& ms = psc->m_ms;

					int nref = ss.m_nref;
					FESurfaceElement* pref = ss.m_pme[nref];

					int n0 = pref->Nodes();
					int nr0[4];
					for (j=0; j<n0; ++j) nr0[j] = pref->m_node[j];

					set(lm, -1);

					lm[0] = ss.Node(nref).m_ID[DOF_X];
					lm[1] = ss.Node(nref).m_ID[DOF_Y];
					lm[2] = ss.Node(nref).m_ID[DOF_Z];
					lm[3] = ss.Node(nref).m_ID[DOF_RU];
					lm[4] = ss.Node(nref).m_ID[DOF_RV];
					lm[5] = ss.Node(nref).m_ID[DOF_RW];

					for (k=0; k<n0; ++k)
					{
						id = fem.m_mesh.Node(nr0[k]).m_ID;
						lm[6*(k+1)  ] = id[DOF_X];
						lm[6*(k+1)+1] = id[DOF_Y];
						lm[6*(k+1)+2] = id[DOF_Z];
						lm[6*(k+1)+3] = id[DOF_RU];
						lm[6*(k+1)+4] = id[DOF_RV];
						lm[6*(k+1)+5] = id[DOF_RW];
					}

					for (j=0; j<ss.Nodes(); ++j)
					{
						FESurfaceElement& me = *ss.m_pme[j];
						int* en = &me.m_node[0];

						set(lm, -1);

						n = me.Nodes();

						lm[0] = ss.Node(j).m_ID[DOF_X];
						lm[1] = ss.Node(j).m_ID[DOF_Y];
						lm[2] = ss.Node(j).m_ID[DOF_Z];
						lm[3] = ss.Node(j).m_ID[DOF_RU];
						lm[4] = ss.Node(j).m_ID[DOF_RV];
						lm[5] = ss.Node(j).m_ID[DOF_RW];

						for (k=0; k<n; ++k)
						{
							id = fem.m_mesh.Node(en[k]).m_ID;
							lm[6*(k+1)  ] = id[DOF_X];
							lm[6*(k+1)+1] = id[DOF_Y];
							lm[6*(k+1)+2] = id[DOF_Z];
							lm[6*(k+1)+3] = id[DOF_RU];
							lm[6*(k+1)+4] = id[DOF_RV];
							lm[6*(k+1)+5] = id[DOF_RW];
						}

						build_add(lm);
					}
				}

				// add rigid wall elements
				FERigidWallInterface* pri = dynamic_cast<FERigidWallInterface*>(pci);
				if (pri)
				{
					vector<int> lm(6);
					FERigidWallSurface& ss = pri->m_ss;

					for (j=0; j<ss.Nodes(); ++j)
					{
						if (ss.gap[j] >= 0)
						{
							lm[0] = ss.Node(j).m_ID[DOF_X];
							lm[1] = ss.Node(j).m_ID[DOF_Y];
							lm[2] = ss.Node(j).m_ID[DOF_Z];
							lm[3] = ss.Node(j).m_ID[DOF_RU];
							lm[4] = ss.Node(j).m_ID[DOF_RV];
							lm[5] = ss.Node(j).m_ID[DOF_RW];

							build_add(lm);
						}
					}
				}
			}
		}
	}
	// All done! We can now finish building the profile and create 
	// the actual sparse matrix. This is done in the following function
	build_end();

	return true;
}

//-----------------------------------------------------------------------------
//! Start building the profile. That is delete the old profile (if there was one)
//! and create a new one. 
void FEStiffnessMatrix::build_begin(int neq)
{
	if (m_pMP) delete m_pMP;
	m_pMP = new SparseMatrixProfile(neq);
	m_nlm = 0;
}

//-----------------------------------------------------------------------------
//! Add an "element" to the matrix profile. The definition of an element is quite
//! general at this point. Any two or more degrees of freedom that are connected 
//! somehow are considered an element. This function takes one argument, namely a
//! list of degrees of freedom that are part of such an "element". Elements are 
//! first copied into a local LM array.  When this array is full it is flushed and
//! all elements in this array are added to the matrix profile.
void FEStiffnessMatrix::build_add(vector<int>& lm)
{
	m_LM[m_nlm++] = lm;
	if (m_nlm >= MAX_LM_SIZE) build_flush();
}

//-----------------------------------------------------------------------------
//! Flush the LM array. The LM array stores a buffer of elements that have to be
//! added to the profile. When this buffer is full it needs to be flushed. This
//! flushin operation causes the actual update of the matrix profile.
void FEStiffnessMatrix::build_flush()
{
	int i, j, n, *lm;

	// Since prescribed dofs have an equation number of < -1 we need to modify that
	// otherwise no storage will be allocated for these dofs (even not diagonal elements!).
	for (i=0; i<m_nlm; ++i)
	{
		n = m_LM[i].size();
		lm = &(m_LM[i])[0];
		for (j=0; j<n; ++j) if (lm[j] < -1) lm[j] = -lm[j]-2;
	}

	m_pMP->UpdateProfile(m_LM, m_nlm);
	m_nlm = 0;
}

//-----------------------------------------------------------------------------
//! This function makes sure the LM buffer is flushed and creates the actual
//! sparse matrix from the matrix profile.
void FEStiffnessMatrix::build_end()
{
	if (m_nlm > 0) build_flush();
	m_pA->Create(*m_pMP);
}
