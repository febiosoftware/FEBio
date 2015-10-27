#include "stdafx.h"
#include "FEMortarSlidingContact.h"
#include "FEStiffnessMatrix.h"
#include <FECore\FEModel.h>

//=============================================================================
// FEMortarSurface
//=============================================================================

//-----------------------------------------------------------------------------
FEMortarSurface::FEMortarSurface(FEMesh* pm) : FEContactSurface(pm) {}

//-----------------------------------------------------------------------------
bool FEMortarSurface::Init()
{
	// always intialize base class first!
	if (FEContactSurface::Init() == false) return false;

	// get the number of nodes
	int NN = Nodes();

	// allocate data structures
	m_p.resize(NN, 0.0);
	m_nu.resize(NN, vec3d(0,0,0));
	m_gap.resize(NN, vec3d(0,0,0));

	return true;
}

//-----------------------------------------------------------------------------
void FEMortarSurface::UpdateNormals()
{
	const int NN = Nodes();

	// reset normals
	m_nu.assign(NN, vec3d(0,0,0));

	// recalculate normals
	const int NF = Elements();
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& el = Element(i);
		int neln = el.Nodes();
		for (int j=0; j<neln; ++j)
		{
			vec3d r0 = Node(el.m_lnode[j]).m_rt;
			vec3d r1 = Node(el.m_lnode[(j+1)%neln]).m_rt;
			vec3d r2 = Node(el.m_lnode[(j+neln-1)%neln]).m_rt;

			vec3d n = (r1 - r0)^(r2 - r0);
			m_nu[el.m_lnode[j]] += n;
		}
	}

	// normalize
	for (int i=0; i<NN; ++i) m_nu[i].unit();
}

//=============================================================================
// FEMortarSlidingContact
//=============================================================================

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_PARAMETER_LIST(FEMortarSlidingContact, FEContactInterface)
	ADD_PARAMETER(m_blaugon      , FE_PARAM_BOOL  , "laugon"       ); 
	ADD_PARAMETER(m_atol         , FE_PARAM_DOUBLE, "tolerance"    );
	ADD_PARAMETER(m_eps          , FE_PARAM_DOUBLE, "penalty"      );
	ADD_PARAMETER(m_naugmin      , FE_PARAM_INT   , "minaug"       );
	ADD_PARAMETER(m_naugmax      , FE_PARAM_INT   , "maxaug"       );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEMortarSlidingContact::FEMortarSlidingContact(FEModel* pfem) : FEContactInterface(pfem), m_ss(&pfem->GetMesh()), m_ms(&pfem->GetMesh())
{
}

//-----------------------------------------------------------------------------
FEMortarSlidingContact::~FEMortarSlidingContact()
{
}

//-----------------------------------------------------------------------------
bool FEMortarSlidingContact::Init()
{
	// initialize surfaces
	if (m_ms.Init() == false) return false;
	if (m_ss.Init() == false) return false;

	// set the integration rule
	m_pT = dynamic_cast<FESurfaceElementTraits*>(FEElementLibrary::GetElementTraits(FE_TRI3G3));

	// allocate sturcture for the integration weights
	int NS = m_ss.Nodes();
	int NM = m_ms.Nodes();
	m_n1.resize(NS,NS);
	m_n2.resize(NS,NM);

	return true;
}

//-----------------------------------------------------------------------------
void FEMortarSlidingContact::Activate()
{
	//! don't forget the base class
	FEContactInterface::Activate();

	// update the normals on the slave surface
	m_ss.UpdateNormals();

	// update the mortar weights
	UpdateMortarWeights();

	// update the nodal gaps
	// (must be done after mortar eights are updated)
	UpdateNodalGaps();
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FEMortarSlidingContact::BuildMatrixProfile(FEStiffnessMatrix& K)
{
	// For now we'll assume that each node on the slave side is connected to the master side
	// This is obviously too much, but we'll worry about improving this later
	int NS = m_ss.Nodes();
	int NM = m_ms.Nodes();
	vector<int> LM(3*(NS+NM));
	for (int i=0; i<NS; ++i)
	{
		FENode& ni = m_ss.Node(i);
		LM[3*i  ] = ni.m_ID[0];
		LM[3*i+1] = ni.m_ID[1];
		LM[3*i+2] = ni.m_ID[2];
	}
	for (int i=0; i<NM; ++i)
	{
		FENode& ni = m_ms.Node(i);
		LM[3*NS + 3*i  ] = ni.m_ID[0];
		LM[3*NS + 3*i+1] = ni.m_ID[1];
		LM[3*NS + 3*i+2] = ni.m_ID[2];
	}
	K.build_add(LM);
}

//-----------------------------------------------------------------------------
//! Update the mortar weights
void FEMortarSlidingContact::UpdateMortarWeights()
{
	// clear weights
	m_n1.zero();
	m_n2.zero();

	// We'll use a patch to store the intersection between two facets
	Patch patch;

	// number of integration points
	const int MAX_INT = 11;
	const int nint = m_pT->nint;
	vector<double>& gw = m_pT->gw;
	vector<double>& gr = m_pT->gr;
	vector<double>& gs = m_pT->gs;

	// These arrays will store the shape function values of the projection points 
	// on the slave and master side when evaluating the integral over a pallet
	double Ns[MAX_INT][4], Nm[MAX_INT][4];

	// loop over all non-mortar facets
	int NSF = m_ss.Elements();
	int NMF = m_ms.Elements();
	for (int i=0; i<NSF; ++i)
	{
		// get the non-mortar surface element
		FESurfaceElement& se = m_ss.Element(i);

		vector<vec3d> x1(FEElement::MAX_NODES);
		for (int m=0; m<se.Nodes(); ++m) x1[m] = m_ss.Node(se.m_lnode[m]).m_rt;

		// loop over all the mortar surface elements
		for (int j=0; j<NMF; ++j)
		{
			// get the next surface element
			FESurfaceElement& me = m_ms.Element(j);

			vector<vec3d> x2(FEElement::MAX_NODES);
			for (int m=0; m<me.Nodes(); ++m) x2[m] = m_ms.Node(me.m_lnode[m]).m_rt;

			// calculate the patch of triangles, representing the intersection
			// of the non-mortar facet with the mortar facet
			if (CalculateIntersection(i, j, patch))
			{
				// loop over all patches
				int np = patch.Size();
				for (int k=0; k<np; ++k)
				{
					// get the next facet
					Patch::FACET& fk = patch.Facet(k);

					// calculate the patch area
					// (We multiply by two because the sum of the integration weights in FEBio sum up to the area
					// of the triangle in natural coordinates (=0.5)).
					double Area = fk.Area()*2.0;

					// loop over integration points
					for (int n=0; n<nint; ++n)
					{
						// evaluate the spatial position of the integration point on the patch
						vec3d xp = fk.Position(gr[n], gs[n]);

						// evaluate the integration points on the slave and master surfaces
						// i.e. determine rs, rm
						double r1 = 0, s1 = 0, r2 = 0, s2 = 0;
						m_ss.ProjectToSurface(se, xp, r1, s1);
						m_ms.ProjectToSurface(me, xp, r2, s2);

						vec3d xa = se.eval(&x1[0], r1, s1);
						vec3d xb = me.eval(&x2[0], r2, s2);

						// evaluate shape functions
						se.shape_fnc(Ns[n], r1, s1);
						me.shape_fnc(Nm[n], r2, s2);
					}

					// Evaluate the contributions to the integrals
					int nA = se.Nodes();
					int nC = me.Nodes();
					for (int A=0; A<nA; ++A)
					{
						int a = se.m_lnode[A];

						// loop over all the nodes on the slave facet
						for (int B=0; B<nA; ++B)
						{
							double n1 = 0;
							for (int n=0; n<nint; ++n)
							{
								n1 += gw[n]*Ns[n][A]*Ns[n][B];
							}
							n1 *= Area;

							int b = se.m_lnode[B];
							m_n1[a][b] += n1;
						}

						// loop over all the nodes on the master facet
						for (int C = 0; C<nC; ++C)
						{
							double n2 = 0;
							for (int n=0; n<nint; ++n)
							{
								n2 += gw[n]*Ns[n][A]*Nm[n][C];
							}
							n2 *= Area;

							int c = me.m_lnode[C];
							m_n2[a][c] += n2;
						}
					}
				}
			}		
		}
	}
}

//-----------------------------------------------------------------------------
//! calculate contact forces
void FEMortarSlidingContact::ContactForces(FEGlobalVector& R)
{
	int NS = m_ss.Nodes();
	int NM = m_ms.Nodes();

	// loop over all slave nodes
	for (int A=0; A<NS; ++A)
	{
		vec3d nuA = m_ss.m_nu[A];
		vec3d gA = m_ss.m_gap[A];
		double pA = gA*nuA;
		vec3d tA = nuA*(pA*m_eps);

		// loop over all slave nodes
		vector<int> en(1);
		vector<int> lm(3);
		vector<double> fe(3);
		for (int B=0; B<NS; ++B)
		{
			FENode& nodeB = m_ss.Node(B);
			en[0] = m_ss.m_node[B];
			lm[0] = nodeB.m_ID[0];
			lm[1] = nodeB.m_ID[1];
			lm[2] = nodeB.m_ID[2];

			double nAB = m_n1[A][B];

			fe[0] = -tA.x*nAB;
			fe[1] = -tA.y*nAB;
			fe[2] = -tA.z*nAB;

			R.Assemble(en, lm, fe);
		}

		// loop over master side
		for (int C=0; C<NM; ++C)
		{
			FENode& nodeC = m_ms.Node(C);
			en[0] = m_ms.m_node[C];
			lm[0] = nodeC.m_ID[0];
			lm[1] = nodeC.m_ID[1];
			lm[2] = nodeC.m_ID[2];

			double nAC = m_n2[A][C];

			fe[0] = tA.x*nAC;
			fe[1] = tA.y*nAC;
			fe[2] = tA.z*nAC;

			R.Assemble(en, lm, fe);
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculate the intersection between two facets
bool FEMortarSlidingContact::CalculateIntersection(int k, int l, Patch& patch)
{
	// clear the patch
	patch.Clear();

	if (k != l) return false;

	// calculate the intersection
	FESurfaceElement& ek = m_ss.Element(k);
	FESurfaceElement& el = m_ms.Element(l);

	vec3d r[3];
	r[0] = m_ss.Node(ek.m_lnode[0]).m_rt;
	r[1] = m_ss.Node(ek.m_lnode[1]).m_rt;
	r[2] = m_ss.Node(ek.m_lnode[2]).m_rt;
	patch.Add(r);

	r[0] = m_ss.Node(ek.m_lnode[2]).m_rt;
	r[1] = m_ss.Node(ek.m_lnode[3]).m_rt;
	r[2] = m_ss.Node(ek.m_lnode[0]).m_rt;
	patch.Add(r);

	// return
	return (patch.Empty() == false);
}

//-----------------------------------------------------------------------------
//! Update the nodal gaps
void FEMortarSlidingContact::UpdateNodalGaps()
{
	// reset nodal gaps
	vector<vec3d>& gap = m_ss.m_gap;
	zero(m_ss.m_gap);

	int NS = m_ss.Nodes();
	int NM = m_ms.Nodes();

	// loop over all slave nodes
	for (int A=0; A<NS; ++A)
	{
		// loop over all slave nodes
		for (int B=0; B<NS; ++B)
		{
			FENode& nodeB = m_ss.Node(B);
			vec3d& xB = nodeB.m_rt;
			double nAB = m_n1[A][B];
			gap[A] += xB*nAB;
		}

		// loop over master side
		for (int C=0; C<NM; ++C)
		{
			FENode& nodeC = m_ms.Node(C);
			vec3d& xC = nodeC.m_rt;
			double nAC = m_n2[A][C];
			gap[A] -= xC*nAC;
		}
	}
}

//-----------------------------------------------------------------------------
//! calculate contact stiffness
void FEMortarSlidingContact::ContactStiffness(FESolver* psolver)
{
	int NS = m_ss.Nodes();
	int NM = m_ms.Nodes();

	vector<int> lmi(3), lmj(3);
	matrix kA(3,3), kG(3,3), ke(3,3);
	for (int A=0; A<NS; ++A)
	{
		vec3d nuA = m_ss.m_nu[A];

		// loop over all slave nodes
		for (int B=0; B<NS; ++B)
		{
			FENode& nodeB = m_ss.Node(B);
			lmi[0] = nodeB.m_ID[0];
			lmi[1] = nodeB.m_ID[1];
			lmi[2] = nodeB.m_ID[2];

			double nAB = m_n1[A][B];

			kA[0][0] = m_eps*nAB*(nuA.x*nuA.x); kA[0][1] = m_eps*nAB*(nuA.x*nuA.y); kA[0][2] = m_eps*nAB*(nuA.x*nuA.z);
			kA[1][0] = m_eps*nAB*(nuA.y*nuA.x); kA[1][1] = m_eps*nAB*(nuA.y*nuA.y); kA[1][2] = m_eps*nAB*(nuA.y*nuA.z);
			kA[2][0] = m_eps*nAB*(nuA.z*nuA.x); kA[2][1] = m_eps*nAB*(nuA.z*nuA.y); kA[2][2] = m_eps*nAB*(nuA.z*nuA.z);

			// loop over slave nodes
			for (int C=0; C<NS; ++C)
			{
				FENode& nodeC = m_ss.Node(C);
				lmj[0] = nodeC.m_ID[0];
				lmj[1] = nodeC.m_ID[1];
				lmj[2] = nodeC.m_ID[2];

				double nAC = m_n1[A][C];

				kG[0][0] = nAC; kG[0][1] = 0.0; kG[0][2] = 0.0;
				kG[1][0] = 0.0; kG[1][1] = nAC; kG[1][2] = 0.0;
				kG[2][0] = 0.0; kG[2][1] = 0.0; kG[2][2] = nAC;

				ke = kA*kG;

				psolver->AssembleStiffness2(lmi, lmj, ke);
			}

			// loop over master nodes
			for (int C=0; C<NM; ++C)
			{
				FENode& nodeC = m_ms.Node(C);
				lmj[0] = nodeC.m_ID[0];
				lmj[1] = nodeC.m_ID[1];
				lmj[2] = nodeC.m_ID[2];

				double nAC = -m_n2[A][C];

				kG[0][0] = nAC; kG[0][1] = 0.0; kG[0][2] = 0.0;
				kG[1][0] = 0.0; kG[1][1] = nAC; kG[1][2] = 0.0;
				kG[2][0] = 0.0; kG[2][1] = 0.0; kG[2][2] = nAC;

				ke = kA*kG;

				psolver->AssembleStiffness2(lmi, lmj, ke);
			}
		}

		// loop over all master nodes
		for (int B=0; B<NM; ++B)
		{
			FENode& nodeB = m_ms.Node(B);
			lmi[0] = nodeB.m_ID[0];
			lmi[1] = nodeB.m_ID[1];
			lmi[2] = nodeB.m_ID[2];

			double nAB = -m_n2[A][B];

			kA[0][0] = m_eps*nAB*(nuA.x*nuA.x); kA[0][1] = m_eps*nAB*(nuA.x*nuA.y); kA[0][2] = m_eps*nAB*(nuA.x*nuA.z);
			kA[1][0] = m_eps*nAB*(nuA.y*nuA.x); kA[1][1] = m_eps*nAB*(nuA.y*nuA.y); kA[1][2] = m_eps*nAB*(nuA.y*nuA.z);
			kA[2][0] = m_eps*nAB*(nuA.z*nuA.x); kA[2][1] = m_eps*nAB*(nuA.z*nuA.y); kA[2][2] = m_eps*nAB*(nuA.z*nuA.z);

			// loop over slave nodes
			for (int C=0; C<NS; ++C)
			{
				FENode& nodeC = m_ss.Node(C);
				lmj[0] = nodeC.m_ID[0];
				lmj[1] = nodeC.m_ID[1];
				lmj[2] = nodeC.m_ID[2];

				double nAC = m_n1[A][C];

				kG[0][0] = nAC; kG[0][1] = 0.0; kG[0][2] = 0.0;
				kG[1][0] = 0.0; kG[1][1] = nAC; kG[1][2] = 0.0;
				kG[2][0] = 0.0; kG[2][1] = 0.0; kG[2][2] = nAC;

				ke = kA*kG;

				psolver->AssembleStiffness2(lmi, lmj, ke);
			}

			// loop over master nodes
			for (int C=0; C<NM; ++C)
			{
				FENode& nodeC = m_ms.Node(C);
				lmj[0] = nodeC.m_ID[0];
				lmj[1] = nodeC.m_ID[1];
				lmj[2] = nodeC.m_ID[2];

				double nAC = -m_n2[A][C];

				kG[0][0] = nAC; kG[0][1] = 0.0; kG[0][2] = 0.0;
				kG[1][0] = 0.0; kG[1][1] = nAC; kG[1][2] = 0.0;
				kG[2][0] = 0.0; kG[2][1] = 0.0; kG[2][2] = nAC;

				ke = kA*kG;

				psolver->AssembleStiffness2(lmi, lmj, ke);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! calculate Lagrangian augmentations
bool FEMortarSlidingContact::Augment(int naug)
{
	return true;
}

//-----------------------------------------------------------------------------
//! update interface data
void FEMortarSlidingContact::Update(int niter)
{
	m_ss.UpdateNormals();
	UpdateMortarWeights();
	UpdateNodalGaps();
}

//-----------------------------------------------------------------------------
//! serialize data to archive
void FEMortarSlidingContact::Serialize(DumpFile& ar)
{
}

//-----------------------------------------------------------------------------
//! shallow copy
void FEMortarSlidingContact::ShallowCopy(DumpStream& dmp, bool bsave)
{
}
