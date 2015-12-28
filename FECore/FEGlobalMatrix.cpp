#include "stdafx.h"
#include "FEGlobalMatrix.h"
#include "FEModel.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
//! Takes a SparseMatrix structure that defines the structure of the global matrix.
FEGlobalMatrix::FEGlobalMatrix(SparseMatrix* pK)
{
	m_pA = pK;
	m_LM.resize(MAX_LM_SIZE);
	m_pMP = 0;
	m_nlm = 0;
}

//-----------------------------------------------------------------------------
//! Deletes the SparseMatrix variable.
FEGlobalMatrix::~FEGlobalMatrix()
{
	delete m_pA;
	m_pA = 0;
	if (m_pMP) delete m_pMP;
}

//-----------------------------------------------------------------------------
void FEGlobalMatrix::Clear()
{ 
	if (m_pA) m_pA->Clear(); 
}

//-----------------------------------------------------------------------------
//! Start building the profile. That is delete the old profile (if there was one)
//! and create a new one. 
void FEGlobalMatrix::build_begin(int neq)
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
void FEGlobalMatrix::build_add(vector<int>& lm)
{
	m_LM[m_nlm++] = lm;
	if (m_nlm >= MAX_LM_SIZE) build_flush();
}

//-----------------------------------------------------------------------------
//! Flush the LM array. The LM array stores a buffer of elements that have to be
//! added to the profile. When this buffer is full it needs to be flushed. This
//! flushin operation causes the actual update of the matrix profile.
void FEGlobalMatrix::build_flush()
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
void FEGlobalMatrix::build_end()
{
	if (m_nlm > 0) build_flush();
	m_pA->Create(*m_pMP);
}

//-----------------------------------------------------------------------------
bool FEGlobalMatrix::Create(FEModel* pfem, int neq, bool breset)
{
    // get nodal DOFS
	FEModel& fem = *pfem;
    DOFS& fedofs = fem.GetDOFS();
    int MAX_NDOFS = fedofs.GetTotalDOFS();

	// keep a pointer to the FEM object
	FEAnalysis* pstep = fem.GetCurrentStep();
	FEMesh& mesh = fem.GetMesh();
	FERigidSystem& rigid = *fem.GetRigidSystem();

	// The first time we come here we build the "static" profile.
	// This static profile stores the contribution to the matrix profile
	// of the "elements" that do not change. Most elements are static except
	// for instance contact elements which can change connectivity in between
	// calls to the Create() function. Storing the static profile instead of
	// reconstructing it every time we come here saves us a lot of time. The 
	// static profile is stored in the variable m_MPs.

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
				d.BuildMatrixProfile(*this);
			}

			// Add rigid bodies to the profile
			if (rigid.Objects())
			{
				vector<int> lm(6);
				int nrb = rigid.Objects();
				for (int i=0; i<nrb; ++i)
				{
					FERigidBody& rb = *rigid.Object(i);
					for (int j=0; j<6; ++j) lm[j] = rb.m_LM[j];
					build_add(lm);
				}
			}

			// Add linear constraints to the profile
			// TODO: we need to add a function build_add(lmi, lmj) for
			// this type of "elements". Now we allocate too much memory
			if (fem.m_LinC.size() > 0)
			{
				int nlin = (int)fem.m_LinC.size();
				vector<int> lm, elm;
				
				// do the cross-term
				// TODO: I have to make this easier. For instance,
				// keep a list that stores for each node the list of
				// elements connected to that node.
				// loop over all solid elements
				for (int nd=0; nd<pstep->Domains(); ++nd)
				{
					FEDomain& dom = *pstep->Domain(nd);
					for (int i=0; i<dom.Elements(); ++i)
					{
						FEElement& el = dom.ElementRef(i);
						dom.UnpackLM(el, elm);
						int ne = (int)elm.size();

						// see if this element connects to the 
						// master node of a linear constraint ...
						int m = el.Nodes();
						for (int j=0; j<m; ++j)
						{
							for (int k=0; k<MAX_NDOFS; ++k)
							{
								int n = fem.m_LCT[el.m_node[j]*MAX_NDOFS + k];
								if (n >= 0)
								{
									// ... it does so we need to connect the 
									// element to the linear constraint
									FELinearConstraint* plc = fem.m_LCA[n];

									int ns = (int)plc->slave.size();

									lm.resize(ne + ns);
									for (int l=0; l<ne; ++l) lm[l] = elm[l];
									
									list<FELinearConstraint::SlaveDOF>::iterator is = plc->slave.begin();
									for (int l=ne; l<ne+ns; ++l, ++is) lm[l] = is->neq;

									build_add(lm);
									break;
								}
							}
						}
					}
				}

				// TODO: do the same thing for shell elements

				// do the constraint term
				int ni;
				list<FELinearConstraint>::iterator ic = fem.m_LinC.begin();
				int n = 0;
				for (int i=0; i<nlin; ++i, ++ic) n += ic->slave.size();
				lm.resize(n);
				ic = fem.m_LinC.begin();
				n = 0;
				for (int i=0; i<nlin; ++i, ++ic)
				{
					ni = (int)ic->slave.size();
					list<FELinearConstraint::SlaveDOF>::iterator is = ic->slave.begin();
					for (int j=0; j<ni; ++j, ++is) lm[n++] = is->neq;
				}
				build_add(lm);
			}

			// do the nonlinear constraints
			int M = fem.NonlinearConstraints();
			for (int m=0; m<M; ++m)
			{
				FENLConstraint* pnlc = fem.NonlinearConstraint(m);
				if (pnlc->IsActive()) pnlc->BuildMatrixProfile(*this);
			}

			// copy the static profile to the MP object
			// Make sure the LM buffer is flushed first.
			build_flush();
			m_MPs = *m_pMP;
		}
		else // --> if (breset)
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
			for (int i=0; i<fem.SurfacePairInteractions(); ++i)
			{
				FESurfacePairInteraction* pci = fem.SurfacePairInteraction(i);
				if (pci->IsActive()) pci->BuildMatrixProfile(*this);
			}
		}
	}
	// All done! We can now finish building the profile and create 
	// the actual sparse matrix. This is done in the following function
	build_end();

	return true;
}

//-----------------------------------------------------------------------------
//! Constructs the stiffness matrix from a FEMesh object. 
bool FEGlobalMatrix::Create(FEMesh& mesh, int neq)
{
	// begin building the profile
	build_begin(neq);
	{
		// Add all elements to the profile
		// Loop over all active domains
		for (int nd=0; nd<mesh.Domains(); ++nd)
		{
			FEDomain& d = mesh.Domain(nd);
			d.BuildMatrixProfile(*this);
		}
	}
	// All done! We can now finish building the profile and create 
	// the actual sparse matrix. This is done in the following function
	build_end();

	return true;
}
