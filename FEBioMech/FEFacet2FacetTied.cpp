#include "stdafx.h"
#include "FEFacet2FacetTied.h"
#include "FECore/FEModel.h"
#include "FECore/log.h"

//-----------------------------------------------------------------------------
// Define tied interface parameters
BEGIN_PARAMETER_LIST(FEFacet2FacetTied, FEContactInterface)
	ADD_PARAMETER(m_blaugon, FE_PARAM_BOOL  , "laugon"          ); 
	ADD_PARAMETER(m_atol   , FE_PARAM_DOUBLE, "tolerance"       );
	ADD_PARAMETER(m_eps    , FE_PARAM_DOUBLE, "penalty"         );
	ADD_PARAMETER(m_naugmin, FE_PARAM_INT   , "minaug"          );
	ADD_PARAMETER(m_naugmax, FE_PARAM_INT   , "maxaug"          );
	ADD_PARAMETER(m_stol   , FE_PARAM_DOUBLE, "search_tolerance");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEFacetTiedSurface::Data::Data()
{
	m_gap = vec3d(0,0,0);
	m_Lm = vec3d(0,0,0);
	m_rs = vec2d(0,0);
	m_pme = (FESurfaceElement*) 0;
}

//-----------------------------------------------------------------------------
FEFacetTiedSurface::FEFacetTiedSurface(FEMesh* pm) : FEContactSurface(pm)
{

}

//-----------------------------------------------------------------------------
bool FEFacetTiedSurface::Init()
{
	// initialize surface data first
	if (FEContactSurface::Init() == false) return false;

	// allocate data structures
	const int NE = Elements();
	m_Data.resize(NE);
	for (int i=0; i<NE; ++i)
	{
		FESurfaceElement& el = Element(i);
		int nint = el.GaussPoints();
		m_Data[i].resize(nint);
	}

	return true;
}

//-----------------------------------------------------------------------------
//! \todo Originally, we only copied Lmd, gap, Ln and reset pme to zero.
//!       Need to check if this achieves the same
void FEFacetTiedSurface::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (bsave)
	{
		for (int i=0; i<(int) m_Data.size(); ++i)
		{
			vector<Data>& di = m_Data[i];
			int nint = (int) di.size();
			for (int j=0; j<nint; ++j)
			{
				Data& d = di[j];
				dmp << d.m_gap;
				dmp << d.m_rs;
				dmp << d.m_Lm;
			}
		}
	}
	else
	{
		for (int i=0; i<(int) m_Data.size(); ++i)
		{
			vector<Data>& di = m_Data[i];
			int nint = (int) di.size();
			for (int j=0; j<nint; ++j)
			{
				Data& d = di[j];
				dmp >> d.m_gap;
				dmp >> d.m_rs;
				dmp >> d.m_Lm;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEFacetTiedSurface::Serialize(DumpFile &ar)
{
	FEContactSurface::Serialize(ar);
	if (ar.IsSaving())
	{
		for (int i=0; i<(int) m_Data.size(); ++i)
		{
			vector<Data>& di = m_Data[i];
			int nint = (int) di.size();
			for (int j=0; j<nint; ++j)
			{
				Data& d = di[j];
				ar << d.m_gap;
				ar << d.m_rs;
				ar << d.m_Lm;
			}
		}
	}
	else
	{
		for (int i=0; i<(int) m_Data.size(); ++i)
		{
			vector<Data>& di = m_Data[i];
			int nint = (int) di.size();
			for (int j=0; j<nint; ++j)
			{
				Data& d = di[j];
				ar >> d.m_gap;
				ar >> d.m_rs;
				ar >> d.m_Lm;
			}
		}
	}
}

//=============================================================================
FEFacet2FacetTied::FEFacet2FacetTied(FEModel* pfem) : FEContactInterface(pfem), m_ss(&pfem->GetMesh()), m_ms(&pfem->GetMesh())
{
	static int count = 1;
	m_ntype = FE_FACET2FACET_TIED;

	// define sibling relationships
	m_ss.SetSibling(&m_ms);
	m_ms.SetSibling(&m_ss);

	// initial parameter values
	m_blaugon = false;
	m_atol    = 0.01;
	m_eps     = 1.0;
	m_naugmin = 0;
	m_naugmax = 10;
	m_stol    = 0.0001;

	// give this interface an ID (TODO: where is this actually used?)
	m_nID = count++;
}

//-----------------------------------------------------------------------------
//! Initialization. This function intializes the surfaces data
bool FEFacet2FacetTied::Init()
{
	// create the surfaces
	if (m_ss.Init() == false) return false;
	if (m_ms.Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
void FEFacet2FacetTied::ShallowCopy(DumpStream& dmp, bool bsave)
{
	m_ss.ShallowCopy(dmp, bsave);
	m_ms.ShallowCopy(dmp, bsave);
}

//-----------------------------------------------------------------------------
//! Interface activation. Also projects slave surface onto master surface
void FEFacet2FacetTied::Activate()
{
	// Don't forget to call base member!
	FEContactInterface::Activate();

	// project slave surface onto master surface
	ProjectSurface(m_ss, m_ms);
}

//-----------------------------------------------------------------------------
void FEFacet2FacetTied::ProjectSurface(FEFacetTiedSurface& ss, FEFacetTiedSurface& ms)
{
	// get the mesh
	FEMesh& mesh = *ss.GetMesh();

	// loop over all slave elements
	for (int i=0; i<ss.Elements(); ++i)
	{
		// get the slave element
		FESurfaceElement& se = ss.Element(i);

		// get nodal coordinates
		int nn = se.Nodes();
		vec3d re[FEElement::MAX_NODES];
		for (int j=0; j<nn; ++j) re[j] = mesh.Node(se.m_node[j]).m_rt;

		// loop over all its integration points
		int nint = se.GaussPoints();
		for (int j=0; j<nint; ++j)
		{
			// get integration point data
			FEFacetTiedSurface::Data& pt = ss.m_Data[i][j];

			// calculate the global coordinates of this integration point
			vec3d x = se.eval(re, j);

			// find the master element
			vec3d q; vec2d rs;
			FESurfaceElement* pme = ms.ClosestPointProjection(x, q, rs, (i==0), m_stol);
			if (pme)
			{
				// store the master element
				pt.m_pme = pme;
				pt.m_rs[0] = rs[0];
				pt.m_rs[1] = rs[1];

				// calculate gap
				pt.m_gap = x - q;
			}
			else pt.m_pme = 0;
		}
	}
}

//-----------------------------------------------------------------------------
//! Update tied interface data. This function re-evaluates the gaps between
//! the slave node and their projections onto the master surface.
//!
void FEFacet2FacetTied::Update(int niter)
{
	// get the mesh
	FEMesh& mesh = *m_ss.GetMesh();

	// loop over all slave elements
	const int NE = m_ss.Elements();
	for (int i=0; i<NE; ++i)
	{
		// next element
		FESurfaceElement& se = m_ss.Element(i);
		int nseln = se.Nodes();

		// get the nodal coordinates
		vec3d rs[FEElement::MAX_NODES];
		for (int j=0; j<nseln; ++j) rs[j] = mesh.Node(se.m_node[j]).m_rt;

		// loop over all integration points
		const int nint = se.GaussPoints();
		for (int n=0; n<nint; ++n)
		{
			// get integration point data
			FEFacetTiedSurface::Data& pt = m_ss.m_Data[i][n];

			FESurfaceElement* pme = pt.m_pme;
			if (pme)
			{
				FESurfaceElement& me = static_cast<FESurfaceElement&>(*pme);

				// get the current slave nodal position
				vec3d rn = se.eval(rs, n);

				// get the natural coordinates of the slave projection
				double r = pt.m_rs[0];
				double s = pt.m_rs[1];

				// get the master nodal coordinates
				int nmeln = me.Nodes();
				vec3d y[FEElement::MAX_NODES];
				for (int l=0; l<nmeln; ++l) y[l] = mesh.Node( me.m_node[l] ).m_rt;

				// calculate the slave node projection
				vec3d q = me.eval(y, r, s);

				// calculate the gap function
				pt.m_gap = rn - q;
			}
		}
	}
}


//-----------------------------------------------------------------------------
//! This function calculates the contact forces for a tied interface.
void FEFacet2FacetTied::ContactForces(FEGlobalVector& R)
{
	vector<int> sLM, mLM, LM, en;
	vector<double> fe;

	// shape functions
	double Hm[FEElement::MAX_NODES];

	// get the mesh
	FEMesh& mesh = *m_ss.GetMesh();

	// loop over all elements
	const int NE = m_ss.Elements();
	for (int i=0; i<NE; ++i)
	{
		// get the next element
		FESurfaceElement& se = m_ss.Element(i);
		int nseln = se.Nodes();

		// integration weights
		double* w = se.GaussWeights();

		// get the element's LM vector
		m_ss.UnpackLM(se, sLM);

		// loop over integration points
		const int nint = se.GaussPoints();
		for (int n=0; n<nint; ++n)
		{
			// get integration point data
			FEFacetTiedSurface::Data& pt = m_ss.m_Data[i][n];

			// get the master element
			FESurfaceElement* pme = pt.m_pme;
			if (pme)
			{
				// get the master element
				FESurfaceElement& me = dynamic_cast<FESurfaceElement&>(*pme);
				m_ms.UnpackLM(me, mLM);
				int nmeln = me.Nodes();

				// get slave contact force
				vec3d tc = pt.m_Lm + pt.m_gap*m_eps;

				// calculate jacobian
				// note that we are integrating over the reference surface
				double detJ = m_ss.jac0(se, n);

				// slave shape functions
				double* Hs = se.H(n);

				// master shape functions
				double r = pt.m_rs[0];
				double s = pt.m_rs[1];
				me.shape_fnc(Hm, r, s);

				// calculate degrees of freedom
				int ndof = 3*(nseln + nmeln);

				// calculate the force vector
				fe.resize(ndof);
				for (int k=0; k<nseln; ++k)
				{
					fe[3*k  ] = -detJ*w[n]*tc.x*Hs[k];
					fe[3*k+1] = -detJ*w[n]*tc.y*Hs[k];
					fe[3*k+2] = -detJ*w[n]*tc.z*Hs[k];
				}
				for (int k=0; k<nmeln; ++k)
				{
					fe[3*(k+nseln)  ] = detJ*w[n]*tc.x*Hm[k];
					fe[3*(k+nseln)+1] = detJ*w[n]*tc.y*Hm[k];
					fe[3*(k+nseln)+2] = detJ*w[n]*tc.z*Hm[k];
				}

				// build the LM vector
				LM.resize(ndof);
				for (int k=0; k<nseln; ++k)
				{
					LM[3*k  ] = sLM[3*k  ];
					LM[3*k+1] = sLM[3*k+1];
					LM[3*k+2] = sLM[3*k+2];
				}

				for (int k=0; k<nmeln; ++k)
				{
					LM[3*(k+nseln)  ] = mLM[3*k  ];
					LM[3*(k+nseln)+1] = mLM[3*k+1];
					LM[3*(k+nseln)+2] = mLM[3*k+2];
				}

				// build the en vector
				en.resize(nseln+nmeln);
				for (int k=0; k<nseln; ++k) en[k      ] = se.m_node[k];
				for (int k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];

				// assemble the global residual
				R.Assemble(en, LM, fe);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculate the stiffness matrix contribution.
void FEFacet2FacetTied::ContactStiffness(FESolver* psolver)
{
	vector<int> sLM, mLM, LM, en;
	matrix ke;

	// shape functions
	double Hm[FEElement::MAX_NODES];

	// loop over all slave elements
	const int NE = m_ss.Elements();
	for (int i=0; i<NE; ++i)
	{
		// get the next element
		FESurfaceElement& se = m_ss.Element(i);
		int nseln = se.Nodes();

		// get the element's LM vector
		m_ss.UnpackLM(se, sLM);

		// integration weights
		double* w = se.GaussWeights();

		// loop over all integration points
		const int nint = se.GaussPoints();
		for (int n=0; n<nint; ++n)
		{
			// get intgration point data
			FEFacetTiedSurface::Data& pt = m_ss.m_Data[i][n];

			// get the master element
			FESurfaceElement* pme = pt.m_pme;
			if (pme)
			{
				// get the master element
				FESurfaceElement& me = dynamic_cast<FESurfaceElement&> (*pme);
				int nmeln = me.Nodes();
				m_ms.UnpackLM(me, mLM);

				// calculate jacobian
				double detJ = m_ss.jac0(se, n);

				// slave shape functions
				double* Hs = se.H(n);

				// master shape functions
				double r = pt.m_rs[0];
				double s = pt.m_rs[1];
				me.shape_fnc(Hm, r, s);

				// calculate degrees of freedom
				int ndof = 3*(nseln + nmeln);

				// create the stiffness matrix
				ke.resize(ndof, ndof);
				ke.zero();
				for (int k=0; k<nseln; ++k)
				{
					for (int l=0; l<nseln; ++l)
					{
						ke[3*k  ][3*l  ] = Hs[k]*Hs[l];
						ke[3*k+1][3*l+1] = Hs[k]*Hs[l];
						ke[3*k+2][3*l+2] = Hs[k]*Hs[l];
					}
				}

				for (int k=0; k<nseln; ++k)
				{
					for (int l=0; l<nmeln; ++l)
					{
						ke[3*k  ][3*(l+nseln)  ] = -Hs[k]*Hm[l];
						ke[3*k+1][3*(l+nseln)+1] = -Hs[k]*Hm[l];
						ke[3*k+2][3*(l+nseln)+2] = -Hs[k]*Hm[l];

						ke[3*(l+nseln)  ][3*k  ] = -Hs[k]*Hm[l];
						ke[3*(l+nseln)+1][3*k+1] = -Hs[k]*Hm[l];
						ke[3*(l+nseln)+2][3*k+2] = -Hs[k]*Hm[l];
					}
				}

				for (int k=0; k<nmeln; ++k)
					for (int l=0; l<nmeln; ++l)
					{
						ke[3*(k+nseln)  ][3*(l+nseln)  ] = Hm[k]*Hm[l];
						ke[3*(k+nseln)+1][3*(l+nseln)+1] = Hm[k]*Hm[l];
						ke[3*(k+nseln)+2][3*(l+nseln)+2] = Hm[k]*Hm[l];
					}

				for (int k=0; k<ndof; ++k)
					for (int l=0; l<ndof; ++l) ke[k][l] *= m_eps*detJ*w[n];

				// build the LM vector
				LM.resize(ndof);
				for (int k=0; k<nseln; ++k)
				{
					LM[3*k  ] = sLM[3*k  ];
					LM[3*k+1] = sLM[3*k+1];
					LM[3*k+2] = sLM[3*k+2];
				}
				for (int k=0; k<nmeln; ++k)
				{
					LM[3*(k+nseln)  ] = mLM[3*k  ];
					LM[3*(k+nseln)+1] = mLM[3*k+1];
					LM[3*(k+nseln)+2] = mLM[3*k+2];
				}

				// build the en vector
				en.resize(nseln+nmeln);
				for (int k=0; k<nseln; ++k) en[k      ] = se.m_node[k];
				for (int k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];

				// assemble the global residual
				psolver->AssembleStiffness(en, LM, ke);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Do an augmentation.
bool FEFacet2FacetTied::Augment(int naug)
{
	// make sure we need to augment
	if (!m_blaugon) return true;

	const int NS = (const int) m_ss.m_Data.size();

	// calculate initial norms
	double normL0 = 0;
	for (int i=0; i<NS; ++i)
	{
		const int nint = (const int) m_ss.m_Data[i].size();
		for (int j=0; j<nint; ++j)
		{
			FEFacetTiedSurface::Data& pt = m_ss.m_Data[i][j];
			vec3d& lm = pt.m_Lm;
			normL0 += lm*lm;
		}
	}
	normL0 = sqrt(normL0);

	// update Lagrange multipliers and calculate current norms
	double normL1 = 0;
	double normgc = 0;
	int N = 0;
	for (int i=0; i<NS; ++i)
	{
		const int nint = (const int) m_ss.m_Data[i].size();
		for (int j=0; j<nint; ++j)
		{
			FEFacetTiedSurface::Data& pt = m_ss.m_Data[i][j];

			vec3d lm = pt.m_Lm + pt.m_gap*m_eps;

			normL1 += lm*lm;
			if (pt.m_pme != 0)
			{
				double g = pt.m_gap.norm();
				normgc += g*g;
				++N;
			}
		}
	}	
	if (N == 0) N=1;
	normL1 = sqrt(normL1);
	normgc = sqrt(normgc / N);

	// check convergence of constraints
	clog.printf(" tied interface # %d\n", m_nID);
	clog.printf("                        CURRENT        REQUIRED\n");
	double pctn = 0;
	if (fabs(normL1) > 1e-10) pctn = fabs((normL1 - normL0)/normL1);
	clog.printf("    normal force : %15le %15le\n", pctn, m_atol);
	clog.printf("    gap function : %15le       ***\n", normgc);
		
	// check convergence
	bool bconv = true;
	if (pctn >= m_atol) bconv = false;
	if (naug < m_naugmin ) bconv = false;
	if (naug >= m_naugmax) bconv = true;

	if (bconv == false) 
	{
		for (int i=0; i<NS; ++i)
		{
			const int nint = (const int) m_ss.m_Data[i].size();
			for (int j=0; j<nint; ++j)
			{
				FEFacetTiedSurface::Data& pt = m_ss.m_Data[i][j];

				// update Lagrange multipliers
				pt.m_Lm = pt.m_Lm + pt.m_gap*m_eps;
			}
		}	
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! Serialize the data to the archive.
void FEFacet2FacetTied::Serialize(DumpFile &ar)
{
	// store contact data
	FEContactInterface::Serialize(ar);

	// store contact surface data
	m_ss.Serialize(ar);
	m_ms.Serialize(ar);
}
