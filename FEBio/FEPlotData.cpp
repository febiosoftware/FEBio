#include "stdafx.h"
#include "FEPlotData.h"
#include "fem.h"
#include "FESlidingInterface2.h"
#include "FEFacet2FacetSliding.h"

//-----------------------------------------------------------------------------
int FEPlotData::VarSize(Var_Type t)
{
	int ndata = 0;
	switch (DataType())
	{
	case FLOAT: ndata = sizeof(float); break;
	case VEC3F: ndata = 3*sizeof(float); break;
	case MAT3FS: ndata = 6*sizeof(float); break;
	}
	assert(ndata);
	return ndata;
}

//-----------------------------------------------------------------------------
void FENodeData::Save(FEM &fem, Archive& ar)
{
	// loop over all node sets
	// write now there is only one, namely the master node set
	// so we just pass the mesh
	int ndata = VarSize(DataType());

	int N = fem.m_mesh.Nodes();
	ar.BeginChunk(0, N*ndata);
	{
		Save(fem.m_mesh, ar);
	}
	ar.EndChunk();
}

//-----------------------------------------------------------------------------
void FEElementData::Save(FEM &fem, Archive& ar)
{
	int ndata = VarSize(DataType());

	// loop over all domains
	FEMesh& m = fem.m_mesh;
	int ND = m.Domains();
	for (int i=0; i<ND; ++i)
	{
		FEDomain& D = m.Domain(i);
		ar.BeginChunk(i+1, ndata*D.Elements());
		{
			Save(D, ar);
		}
		ar.EndChunk();
	}
}

//-----------------------------------------------------------------------------
void FEFaceData::Save(FEM &fem, Archive& ar)
{
	int ndata = VarSize(DataType());

	// loop over all surfaces
	FEMesh& m = fem.m_mesh;
	int NS = m.Surfaces();
	for (int i=0; i<NS; ++i)
	{
		FESurface& S = m.Surface(i);
		ar.BeginChunk(i+1, ndata*m.Elements());
		{
			Save(S, ar);
		}
		ar.EndChunk();
	}
}

//=============================================================================
//                            N O D E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Store the nodal displacements
void FEPlotNodeDisplacement::Save(FEMesh& m, Archive& ar)
{
	float xf[3];
	for (int i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		xf[0] = (float) (node.m_rt.x - node.m_r0.x);
		xf[1] = (float) (node.m_rt.y - node.m_r0.y);
		xf[2] = (float) (node.m_rt.z - node.m_r0.z);

		ar << xf[0] << xf[1] << xf[2];
	}
}

//-----------------------------------------------------------------------------
void FEPlotNodeVelocity::Save(FEMesh& m, Archive& ar)
{
	float xf[3];
	for (int i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		xf[0] = (float) node.m_vt.x;
		xf[1] = (float) node.m_vt.y;
		xf[2] = (float) node.m_vt.z;

		ar << xf[0] << xf[1] << xf[2];
	}
}

//-----------------------------------------------------------------------------
void FEPlotNodeAcceleration::Save(FEMesh& m, Archive& ar)
{
	float xf[3];
	for (int i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		xf[0] = (float) node.m_at.x;
		xf[1] = (float) node.m_at.y;
		xf[2] = (float) node.m_at.z;

		ar << xf[0] << xf[1] << xf[2];
	}
}

//-----------------------------------------------------------------------------
void FEPlotFluidPressure::Save(FEMesh &m, Archive& ar)
{
	float t;
	for (int i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);
		t = (float) node.m_pt;
		ar << t;
	}
}

//=============================================================================
//                           E L E M E N T   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
void FEPlotElementStress::Save(FEDomain& dom, Archive& ar)
{
	int i, j;

	// write solid element data
	float s[6] = {0};
	double f;
	int nint;
	FESolidDomain* pbd = dynamic_cast<FESolidDomain*>(&dom);
	if (pbd)
	{
		for (i=0; i<pbd->Elements(); ++i)
		{
			FESolidElement& el = pbd->Element(i);

			for (j=0; j<6; ++j) s[j] = 0;

			nint = el.GaussPoints();

			f = 1.0 / (double) nint;

			// since the PLOT file requires floats we need to convert
			// the doubles to single precision
			// we output the average stress values of the gauss points
			for (j=0; j<nint; ++j)
			{
				FEElasticMaterialPoint* ppt = (el.m_State[j]->ExtractData<FEElasticMaterialPoint>());

				if (ppt)
				{
					FEElasticMaterialPoint& pt = *ppt;
					s[0] += (float) (f*pt.s.xx());
					s[1] += (float) (f*pt.s.yy());
					s[2] += (float) (f*pt.s.zz());
					s[3] += (float) (f*pt.s.xy());
					s[4] += (float) (f*pt.s.yz());
					s[5] += (float) (f*pt.s.xz());
				}
			}

			ar << s[0] << s[1] << s[2] << s[3] << s[4] << s[5];
		}
	}
}


//-----------------------------------------------------------------------------
void FEPlotFluidFlux::Save(FEDomain &dom, Archive& ar)
{
	int i, j;
	float af[3];
	vec3d ew;
	FESolidDomain* pbd = dynamic_cast<FESolidDomain*>(&dom);
	if (pbd)
	{
		for (i=0; i<pbd->Elements(); ++i)
		{
			FESolidElement& el = pbd->Element(i);

			// calculate average flux
			ew = vec3d(0,0,0);
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FEPoroElasticMaterialPoint* pt = (mp.ExtractData<FEPoroElasticMaterialPoint>());

				if (pt) ew += pt->m_w;
			}

			ew /= el.GaussPoints();

			af[0] = (float) ew.x;
			af[1] = (float) ew.y;
			af[2] = (float) ew.z;

			ar << af[0] << af[1] << af[2];
		}
	}
}

//-----------------------------------------------------------------------------
void FEPlotFiberVector::Save(FEDomain &dom, Archive& ar)
{
	int i, j, n;
	float a[3];
	vec3d r;
	FESolidDomain* pbd = dynamic_cast<FESolidDomain*>(&dom);
	if (pbd)
	{
		int BE = pbd->Elements();
		for (i=0; i<BE; ++i)
		{
			FESolidElement& el = pbd->Element(i);
			n = el.GaussPoints();
			r = vec3d(0,0,0);
			for (j=0; j<n; ++j)
			{
				FEElasticMaterialPoint& pt = *el.m_State[j]->ExtractData<FEElasticMaterialPoint>();
				r.x += pt.Q[0][0];
				r.y += pt.Q[1][0];
				r.z += pt.Q[2][0];
			}
			r /= n;
			a[0] = (float) r.x;
			a[1] = (float) r.y;
			a[2] = (float) r.z;

			ar << a[0] << a[1] << a[2];
		}
	}
}

//=============================================================================
//                           S U R F A C E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
void FEPlotContactGap::Save(FESurface& surf, Archive& ar)
{
/*
	FEMesh& mesh = fem.m_mesh;

	vector<float> t(mesh.Nodes());
	zero(t);

	int i, j;

	for (i=0; i<fem.ContactInterfaces(); ++i)
	{
		FESlidingInterface* psi = dynamic_cast<FESlidingInterface*>(fem.m_CI[i]);
		if (psi)
		{
			FESlidingSurface& ms = psi->m_ms;
			FESlidingSurface& ss = psi->m_ss;

			for (j=0; j<ms.Nodes(); ++j) t[ms.node[j]] += (float) ms.gap[j];
			for (j=0; j<ss.Nodes(); ++j) t[ss.node[j]] += (float) ss.gap[j];
		}

		FETiedInterface* pti = dynamic_cast<FETiedInterface*>(fem.m_CI[i]);
		if (pti)
		{
			FETiedContactSurface& ms = pti->ms;
			FETiedContactSurface& ss = pti->ss;

			for (j=0; j<ms.Nodes(); ++j) t[ms.node[j]] += (float) ms.gap[j].norm();
			for (j=0; j<ss.Nodes(); ++j) t[ss.node[j]] += (float) ss.gap[j].norm();
		}

		FERigidWallInterface* pri = dynamic_cast<FERigidWallInterface*>(fem.m_CI[i]);
		if (pri)
		{
			FERigidWallSurface& ss = pri->m_ss;
			for (j=0; j<ss.Nodes(); ++j) t[ss.node[j]] += (float) ss.gap[j];
		}

		FEFacet2FacetSliding* pf = dynamic_cast<FEFacet2FacetSliding*>(fem.m_CI[i]);
		if (pf)
		{
			vector<int> val(fem.m_mesh.Nodes()); zero(val);
			double gi[4], gn[4];
			int ni, ne, n, k;

			for (n=0; n<pf->m_npass; ++n)
			{
				FEFacetSlidingSurface& s = (n==0?pf->m_ss:pf->m_ms);

				int nint = 0;
				for (j=0; j<s.Elements(); ++j)
				{
					FESurfaceElement& el = s.Element(j);
					ne = el.Nodes();
					ni = el.GaussPoints();
					for (k=0; k<ni; ++k, ++nint)
					{
						gi[k] = s.m_gap[nint];
					}

					el.project_to_nodes(gi, gn);

					for (k=0; k<ne; ++k)
					{
						int m = el.m_node[k];
						t[m] += (float) gn[k];
						val[m]++;
					}
				}
			}

			for (j=0; j<fem.m_mesh.Nodes(); ++j) if (val[j] > 1) t[j] /= (float) val[j];
		}

		FESlidingInterface2* ps2 = dynamic_cast<FESlidingInterface2*>(fem.m_CI[i]);
		if (ps2)
		{
			vector<int> val(fem.m_mesh.Nodes()); zero(val);
			double gi[4], gn[4];
			int ni, ne, n, k;

			for (n=0; n<ps2->m_npass; ++n)
			{
				FESlidingSurface2& s = (n==0?ps2->m_ss:ps2->m_ms);

				int nint = 0;
				for (j=0; j<s.Elements(); ++j)
				{
					FESurfaceElement& el = s.Element(j);
					ne = el.Nodes();
					ni = el.GaussPoints();
					for (k=0; k<ni; ++k, ++nint)
					{
						gi[k] = s.m_gap[nint];
					}

					el.project_to_nodes(gi, gn);

					for (k=0; k<ne; ++k)
					{
						int m = el.m_node[k];
						t[m] += (float) gn[k];
						val[m]++;
					}
				}
			}

			for (j=0; j<fem.m_mesh.Nodes(); ++j) if (val[j] > 1) t[j] /= (float) val[j];
		}
	}

	// store the data to file
	// Note that we only save gap values of nodes that are actually in contact
	for (i=0; i<mesh.Nodes(); ++i) t[i] = (t[i]<0? 0.f : t[i]);
	fwrite(&t[0], sizeof(float), mesh.Nodes(), fp);
*/
}

//-----------------------------------------------------------------------------
void FEPlotContactTraction::Save(FESurface &surf, Archive& ar)
{
/*
	int i, j, k, n;

	vector<float> acc(3*fem.m_mesh.Nodes()); zero(acc);
	for (i=0; i<(int) fem.m_CI.size(); ++i)
	{
		FESlidingInterface* psi = dynamic_cast<FESlidingInterface*> (fem.m_CI[i]);
		if (psi)
		{
			for (n=0; n<psi->m_npass; ++n)
			{
				FESlidingSurface& ss = (n==0?psi->m_ss:psi->m_ms);
				FESlidingSurface& ms = (n==0?psi->m_ms:psi->m_ss);
				for (j=0; j<ss.Nodes(); ++j)
				{
					int m = ss.node[j];
					vec3d t = ss.traction(j);

					acc[3*m  ] += (float) t.x;
					acc[3*m+1] += (float) t.y;
					acc[3*m+2] += (float) t.z;
				}
			}
		}

		FEFacet2FacetSliding* pf = dynamic_cast<FEFacet2FacetSliding*>(fem.m_CI[i]);
		if (pf)
		{
			vector<int> val(fem.m_mesh.Nodes()); zero(val);
			double ti[4], tn[4], gi[4], gn[4], li[4], ln[4];
			int ni, ne;

			for (n=0; n<pf->m_npass; ++n)
			{
				FEFacetSlidingSurface& s = (n==0?pf->m_ss:pf->m_ms);

				int nint = 0;
				for (j=0; j<s.Elements(); ++j)
				{
					FESurfaceElement& el = s.Element(j);
					ne = el.Nodes();
					ni = el.GaussPoints();
					for (k=0; k<ni; ++k, ++nint)
					{
						li[k] = s.m_Lm[nint];
						gi[k] = s.m_gap[nint];
						ti[k] = li[k] + pf->m_epsn*gi[k];

						gi[k] = (gi[k]>=0?gi[k] : 0);
						ti[k] = (ti[k]>=0?ti[k] : 0);
					}

					el.project_to_nodes(li, ln);
					el.project_to_nodes(gi, gn);
					el.project_to_nodes(ti, tn);

					for (k=0; k<ne; ++k)
					{
						int m = el.m_node[k];
						acc[3*m  ] += (float) (ln[k]>=0?ln[k]:0);
						acc[3*m+1] += (float) (gn[k]>=0?gn[k]:0);
						acc[3*m+2] += (float) (tn[k]>=0?tn[k]:0);
						val[m]++;
					}
				}
			}

			for (j=0; j<fem.m_mesh.Nodes(); ++j) if (val[j] > 1) 
			{ 
				acc[3*j  ] /= (float) val[j]; 
				acc[3*j+1] /= (float) val[j]; 
				acc[3*j+2] /= (float) val[j]; 
			}
		}

		FESlidingInterface2* ps2 = dynamic_cast<FESlidingInterface2*>(fem.m_CI[i]);
		if (ps2)
		{
			vector<int> val(fem.m_mesh.Nodes()); zero(val);
			double ti[4], tn[4], gi[4], gn[4], li[4], ln[4];
			int ni, ne;

			for (n=0; n<ps2->m_npass; ++n)
			{
				FESlidingSurface2& s = (n==0?ps2->m_ss:ps2->m_ms);

				int nint = 0;
				for (j=0; j<s.Elements(); ++j)
				{
					FESurfaceElement& el = s.Element(j);
					ne = el.Nodes();
					ni = el.GaussPoints();
					for (k=0; k<ni; ++k, ++nint)
					{
						li[k] = s.m_Lmd[nint];
						gi[k] = s.m_gap[nint];
						ti[k] = li[k] + ps2->m_epsn*gi[k];

						gi[k] = (gi[k]>=0?gi[k] : 0);
						ti[k] = (ti[k]>=0?ti[k] : 0);
					}

					el.project_to_nodes(li, ln);
					el.project_to_nodes(gi, gn);
					el.project_to_nodes(ti, tn);

					for (k=0; k<ne; ++k)
					{
						int m = el.m_node[k];
						acc[3*m  ] += (float) (ln[k]>=0?ln[k]:0);
						acc[3*m+1] += (float) (gn[k]>=0?gn[k]:0);
						acc[3*m+2] += (float) (tn[k]>=0?tn[k]:0);
						val[m]++;
					}
				}
			}

			for (j=0; j<fem.m_mesh.Nodes(); ++j) if (val[j] > 1) 
			{ 
				acc[3*j  ] /= (float) val[j]; 
				acc[3*j+1] /= (float) val[j]; 
				acc[3*j+2] /= (float) val[j]; 
			}
		}
	}

	fwrite(&acc[0], sizeof(float)*3, fem.m_mesh.Nodes(), fp);
*/
}
