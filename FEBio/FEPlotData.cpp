#include "stdafx.h"
#include "FEPlotData.h"
#include "fem.h"
#include "FESlidingInterface2.h"
#include "FEFacet2FacetSliding.h"

//-----------------------------------------------------------------------------
void FEPlotNodeDisplacement::Save(FEM& fem, Archive& ar)
{
	float xf[3];
	for (int i=0; i<fem.m_mesh.Nodes(); ++i)
	{
		FENode& node = fem.m_mesh.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		xf[0] = (float) (node.m_rt.x - node.m_r0.x);
		xf[1] = (float) (node.m_rt.y - node.m_r0.y);
		xf[2] = (float) (node.m_rt.z - node.m_r0.z);

		ar.write(xf, sizeof(float), 3);
	}
}

//-----------------------------------------------------------------------------
void FEPlotNodeVelocity::Save(FEM& fem, Archive& ar)
{
	float xf[3];
	for (int i=0; i<fem.m_mesh.Nodes(); ++i)
	{
		FENode& node = fem.m_mesh.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		xf[0] = (float) node.m_vt.x;
		xf[1] = (float) node.m_vt.y;
		xf[2] = (float) node.m_vt.z;

		ar.write(xf, sizeof(float), 3);
	}
}

//-----------------------------------------------------------------------------
void FEPlotNodeAcceleration::Save(FEM& fem, Archive& ar)
{
	float xf[3];
	for (int i=0; i<fem.m_mesh.Nodes(); ++i)
	{
		FENode& node = fem.m_mesh.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		xf[0] = (float) node.m_at.x;
		xf[1] = (float) node.m_at.y;
		xf[2] = (float) node.m_at.z;

		ar.write(xf, sizeof(float), 3);
	}
}

//-----------------------------------------------------------------------------
void FEPlotElementStress::Save(FEM& fem, Archive& ar)
{
	int i, j;

	FEMesh& mesh = fem.m_mesh;

	// write solid element data
	float s[6] = {0};
	double f;
	int nint;
	for (int nd=0; nd < mesh.Domains(); ++nd)
	{
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(nd));
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

				ar.write(s, sizeof(float), 6);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEPlotContactGap::Save(FEM &fem, Archive &ar)
{
	FEMesh& mesh = fem.m_mesh;

	vector<float> t(mesh.Nodes());
	t.zero();

	int i, j;

	for (i=0; i<fem.ContactInterfaces(); ++i)
	{
		FESlidingInterface* psi = dynamic_cast<FESlidingInterface*>(&fem.m_CI[i]);
		if (psi)
		{
			FESlidingSurface& ms = psi->m_ms;
			FESlidingSurface& ss = psi->m_ss;

			for (j=0; j<ms.Nodes(); ++j) t[ms.node[j]] += (float) ms.gap[j];
			for (j=0; j<ss.Nodes(); ++j) t[ss.node[j]] += (float) ss.gap[j];
		}

		FETiedInterface* pti = dynamic_cast<FETiedInterface*>(&fem.m_CI[i]);
		if (pti)
		{
			FETiedContactSurface& ms = pti->ms;
			FETiedContactSurface& ss = pti->ss;

			for (j=0; j<ms.Nodes(); ++j) t[ms.node[j]] += (float) ms.gap[j].norm();
			for (j=0; j<ss.Nodes(); ++j) t[ss.node[j]] += (float) ss.gap[j].norm();
		}

		FERigidWallInterface* pri = dynamic_cast<FERigidWallInterface*>(&fem.m_CI[i]);
		if (pri)
		{
			FERigidWallSurface& ss = pri->m_ss;
			for (j=0; j<ss.Nodes(); ++j) t[ss.node[j]] += (float) ss.gap[j];
		}

		FEFacet2FacetSliding* pf = dynamic_cast<FEFacet2FacetSliding*>(&fem.m_CI[i]);
		if (pf)
		{
			vector<int> val(fem.m_mesh.Nodes()); val.zero();
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

		FESlidingInterface2* ps2 = dynamic_cast<FESlidingInterface2*>(&fem.m_CI[i]);
		if (ps2)
		{
			vector<int> val(fem.m_mesh.Nodes()); val.zero();
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
	for (i=0; i<mesh.Nodes(); ++i) ar << (t[i]<0? 0.f : t[i]);
}

//-----------------------------------------------------------------------------
void FEPlotContactTraction::Save(FEM &fem, Archive &ar)
{
	int i, j, k, n;

	vector<float[3]> acc(fem.m_mesh.Nodes());
	for (i=0; i<fem.m_mesh.Nodes(); ++i) acc[i][0] = acc[i][1] = acc[i][2] = 0;
	for (i=0; i<fem.m_CI.size(); ++i)
	{
		FESlidingInterface* psi = dynamic_cast<FESlidingInterface*> (&fem.m_CI[i]);
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

					acc[m][0] += (float) t.x;
					acc[m][1] += (float) t.y;
					acc[m][2] += (float) t.z;
						
//					acc[m][0] = (float) ss.Lt[j][0];
//					acc[m][1] = (float) ss.Lt[j][1];
//					acc[m][2] = (float) ss.Lm[j];
				}
			}
		}

		FEFacet2FacetSliding* pf = dynamic_cast<FEFacet2FacetSliding*>(&fem.m_CI[i]);
		if (pf)
		{
			vector<int> val(fem.m_mesh.Nodes()); val.zero();
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
						acc[m][0] += (float) (ln[k]>=0?ln[k]:0);
						acc[m][1] += (float) (gn[k]>=0?gn[k]:0);
						acc[m][2] += (float) (tn[k]>=0?tn[k]:0);
						val[m]++;
					}
				}
			}

			for (j=0; j<fem.m_mesh.Nodes(); ++j) if (val[j] > 1) 
			{ 
				acc[j][0] /= (float) val[j]; 
				acc[j][1] /= (float) val[j]; 
				acc[j][2] /= (float) val[j]; 
			}
		}

		FESlidingInterface2* ps2 = dynamic_cast<FESlidingInterface2*>(&fem.m_CI[i]);
		if (ps2)
		{
			vector<int> val(fem.m_mesh.Nodes()); val.zero();
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
						acc[m][0] += (float) (ln[k]>=0?ln[k]:0);
						acc[m][1] += (float) (gn[k]>=0?gn[k]:0);
						acc[m][2] += (float) (tn[k]>=0?tn[k]:0);
						val[m]++;
					}
				}
			}

			for (j=0; j<fem.m_mesh.Nodes(); ++j) if (val[j] > 1) 
			{ 
				acc[j][0] /= (float) val[j]; 
				acc[j][1] /= (float) val[j]; 
				acc[j][2] /= (float) val[j]; 
			}
		}
	}

	ar.write(acc, sizeof(float)*3, fem.m_mesh.Nodes() );
}

//-----------------------------------------------------------------------------
void FEPlotFluidPressure::Save(FEM &fem, Archive &ar)
{
	float t;
	for (int i=0; i<fem.m_mesh.Nodes(); ++i)
	{
		FENode& node = fem.m_mesh.Node(i);
		t = (float) node.m_pt;
		ar << t;
	}
}

//-----------------------------------------------------------------------------
void FEPlotFluidFlux::Save(FEM &fem, Archive &ar)
{
	FEMesh& mesh = fem.m_mesh;

	int i, j;

	float af[3];
	vec3d ew;
	for (int nd=0; nd < mesh.Domains(); ++nd)
	{
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(nd));
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

				ar.write(af, sizeof(float), 3);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEPlotFiberVector::Save(FEM &fem, Archive &ar)
{
	int i, j, n;
	FEMesh& mesh = fem.m_mesh;


	float a[3];
	vec3d r;

	for (int nd=0; nd<mesh.Domains(); ++nd)
	{
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(nd));
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
				ar.write(a, sizeof(float), 3);
			}
		}
	}
}
