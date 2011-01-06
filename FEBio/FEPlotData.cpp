#include "stdafx.h"
#include "FEPlotData.h"
#include "fem.h"
#include "FESlidingInterface.h"
#include "FESlidingInterface2.h"
#include "FEFacet2FacetSliding.h"
#include "FEPlotDataFactory.h"

//-----------------------------------------------------------------------------
REGISTER_PLOTDATA(FEPlotNodeDisplacement, "displacement"    );
REGISTER_PLOTDATA(FEPlotNodeVelocity    , "velocity"        );
REGISTER_PLOTDATA(FEPlotNodeAcceleration, "acceleration"    );
REGISTER_PLOTDATA(FEPlotFluidPressure   , "fluid pressure"  );
REGISTER_PLOTDATA(FEPlotElementStress   , "stress"          );
REGISTER_PLOTDATA(FEPlotFluidFlux       , "fluid flux"      );
REGISTER_PLOTDATA(FEPlotFiberVector     , "fiber vector"    );
REGISTER_PLOTDATA(FEPlotContactGap      , "contact gap"     );
REGISTER_PLOTDATA(FEPlotContactTraction , "contact traction");
REGISTER_PLOTDATA(FEPlotShellThickness  , "shell thickness" );

//-----------------------------------------------------------------------------
int FEPlotData::VarSize(Var_Type t)
{
	int ndata = 0;
	switch (DataType())
	{
	case FLOAT: ndata = 1; break;
	case VEC3F: ndata = 3; break;
	case MAT3FS: ndata = 6; break;
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
	vector<float> a; a.reserve(ndata*N);
	if (Save(fem.m_mesh, a))
	{
		assert(a.size() == N*ndata);
		ar.WriteChunk(0, a);
	}
}

//-----------------------------------------------------------------------------
void FEElementData::Save(FEM &fem, Archive& ar)
{
	int ndata = VarSize(DataType());
	int fmt = m_sfmt;

	// loop over all domains
	FEMesh& m = fem.m_mesh;
	int ND = m.Domains();
	for (int i=0; i<ND; ++i)
	{
		FEDomain& D = m.Domain(i);
		int nsize = ndata*D.Elements();
		vector<float> a; 

		if (fmt == FMT_MULT)
		{
			// since all elements have the same type within a domain
			// we just grab the number of nodes of the first element 
			// to figure out how much storage we need
			FEElement& e = D.ElementRef(0);
			int n = e.Nodes();
			nsize *= n;
		}
	
		a.reserve(nsize);
		if (Save(D, a))
		{
			assert(a.size() == nsize);
			ar.WriteChunk(i+1, a);
		}
	}
}

//-----------------------------------------------------------------------------
void FEFaceData::Save(FEM &fem, Archive& ar)
{
	// loop over all surfaces
	FEMesh& m = fem.m_mesh;
	int NS = m.Surfaces();
	for (int i=0; i<NS; ++i)
	{
		FESurface& S = m.Surface(i);

		// determine data size
		int nsize = VarSize(DataType());
		switch (m_sfmt)
		{
		case FMT_NODE: nsize *= S.Nodes(); break;
		case FMT_ITEM: nsize *= S.Elements(); break;
		case FMT_MULT: nsize *= 4*S.Elements(); break;
		default:
			assert(false);
		}

		// save data
		vector<float> a; a.reserve(nsize);
		if (Save(S, a))
		{
			assert(a.size() == nsize);
			ar.WriteChunk(i+1, a);
		}
	}
}

//=============================================================================
//                            N O D E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Store the nodal displacements
bool FEPlotNodeDisplacement::Save(FEMesh& m, vector<float>& a)
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

		a.push_back(xf[0]);
		a.push_back(xf[1]);
		a.push_back(xf[2]);
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotNodeVelocity::Save(FEMesh& m, vector<float>& a)
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

		a.push_back(xf[0]);
		a.push_back(xf[1]);
		a.push_back(xf[2]);
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotNodeAcceleration::Save(FEMesh& m, vector<float>& a)
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

		a.push_back(xf[0]);
		a.push_back(xf[1]);
		a.push_back(xf[2]);
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidPressure::Save(FEMesh &m, vector<float>& a)
{
	for (int i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);
		a.push_back((float) node.m_pt);
	}
	return true;
}

//=============================================================================
//                           E L E M E N T   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
bool FEPlotElementStress::Save(FEDomain& dom, vector<float>& a)
{
	// write solid stresses
	FESolidDomain* pbd = dynamic_cast<FESolidDomain*>(&dom);
	if (pbd) { WriteSolidStress(*pbd, a); return true; }

	// write shell stresses
	FEShellDomain* pbs = dynamic_cast<FEShellDomain*>(&dom);
	if (pbs) { WriteShellStress(*pbs, a); return true; }

	return false;
}

//-----------------------------------------------------------------------------
void FEPlotElementStress::WriteSolidStress(FESolidDomain& d, vector<float>& a)
{
	// write solid element data
	for (int i=0; i<d.Elements(); ++i)
	{
		FESolidElement& el = d.Element(i);

		float s[6] = {0};
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
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

		a.push_back(s[0]);
		a.push_back(s[1]);
		a.push_back(s[2]);
		a.push_back(s[3]);
		a.push_back(s[4]);
		a.push_back(s[5]);
	}
}

//-----------------------------------------------------------------------------
void FEPlotElementStress::WriteShellStress(FEShellDomain& d, vector<float>& a)
{
	// write shell element data
	for (int i=0; i<d.Elements(); ++i)
	{
		FEShellElement& el = d.Element(i);

		float s[6] = {0};
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
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

		a.push_back(s[0]);
		a.push_back(s[1]);
		a.push_back(s[2]);
		a.push_back(s[3]);
		a.push_back(s[4]);
		a.push_back(s[5]);
	}
}

//-----------------------------------------------------------------------------
bool FEPlotFluidFlux::Save(FEDomain &dom, vector<float>& a)
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

			a.push_back(af[0]);
			a.push_back(af[1]);
			a.push_back(af[2]);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFiberVector::Save(FEDomain &dom, vector<float>& a)
{
	int i, j, n;
	float f[3];
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
			f[0] = (float) r.x;
			f[1] = (float) r.y;
			f[2] = (float) r.z;

			a.push_back(a[0]);
			a.push_back(a[1]);
			a.push_back(a[2]);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Store shell thicknesses
bool FEPlotShellThickness::Save(FEDomain &dom, vector<float> &a)
{
	FEShellDomain* pbs = dynamic_cast<FEShellDomain*>(&dom);
	if (pbs)
	{
		int NS = pbs->Elements();
		for (int i=0; i<NS; ++i)
		{
			FEShellElement& e = pbs->Element(i);
			dom.UnpackElement(e);
			int n = e.Nodes();
			vec3d* D = e.Dt();
			for (int j=0; j<n; ++j)
			{
				double h = e.m_h0[j] * D[j].norm();
				a.push_back((float) h);
			}
		}
		return true;
	}
	return false;
}

//=============================================================================
//                           S U R F A C E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
bool FEPlotContactGap::Save(FESurface& surf, vector<float>& a)
{
	FESlidingSurface* ps = dynamic_cast<FESlidingSurface*>(&surf);
	if (ps)
	{
		int NN = ps->Nodes();
		a.assign(NN, 0.f);
		for (int i=0; i<NN; ++i) a[i] = (float) ps->gap[i];
		return true;
	}

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
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotContactTraction::Save(FESurface &surf, vector<float>& a)
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
	return false;
}
