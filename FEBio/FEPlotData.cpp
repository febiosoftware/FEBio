#include "stdafx.h"
#include "FEPlotData.h"
#include "fem.h"
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
REGISTER_PLOTDATA(FEPlotContactPressure , "contact pressure");
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
void FEDomainData::Save(FEM &fem, Archive& ar)
{
	// loop over all domains
	FEMesh& m = fem.m_mesh;
	int ND = m.Domains();
	for (int i=0; i<ND; ++i)
	{
		// get the domain
		FEDomain& D = m.Domain(i);

		// calculate the size of the data vector
		int nsize = VarSize(DataType());
		switch (m_sfmt)
		{
		case FMT_NODE: nsize *= D.Nodes(); break;
		case FMT_ITEM: nsize *= D.Elements(); break;
		case FMT_MULT:
			{
				// since all elements have the same type within a domain
				// we just grab the number of nodes of the first element 
				// to figure out how much storage we need
				FEElement& e = D.ElementRef(0);
				int n = e.Nodes();
				nsize *= n*D.Elements();
			}
			break;
		default:
			assert(false);
		}

		// fill data vector and save
		vector<float> a; 
		a.reserve(nsize);
		if (Save(D, a))
		{
			assert(a.size() == nsize);
			ar.WriteChunk(i+1, a);
		}
	}
}

//-----------------------------------------------------------------------------
void FESurfaceData::Save(FEM &fem, Archive& ar)
{
	// loop over all surfaces
	FEMesh& m = fem.m_mesh;
	int NS = m.Surfaces();
	for (int i=0; i<NS; ++i)
	{
		FESurface& S = m.Surface(i);

		// Determine data size.
		// Note that for the FMT_MULT case we are 
		// assuming four data entries per facet
		// regardless of the nr of nodes a facet really has
		// this is because for surfaces, all elements are not
		// necessarily of the same type
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

//=============================================================================
//                           E L E M E N T   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
bool FEPlotElementStress::Save(FEDomain& dom, vector<float>& a)
{
	// write solid stresses
	FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&dom);
	if (pbd) return WriteSolidStress(*pbd, a);

	// write shell stresses
	FEElasticShellDomain* pbs = dynamic_cast<FEElasticShellDomain*>(&dom);
	if (pbs) return WriteShellStress(*pbs, a);

	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotElementStress::WriteSolidStress(FEElasticSolidDomain& d, vector<float>& a)
{
	// make sure this is not a rigid body
	if (dynamic_cast<FERigidSolidDomain*>(&d)) return false;

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

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotElementStress::WriteShellStress(FEElasticShellDomain& d, vector<float>& a)
{
	// make sure this is not a rigid body
	if (dynamic_cast<FERigidShellDomain*>(&d)) return false;

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

	return true;
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
	FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&dom);
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

			a.push_back(f[0]);
			a.push_back(f[1]);
			a.push_back(f[2]);
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

//-----------------------------------------------------------------------------
bool FEPlotFluidPressure::Save(FEDomain &dom, vector<float>& a)
{
	FEPoroSolidDomain* pd = dynamic_cast<FEPoroSolidDomain*>(&dom);
	if (pd)
	{
		int N = pd->Nodes();
		for (int i=0; i<N; ++i)
		{
			FENode& node = pd->Node(i);
			a.push_back((float) node.m_pt);
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
	if (ps) return SaveSliding(*ps, a);

	FEFacetSlidingSurface* pf = dynamic_cast<FEFacetSlidingSurface*>(&surf);
	if (pf) return SaveFacetSliding(*pf, a);

	FESlidingSurface2* ps2 = dynamic_cast<FESlidingSurface2*>(&surf);
	if (ps2) return SaveSliding2(*ps2, a);

	FETiedContactSurface* pt = dynamic_cast<FETiedContactSurface*>(&surf);
	if (pt) return SaveTied(*pt, a);

	return false;
}

//-----------------------------------------------------------------------------
// Store the gap values for a sliding surface. Although the gap values are
// stored at the nodes, we use the FMT_MULT format to be consistent with
// the other contact surfaces.
bool FEPlotContactGap::SaveSliding(FESlidingSurface& s, vector<float>& a)
{
	int NF = s.Elements();
	a.assign(4*NF, 0.f);
	for (int i=0; i<NF; ++i) 
	{
		FESurfaceElement& f = s.Element(i);
		a[4*i  ] = (float) s.gap[f.m_lnode[0]];
		a[4*i+1] = (float) s.gap[f.m_lnode[1]];
		a[4*i+2] = (float) s.gap[f.m_lnode[2]];
		a[4*i+3] = (float) s.gap[f.m_lnode[3]];
	}
	return true;
}

//-----------------------------------------------------------------------------
// Store the gap values for a tied surface. Although the gap values are
// stored at the nodes, we use the FMT_MULT format to be consistent with
// the other contact surfaces.
bool FEPlotContactGap::SaveTied(FETiedContactSurface& s, vector<float>& a)
{
	int NF = s.Elements();
	a.assign(4*NF, 0.f);
	for (int i=0; i<NF; ++i) 
	{
		FESurfaceElement& f = s.Element(i);
		a[4*i  ] = (float) s.gap[f.m_lnode[0]].norm();
		a[4*i+1] = (float) s.gap[f.m_lnode[1]].norm();
		a[4*i+2] = (float) s.gap[f.m_lnode[2]].norm();
		a[4*i+3] = (float) s.gap[f.m_lnode[3]].norm();
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotContactGap::SaveFacetSliding(FEFacetSlidingSurface& s, vector<float>& a)
{
	int NF = s.Elements();
	a.assign(4*NF, 0.f);
	int nint = 0;
	double gi[4], gn[4];
	int i, k;
	for (i=0; i<NF; ++i)
	{
		FESurfaceElement& el = s.Element(i);
		int ne = el.Nodes();
		int ni = el.GaussPoints();
		for (k=0; k<ni; ++k, ++nint) gi[k] = s.m_gap[nint];

		for (k=0; k<ni; ++k) if (gi[k] < 0) gi[k] = 0;
		el.project_to_nodes(gi, gn);

		for (k=0; k<ni; ++k) if (gn[k] < 0) gn[k] = 0;

		a[4*i  ] = (float) gn[0];
		a[4*i+1] = (float) gn[1];
		a[4*i+2] = (float) gn[2];
		a[4*i+3] = (float) (ne == 4? gn[3] : gn[2]);
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotContactGap::SaveSliding2(FESlidingSurface2& s, vector<float>& a)
{
	int NF = s.Elements();
	a.assign(4*NF, 0.f);
	int nint = 0;
	double gi[4], gn[4];
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& el = s.Element(i);
		int ne = el.Nodes();
		int ni = el.GaussPoints();
		for (int k=0; k<ni; ++k, ++nint) gi[k] = s.m_gap[nint];

		el.project_to_nodes(gi, gn);

		a[4*i  ] = (float) gn[0];
		a[4*i+1] = (float) gn[1];
		a[4*i+2] = (float) gn[2];
		a[4*i+3] = (float) (ne == 4? gn[3] : gn[2]);
	}
	return true;
}

//-----------------------------------------------------------------------------
// Store the contact pressures
bool FEPlotContactPressure::Save(FESurface &surf, vector<float>& a)
{
	FESlidingSurface* ps = dynamic_cast<FESlidingSurface*>(&surf);
	if (ps) return SaveSliding(*ps, a);

	FEFacetSlidingSurface* pf = dynamic_cast<FEFacetSlidingSurface*>(&surf);
	if (pf) return SaveFacetSliding(*pf, a);

	FESlidingSurface2* ps2 = dynamic_cast<FESlidingSurface2*>(&surf);
	if (ps2) return SaveSliding2(*ps2, a);

	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotContactPressure::SaveSliding(FESlidingSurface &s, vector<float>& a)
{
	int NF = s.Elements();
	a.assign(4*NF, 0.f);
	for (int i=0; i<NF; ++i) 
	{
		FESurfaceElement& f = s.Element(i);
		for (int j=0; j<4; ++j)
		{
			double g = s.gap[f.m_lnode[j]];
			double L = s.Lm[f.m_lnode[j]];
			double e = s.eps[f.m_lnode[j]];
			double t = MBRACKET(L + g*e);
			a[4*i+j] = (float) t;
		}
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotContactPressure::SaveFacetSliding(FEFacetSlidingSurface &s, vector<float>& a)
{
	int NF = s.Elements();
	a.assign(4*NF, 0.f);
	int nint = 0;
	double ti[4], tn[4];
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& el = s.Element(i);
		int ne = el.Nodes();
		int ni = el.GaussPoints();
		for (int k=0; k<ni; ++k, ++nint)
		{
			double L = s.m_Lm[nint];
			ti[k] = L;// + pf->m_epsn*gi[k];
			ti[k] = (ti[k]>=0?ti[k] : 0);		
		}

		el.project_to_nodes(ti, tn);

		a[4*i  ] = (float) tn[0];
		a[4*i+1] = (float) tn[1];
		a[4*i+2] = (float) tn[2];
		a[4*i+3] = (float) (ne == 4? tn[3] : tn[2]);
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotContactPressure::SaveSliding2(FESlidingSurface2 &s, vector<float>& a)
{
	int NF = s.Elements();
	a.assign(4*NF, 0.f);
	int nint = 0;
	double ti[4], tn[4];
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& el = s.Element(i);
		int ne = el.Nodes();
		int ni = el.GaussPoints();
		for (int k=0; k<ni; ++k, ++nint)
		{
			double L = s.m_Lmd[nint];
			ti[k] = L;// + pf->m_epsn*gi[k];
			ti[k] = (ti[k]>=0?ti[k] : 0);		
		}

		el.project_to_nodes(ti, tn);

		a[4*i  ] = (float) tn[0];
		a[4*i+1] = (float) tn[1];
		a[4*i+2] = (float) tn[2];
		a[4*i+3] = (float) (ne == 4? tn[3] : tn[2]);
	}
	return true;
}
