#include "stdafx.h"
#include "FEPlotSurfaceData.h"
#include "FEElasticSolidDomain.h"
#include "fem.h"

//-----------------------------------------------------------------------------
REGISTER_FEBIO_CLASS(FEPlotContactGap      , FEPlotData, "contact gap"     );
REGISTER_FEBIO_CLASS(FEPlotContactPressure , FEPlotData, "contact pressure");
REGISTER_FEBIO_CLASS(FEPlotContactTraction , FEPlotData, "contact traction");

//=============================================================================
// Contact Gap
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

	FESlidingSurface3* ps3 = dynamic_cast<FESlidingSurface3*>(&surf);
	if (ps3) return SaveSliding3(*ps3, a);
	
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
bool FEPlotContactGap::SaveSliding3(FESlidingSurface3& s, vector<float>& a)
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

//=============================================================================
// Contact Pressure
//=============================================================================
//
// Store the contact pressures
bool FEPlotContactPressure::Save(FESurface &surf, vector<float>& a)
{
	FESlidingSurface* ps = dynamic_cast<FESlidingSurface*>(&surf);
	if (ps) return SaveSliding(*ps, a);

	FEFacetSlidingSurface* pf = dynamic_cast<FEFacetSlidingSurface*>(&surf);
	if (pf) return SaveFacetSliding(*pf, a);

	FESlidingSurface2* ps2 = dynamic_cast<FESlidingSurface2*>(&surf);
	if (ps2) return SaveSliding2(*ps2, a);

	FESlidingSurface3* ps3 = dynamic_cast<FESlidingSurface3*>(&surf);
	if (ps3) return SaveSliding3(*ps3, a);
	
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
			double L = s.m_Ln[nint];
			ti[k] = L;// + pf->m_epsn*gi[k];
			ti[k] = (ti[k]>=0?ti[k] : 0);		
		}

		el.project_to_nodes(ti, tn);
		for (int k=0; k<ni; ++k)
			tn[k] = (tn[k]>=0?tn[k] : 0);		

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
		for (int k=0; k<ni; ++k, ++nint) ti[k] = s.m_Ln[nint];

		el.project_to_nodes(ti, tn);

		a[4*i  ] = (float) tn[0];
		a[4*i+1] = (float) tn[1];
		a[4*i+2] = (float) tn[2];
		a[4*i+3] = (float) (ne == 4? tn[3] : tn[2]);
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotContactPressure::SaveSliding3(FESlidingSurface3 &s, vector<float>& a)
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
		for (int k=0; k<ni; ++k, ++nint) ti[k] = s.m_Ln[nint];
		
		el.project_to_nodes(ti, tn);
		
		a[4*i  ] = (float) tn[0];
		a[4*i+1] = (float) tn[1];
		a[4*i+2] = (float) tn[2];
		a[4*i+3] = (float) (ne == 4? tn[3] : tn[2]);
	}
	return true;
}

//=============================================================================
// Contact traction
//=============================================================================

bool FEPlotContactTraction::Save(FESurface &surf, std::vector<float> &a)
{
	FESlidingSurface* ps = dynamic_cast<FESlidingSurface*>(&surf);
	if (ps) return SaveSliding(*ps, a);

	FEFacetSlidingSurface* pf = dynamic_cast<FEFacetSlidingSurface*>(&surf);
	if (pf) return SaveFacetSliding(*pf, a);

	FESlidingSurface2* ps2 = dynamic_cast<FESlidingSurface2*>(&surf);
	if (ps2) return SaveSliding2(*ps2, a);

	FESlidingSurface3* ps3 = dynamic_cast<FESlidingSurface3*>(&surf);
	if (ps3) return SaveSliding3(*ps3, a);
	
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotContactTraction::SaveSliding(FESlidingSurface &s, std::vector<float> &a)
{
	vec3d t[4];
	for (int i=0; i<s.Elements(); ++i)
	{
		FESurfaceElement& e = s.Element(i);
		int ne = e.Nodes();

		for (int j=0; j<ne; ++j)
		{
			int nj = e.m_lnode[j];
			double gi = s.gap[nj];
			double Li = s.m_Ln[nj];
			vec3d ti = s.nu[nj];
			if (gi > 0) t[j] = ti*Li; else t[j] = vec3d(0,0,0);
		}
		if (ne == 3) t[3] = t[2];

		a.push_back((float) t[0].x); a.push_back((float) t[0].y); a.push_back((float) t[0].z);
		a.push_back((float) t[1].x); a.push_back((float) t[1].y); a.push_back((float) t[1].z);
		a.push_back((float) t[2].x); a.push_back((float) t[2].y); a.push_back((float) t[2].z);
		a.push_back((float) t[3].x); a.push_back((float) t[3].y); a.push_back((float) t[3].z);
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotContactTraction::SaveFacetSliding(FEFacetSlidingSurface &s, std::vector<float> &a)
{
	int nint = 0;
	for (int j=0; j<s.Elements(); ++j)
	{
		FESurfaceElement& el = s.Element(j);
		int ne = el.Nodes();
		int ni = el.GaussPoints();

		// calculate the average traction
		vec3d t(0,0,0);
		for (int k=0; k<ni; ++k, ++nint)
		{
			double gi = s.m_gap[nint];
			double Li = s.m_Ln[nint];
			vec3d  ti = s.m_nu[nint];
			if (gi > 0) t += ti*(Li);
		}
		t /= (double) ni;

		// project average traction to nodes
		a.push_back((float) t.x); a.push_back((float) t.y); a.push_back((float) t.z);
		a.push_back((float) t.x); a.push_back((float) t.y); a.push_back((float) t.z);
		a.push_back((float) t.x); a.push_back((float) t.y); a.push_back((float) t.z);
		a.push_back((float) t.x); a.push_back((float) t.y); a.push_back((float) t.z);
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotContactTraction::SaveSliding2(FESlidingSurface2 &s, std::vector<float> &a)
{
	int nint = 0;
	for (int j=0; j<s.Elements(); ++j)
	{
		FESurfaceElement& el = s.Element(j);
		int ne = el.Nodes();
		int ni = el.GaussPoints();

		// calculate the average traction
		vec3d t(0,0,0);
		for (int k=0; k<ni; ++k, ++nint)
		{
			double gi = s.m_gap[nint];
			double Li = s.m_Ln[nint];
			vec3d  ti = s.m_nu[nint];
			if (gi > 0) t += ti*(Li);
		}
		t /= (double) ni;

		// project average traction to nodes
		a.push_back((float) t.x); a.push_back((float) t.y); a.push_back((float) t.z);
		a.push_back((float) t.x); a.push_back((float) t.y); a.push_back((float) t.z);
		a.push_back((float) t.x); a.push_back((float) t.y); a.push_back((float) t.z);
		a.push_back((float) t.x); a.push_back((float) t.y); a.push_back((float) t.z);
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotContactTraction::SaveSliding3(FESlidingSurface3 &s, std::vector<float> &a)
{
	int nint = 0;
	for (int j=0; j<s.Elements(); ++j)
	{
		FESurfaceElement& el = s.Element(j);
		int ne = el.Nodes();
		int ni = el.GaussPoints();

		// calculate the average traction
		vec3d t(0,0,0);
		for (int k=0; k<ni; ++k, ++nint)
		{
			double gi = s.m_gap[nint];
			double Li = s.m_Ln[nint];
			vec3d  ti = s.m_nu[nint];
			if (gi > 0) t += ti*(Li);
		}
		t /= (double) ni;

		// project average traction to nodes
		a.push_back((float) t.x); a.push_back((float) t.y); a.push_back((float) t.z);
		a.push_back((float) t.x); a.push_back((float) t.y); a.push_back((float) t.z);
		a.push_back((float) t.x); a.push_back((float) t.y); a.push_back((float) t.z);
		a.push_back((float) t.x); a.push_back((float) t.y); a.push_back((float) t.z);
	}

	return true;
}
