#include "stdafx.h"
#include "FEPlotSurfaceData.h"
#include "FEBioMech/FEElasticSolidDomain.h"
#include "FECore/febio.h"
#include "FEBioPlotFile.h"

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
	int MFN = FEBioPlotFile::PLT_MAX_FACET_NODES;
	a.assign(MFN*NF, 0.f);
	for (int i=0; i<NF; ++i) 
	{
		FESurfaceElement& f = s.Element(i);
		int ne = f.m_lnode.size();
		for (int j= 0; j< ne; ++j) a[MFN*i + j] = (float) s.m_gap[f.m_lnode[j]];
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
	int MFN = FEBioPlotFile::PLT_MAX_FACET_NODES;
	a.assign(MFN*NF, 0.f);
	for (int i=0; i<NF; ++i) 
	{
		FESurfaceElement& f = s.Element(i);
		int ne = f.m_lnode.size();
		for (int j= 0; j< ne; ++j) a[MFN*i + j] = (float) s.m_gap[f.m_lnode[j]].norm();
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotContactGap::SaveFacetSliding(FEFacetSlidingSurface& s, vector<float>& a)
{
	int NF = s.Elements();
	int NFM = FEBioPlotFile::PLT_MAX_FACET_NODES;
	a.assign(NFM*NF, 0.f);
	int nint = 0;
	double gi[FEElement::MAX_NODES], gn[FEElement::MAX_NODES];
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& el = s.Element(i);
		int ne = el.Nodes();
		int ni = el.GaussPoints();
		for (int k=0; k<ni; ++k, ++nint) gi[k] = s.m_Data[i][k].m_gap;

		for (int k=0; k<ni; ++k) if (gi[k] < 0) gi[k] = 0;
		el.project_to_nodes(gi, gn);

		for (int k=0; k<ne; ++k) if (gn[k] < 0) gn[k] = 0;

		for (int k=0; k<ne; ++k) a[NFM*i + k] = (float) gn[k];
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotContactGap::SaveSliding2(FESlidingSurface2& s, vector<float>& a)
{
	int NF = s.Elements();
	const int NFM = FEBioPlotFile::PLT_MAX_FACET_NODES;
	a.assign(NFM*NF, 0.f);
	double gi[NFM], gn[NFM];
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& el = s.Element(i);
		int ne = el.Nodes();
		int ni = el.GaussPoints();
		for (int k=0; k<ni; ++k) gi[k] = s.m_Data[i][k].m_gap;

		el.project_to_nodes(gi, gn);

		for (int k=0; k<ne; ++k) a[NFM*i + k] = (float) gn[k];
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotContactGap::SaveSliding3(FESlidingSurface3& s, vector<float>& a)
{
	int NF = s.Elements();
	const int NFM = FEBioPlotFile::PLT_MAX_FACET_NODES;
	double gi[NFM], gn[NFM];
	a.assign(NFM*NF, 0.f);
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& el = s.Element(i);
		int ne = el.Nodes();
		int ni = el.GaussPoints();
		for (int k=0; k<ni; ++k) gi[k] = s.m_Data[i][k].m_gap;
		
		el.project_to_nodes(gi, gn);
		
		for (int k=0; k<ne; ++k) a[NFM*i + k] = (float) gn[k];
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
	const int MFN = FEBioPlotFile::PLT_MAX_FACET_NODES;
	a.assign(MFN*NF, 0.f);
	for (int i=0; i<NF; ++i) 
	{
		FESurfaceElement& f = s.Element(i);
		int ne = f.Nodes();
		for (int j=0; j<ne; ++j) a[MFN*i+j] = (float) s.m_Ln[f.m_lnode[j]];
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotContactPressure::SaveFacetSliding(FEFacetSlidingSurface &s, vector<float>& a)
{
	int NF = s.Elements();
	const int MFN = FEBioPlotFile::PLT_MAX_FACET_NODES;
	a.assign(MFN*NF, 0.f);
	double ti[MFN], tn[MFN];
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& el = s.Element(i);
		int ne = el.Nodes();
		int ni = el.GaussPoints();
		for (int k=0; k<ni; ++k)
		{
			double L = s.m_Data[i][k].m_Ln;
			ti[k] = L;// + pf->m_epsn*gi[k];
			ti[k] = (ti[k]>=0?ti[k] : 0);		
		}

		el.project_to_nodes(ti, tn);
		for (int k=0; k<ni; ++k)
			tn[k] = (tn[k]>=0?tn[k] : 0);		

		for (int k=0; k<ne; ++k) a[MFN*i + k] = (float) tn[k];
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotContactPressure::SaveSliding2(FESlidingSurface2 &s, vector<float>& a)
{
	int NF = s.Elements();
	const int MFN = FEBioPlotFile::PLT_MAX_FACET_NODES;
	a.assign(MFN*NF, 0.f);
	double ti[MFN], tn[MFN];
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& el = s.Element(i);
		int ne = el.Nodes();
		int ni = el.GaussPoints();
		for (int k=0; k<ni; ++k) ti[k] = s.m_Data[i][k].m_Ln;

		el.project_to_nodes(ti, tn);

		for (int k=0; k<ne; ++k) a[MFN*i + k] = (float) tn[k];
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotContactPressure::SaveSliding3(FESlidingSurface3 &s, vector<float>& a)
{
	int NF = s.Elements();
	const int MFN = FEBioPlotFile::PLT_MAX_FACET_NODES;
	a.assign(MFN*NF, 0.f);
	double ti[MFN], tn[MFN];
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& el = s.Element(i);
		int ne = el.Nodes();
		int ni = el.GaussPoints();
		for (int k=0; k<ni; ++k) ti[k] = s.m_Data[i][k].m_Ln;
		
		el.project_to_nodes(ti, tn);
		
		for (int k=0; k<ne; ++k) a[MFN*i + k] = (float) tn[k];
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
	const int MFN = FEBioPlotFile::PLT_MAX_FACET_NODES;
	vec3d t[MFN];
	for (int i=0; i<s.Elements(); ++i)
	{
		FESurfaceElement& e = s.Element(i);
		int ne = e.Nodes();

		for (int j=0; j<ne; ++j)
		{
			int nj = e.m_lnode[j];
			double gi = s.m_gap[nj];
			double Li = s.m_Ln[nj];
			vec3d ti = s.m_nu[nj];
			if (gi > 0) t[j] = ti*Li; else t[j] = vec3d(0,0,0);
		}

		for (int j=0; j<MFN; ++j)
		{
			a.push_back((float) t[j].x);
			a.push_back((float) t[j].y);
			a.push_back((float) t[j].z); 
		}
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotContactTraction::SaveFacetSliding(FEFacetSlidingSurface &s, std::vector<float> &a)
{
	const int MFN = FEBioPlotFile::PLT_MAX_FACET_NODES;
	a.assign(3*MFN*s.Elements(), 0.f);
	double tix[MFN], tiy[MFN], tiz[MFN];
	double tnx[MFN], tny[MFN], tnz[MFN];
	for (int j=0; j<s.Elements(); ++j)
	{
		FESurfaceElement& el = s.Element(j);
		int ne = el.Nodes();
		int ni = el.GaussPoints();

		vec3d t;
		for (int k=0; k<ni; ++k)
		{
			FEFacetSlidingSurface::Data& pt = s.m_Data[j][k];
			double gi = pt.m_gap;
			double Li = pt.m_Ln;
			vec3d  ti = pt.m_nu;
			if (gi > 0) t = ti*(Li); else t = vec3d(0,0,0);
			tix[k] = t.x; tiy[k] = t.y; tiz[k] = t.z;
		}

		// project traction to nodes
		el.project_to_nodes(tix, tnx);
		el.project_to_nodes(tiy, tny);
		el.project_to_nodes(tiz, tnz);

		// store in archive
		for (int k=0; k<ne; ++k)
		{
			a[3*MFN*j +3*k   ] = (float) tnx[k];
			a[3*MFN*j +3*k +1] = (float) tny[k];
			a[3*MFN*j +3*k +2] = (float) tnz[k];
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotContactTraction::SaveSliding2(FESlidingSurface2 &s, std::vector<float> &a)
{
	const int MFN = FEBioPlotFile::PLT_MAX_FACET_NODES;
	a.assign(3*MFN*s.Elements(), 0.f);
	double tix[MFN], tiy[MFN], tiz[MFN];
	double tnx[MFN], tny[MFN], tnz[MFN];
	for (int j=0; j<s.Elements(); ++j)
	{
		FESurfaceElement& el = s.Element(j);
		int ne = el.Nodes();
		int ni = el.GaussPoints();

		vec3d t;
		for (int k=0; k<ni; ++k)
		{
			FESlidingSurface2::Data& pt = s.m_Data[j][k];
			double gi = pt.m_gap;
			double Li = pt.m_Ln;
			vec3d  ti = pt.m_nu;
			if (gi > 0) t = ti*(Li); else t = vec3d(0,0,0);
			tix[k] = t.x; tiy[k] = t.y; tiz[k] = t.z;
		}

		// project traction to nodes
		el.project_to_nodes(tix, tnx);
		el.project_to_nodes(tiy, tny);
		el.project_to_nodes(tiz, tnz);

		// store in archive
		for (int k=0; k<ne; ++k)
		{
			a[3*MFN*j +3*k   ] = (float) tnx[k];
			a[3*MFN*j +3*k +1] = (float) tny[k];
			a[3*MFN*j +3*k +2] = (float) tnz[k];
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotContactTraction::SaveSliding3(FESlidingSurface3 &s, std::vector<float> &a)
{
	const int MFN = FEBioPlotFile::PLT_MAX_FACET_NODES;
	a.assign(3*MFN*s.Elements(), 0.f);
	double tix[MFN], tiy[MFN], tiz[MFN];
	double tnx[MFN], tny[MFN], tnz[MFN];
	for (int j=0; j<s.Elements(); ++j)
	{
		FESurfaceElement& el = s.Element(j);
		int ne = el.Nodes();
		int ni = el.GaussPoints();

		vec3d t;
		for (int k=0; k<ni; ++k)
		{
			FESlidingSurface3::Data& pt = s.m_Data[j][k];
			double gi = pt.m_gap;
			double Li = pt.m_Ln;
			vec3d  ti = pt.m_nu;
			if (gi > 0) t = ti*(Li); else t = vec3d(0,0,0);
			tix[k] = t.x; tiy[k] = t.y; tiz[k] = t.z;
		}

		// project traction to nodes
		el.project_to_nodes(tix, tnx);
		el.project_to_nodes(tiy, tny);
		el.project_to_nodes(tiz, tnz);

		// store in archive
		for (int k=0; k<ne; ++k)
		{
			a[3*MFN*j +3*k   ] = (float) tnx[k];
			a[3*MFN*j +3*k +1] = (float) tny[k];
			a[3*MFN*j +3*k +2] = (float) tnz[k];
		}
	}

	return true;
}
