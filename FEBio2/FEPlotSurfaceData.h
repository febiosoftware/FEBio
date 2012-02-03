#pragma once
#include "FECore/FEPlotData.h"
#include "FEBioLib/FESlidingInterface.h"
#include "FEBioLib/FEFacet2FacetSliding.h"
#include "FESlidingInterface2.h"
#include "FESlidingInterface3.h"
#include "FEBioLib/FETiedInterface.h"

//=============================================================================
//                         S U R F A C E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Contact gap
//!
class FEPlotContactGap : public FESurfaceData
{
public:
	FEPlotContactGap(FEModel* pfem) : FESurfaceData(PLT_FLOAT, FMT_MULT){}
	bool Save(FESurface& surf, vector<float>& a);

protected:
	bool SaveSliding     (FESlidingSurface&      s, vector<float>& a);
	bool SaveFacetSliding(FEFacetSlidingSurface& s, vector<float>& a);
	bool SaveSliding2    (FESlidingSurface2&	 s, vector<float>& a);
	bool SaveSliding3    (FESlidingSurface3&	 s, vector<float>& a);
	bool SaveTied        (FETiedContactSurface&	 s, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Contact pressure
//!
class FEPlotContactPressure : public FESurfaceData
{
public:
	FEPlotContactPressure(FEModel* pfem) : FESurfaceData(PLT_FLOAT, FMT_MULT){}
	bool Save(FESurface& surf, vector<float>& a);

protected:
	bool SaveSliding     (FESlidingSurface&      s, vector<float>& a);
	bool SaveFacetSliding(FEFacetSlidingSurface& s, vector<float>& a);
	bool SaveSliding2    (FESlidingSurface2&	 s, vector<float>& a);
	bool SaveSliding3    (FESlidingSurface3&	 s, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Contact traction
//!
class FEPlotContactTraction : public FESurfaceData
{
public:
	FEPlotContactTraction(FEModel* pfem) : FESurfaceData(PLT_VEC3F, FMT_MULT){}
	bool Save(FESurface& surf, vector<float>& a);

protected:
	bool SaveSliding(FESlidingSurface& s, vector<float>& a);
	bool SaveFacetSliding(FEFacetSlidingSurface& s, vector<float>& a);
	bool SaveSliding2    (FESlidingSurface2&	 s, vector<float>& a);
	bool SaveSliding3    (FESlidingSurface3&	 s, vector<float>& a);
};
