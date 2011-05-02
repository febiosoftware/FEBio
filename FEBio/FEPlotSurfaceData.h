#pragma once
#include "FEPlotData.h"

//=============================================================================
//                         S U R F A C E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! (this is just test data: don't use it)
class FETestData : public FEDomainData
{
public:
	FETestData() : FEDomainData(FLOAT, FMT_NODE){}
	bool Save(FEDomain& dom, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Contact gap
//!
class FEPlotContactGap : public FESurfaceData
{
public:
	FEPlotContactGap() : FESurfaceData(FLOAT, FMT_MULT){}
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
	FEPlotContactPressure() : FESurfaceData(FLOAT, FMT_MULT){}
	bool Save(FESurface& surf, vector<float>& a);

protected:
	bool SaveSliding     (FESlidingSurface&      s, vector<float>& a);
	bool SaveFacetSliding(FEFacetSlidingSurface& s, vector<float>& a);
	bool SaveSliding2    (FESlidingSurface2&	 s, vector<float>& a);
	bool SaveSliding3    (FESlidingSurface3&	 s, vector<float>& a);
};
