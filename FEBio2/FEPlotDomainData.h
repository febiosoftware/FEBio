#pragma once
#include "FEPlotData.h"
#include "FEElasticShellDomain.h"
#include "FEElasticSolidDomain.h"

//=============================================================================
// The following classes define data fields for domains. 
//=============================================================================

//-----------------------------------------------------------------------------
//! Element stresses
class FEPlotElementStress : public FEDomainData
{
public:
	FEPlotElementStress() : FEDomainData(MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& dom, vector<float>& a);

protected:
	bool WriteSolidStress(FEElasticSolidDomain& d, vector<float>& a);
	bool WriteShellStress(FEElasticShellDomain& d, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Relative volume
class FEPlotRelativeVolume : public FEDomainData
	{
	public:
		FEPlotRelativeVolume() : FEDomainData(FLOAT, FMT_ITEM){}
		bool Save(FEDomain& dom, vector<float>& a);
	};

//-----------------------------------------------------------------------------
//! Actual fluid pressure
class FEPlotActualFluidPressure : public FEDomainData
	{
	public:
		FEPlotActualFluidPressure() : FEDomainData(FLOAT, FMT_ITEM){}
		bool Save(FEDomain& dom, vector<float>& a);
	};

//-----------------------------------------------------------------------------
//! Fluid flux
class FEPlotFluidFlux : public FEDomainData
{
public:
	FEPlotFluidFlux() : FEDomainData(VEC3F, FMT_ITEM){}
	bool Save(FEDomain& dom, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Actual solute concentration
class FEPlotActualSoluteConcentration : public FEDomainData
	{
	public:
		FEPlotActualSoluteConcentration() : FEDomainData(FLOAT, FMT_ITEM){}
		bool Save(FEDomain& dom, vector<float>& a);
	};

//-----------------------------------------------------------------------------
//! Solute flux
class FEPlotSoluteFlux : public FEDomainData
	{
	public:
		FEPlotSoluteFlux() : FEDomainData(VEC3F, FMT_ITEM){}
		bool Save(FEDomain& dom, vector<float>& a);
	};

//-----------------------------------------------------------------------------
//! Material fibers
class FEPlotFiberVector : public FEDomainData
{
public:
	FEPlotFiberVector() : FEDomainData(VEC3F, FMT_ITEM){}
	bool Save(FEDomain& dom, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Shell thicknesses
class FEPlotShellThickness : public FEDomainData
{
public:
	FEPlotShellThickness() : FEDomainData(FLOAT, FMT_MULT){}
	bool Save(FEDomain& dom, vector<float>& a);
};


//-----------------------------------------------------------------------------
//! Nodal effective fluid pressures
class FEPlotEffectiveFluidPressure : public FEDomainData
{
public:
	FEPlotEffectiveFluidPressure() : FEDomainData(FLOAT, FMT_NODE){}
	bool Save(FEDomain& m, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Nodal effective solute concentrations
class FEPlotEffectiveSoluteConcentration : public FEDomainData
{
public:
	FEPlotEffectiveSoluteConcentration() : FEDomainData(FLOAT, FMT_NODE){}
	bool Save(FEDomain& m, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Damage reduction factor
class FEPlotDamage : public FEDomainData
{
public:
	FEPlotDamage() : FEDomainData(FLOAT, FMT_ITEM){}
	bool Save(FEDomain& m, vector<float>& a);
};
