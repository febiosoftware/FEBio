#pragma once
#include "FECore/FEPlotData.h"
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
	FEPlotElementStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
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
		FEPlotRelativeVolume(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
		bool Save(FEDomain& dom, vector<float>& a);
	};

//-----------------------------------------------------------------------------
//! Actual fluid pressure
class FEPlotActualFluidPressure : public FEDomainData
	{
	public:
		FEPlotActualFluidPressure(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
		bool Save(FEDomain& dom, vector<float>& a);
	};

//-----------------------------------------------------------------------------
//! Fluid flux
class FEPlotFluidFlux : public FEDomainData
{
public:
	FEPlotFluidFlux(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM){}
	bool Save(FEDomain& dom, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Actual solute concentration
class FEPlotActualSoluteConcentration : public FEDomainData
	{
	public:
		FEPlotActualSoluteConcentration(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
		bool Save(FEDomain& dom, vector<float>& a);
	};

//-----------------------------------------------------------------------------
//! Solute flux
class FEPlotSoluteFlux : public FEDomainData
	{
	public:
		FEPlotSoluteFlux(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM){}
		bool Save(FEDomain& dom, vector<float>& a);
	};

//-----------------------------------------------------------------------------
//! Material fibers
class FEPlotFiberVector : public FEDomainData
{
public:
	FEPlotFiberVector(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM){}
	bool Save(FEDomain& dom, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Shell thicknesses
class FEPlotShellThickness : public FEDomainData
{
public:
	FEPlotShellThickness(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_MULT){}
	bool Save(FEDomain& dom, vector<float>& a);
};


//-----------------------------------------------------------------------------
//! Nodal effective fluid pressures
class FEPlotEffectiveFluidPressure : public FEDomainData
{
public:
	FEPlotEffectiveFluidPressure(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_NODE){}
	bool Save(FEDomain& m, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Nodal effective solute concentrations
class FEPlotEffectiveSoluteConcentration : public FEDomainData
{
public:
	FEPlotEffectiveSoluteConcentration(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_NODE){}
	bool Save(FEDomain& m, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Damage reduction factor
class FEPlotDamage : public FEDomainData
{
public:
	FEPlotDamage(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& m, vector<float>& a);
};
