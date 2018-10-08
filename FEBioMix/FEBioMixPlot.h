#pragma once
#include <FECore/FEPlotData.h>
#include <FECore/FEElement.h>

//=============================================================================
//                         S U R F A C E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Fluid flow rate
//!
class FEPlotMixtureFluidFlowRate : public FEPlotSurfaceData
{
private:
    FEModel*            m_pfem;
    bool                m_binit;
    vector<FEElement*>  m_elem;
    vector<vec3d>       m_area;
    
public:
	FEPlotMixtureFluidFlowRate(FEModel* pfem) : FEPlotSurfaceData(PLT_FLOAT, FMT_REGION){ m_pfem = pfem; m_binit = true; }
    bool Save(FESurface& surf, FEDataStream& a);
};

//=============================================================================
//							D O M A I N   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Actual fluid pressure
class FEPlotActualFluidPressure : public FEPlotDomainData
{
public:
	FEPlotActualFluidPressure(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid flux
class FEPlotFluidFlux : public FEPlotDomainData
{
public:
	FEPlotFluidFlux(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal Fluid flux
class FEPlotNodalFluidFlux : public FEPlotDomainData
{
public:
	FEPlotNodalFluidFlux(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_MULT){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Actual solute concentration
class FEPlotActualSoluteConcentration : public FEPlotDomainData
{
public:
	FEPlotActualSoluteConcentration(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
protected:
	FEModel*	m_pfem;
	vector<int>	m_sol;
};

//-----------------------------------------------------------------------------
//! Solute flux (for biphasic solute problems)
class FEPlotSoluteFlux : public FEPlotDomainData
{
public:
	FEPlotSoluteFlux(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);

protected:
	vector<int>	m_sol;
};

//-----------------------------------------------------------------------------
//! Osmolarity
class FEPlotOsmolarity : public FEPlotDomainData
{
public:
    FEPlotOsmolarity(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
class FEPlotSBMConcentration : public FEPlotDomainData
{
public:
	FEPlotSBMConcentration(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
protected:
	FEModel*	m_pfem;
	vector<int>	m_sbm;
};

//-----------------------------------------------------------------------------
//! Electric potential
class FEPlotElectricPotential : public FEPlotDomainData
{
public:
	FEPlotElectricPotential(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Current density
class FEPlotCurrentDensity : public FEPlotDomainData
{
public:
	FEPlotCurrentDensity(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Referential solid volume fraction
class FEPlotReferentialSolidVolumeFraction : public FEPlotDomainData
{
public:
    FEPlotReferentialSolidVolumeFraction(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fixed charge density
class FEPlotFixedChargeDensity : public FEPlotDomainData
{
public:
	FEPlotFixedChargeDensity(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Referential fixed charge density
class FEPlotReferentialFixedChargeDensity : public FEPlotDomainData
{
public:
	FEPlotReferentialFixedChargeDensity(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal effective fluid pressures
class FEPlotEffectiveFluidPressure : public FEPlotDomainData
{
public:
	FEPlotEffectiveFluidPressure(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_NODE){}
	bool Save(FEDomain& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal effective downstream fluid pressures
class FEPlotEffectiveShellFluidPressure : public FEPlotDomainData
{
public:
	FEPlotEffectiveShellFluidPressure(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_NODE){}
	bool Save(FEDomain& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal effective solute concentrations (for biphasic-solute problems)
class FEPlotEffectiveSoluteConcentration : public FEPlotDomainData
{
public:
	FEPlotEffectiveSoluteConcentration(FEModel* pfem);
	bool Save(FEDomain& m, FEDataStream& a);
protected:
	FEModel*	m_pfem;
	vector<int>	m_sol;
};

//-----------------------------------------------------------------------------
//! Nodal effective solute concentrations (for biphasic-solute problems)
class FEPlotEffectiveShellSoluteConcentration : public FEPlotDomainData
{
public:
	FEPlotEffectiveShellSoluteConcentration(FEModel* pfem);
	bool SetFilter(const char* sz);
	bool SetFilter(int nsol);
	bool Save(FEDomain& m, FEDataStream& a);
protected:
	int			m_nsol;
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
//! Receptor-ligand complex concentration
class FEPlotReceptorLigandConcentration : public FEPlotDomainData
{
public:
	FEPlotReceptorLigandConcentration(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Base class for solid-bound molecule referential apparent density
class FEPlotSBMRefAppDensity : public FEPlotDomainData
{
public:
	FEPlotSBMRefAppDensity(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
protected:
	vector<int>	m_sbm;
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
//! effective elasticity
class FEPlotEffectiveElasticity : public FEPlotDomainData
{
public:
	FEPlotEffectiveElasticity(FEModel* pfem) : FEPlotDomainData(PLT_TENS4FS, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//=============================================================================
//                         S U R F A C E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Fluid force
//!
class FEPlotFluidForce : public FEPlotSurfaceData
{
public:
	FEPlotFluidForce(FEModel* pfem) : FEPlotSurfaceData(PLT_VEC3F, FMT_REGION){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid force
//! 
class FEPlotFluidForce2 : public FEPlotSurfaceData
{
public:
	FEPlotFluidForce2(FEModel* pfem) : FEPlotSurfaceData(PLT_VEC3F, FMT_REGION){}
	bool Save(FESurface& surf, FEDataStream& a);
};


//-----------------------------------------------------------------------------
//! Fluid pressure gap
//!
class FEPlotPressureGap : public FEPlotSurfaceData
{
public:
	FEPlotPressureGap(FEModel* pfem) : FEPlotSurfaceData(PLT_FLOAT, FMT_MULT){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid load support
//!
class FEPlotFluidLoadSupport : public FEPlotSurfaceData
{
public:
    FEPlotFluidLoadSupport(FEModel* pfem) : FEPlotSurfaceData(PLT_FLOAT, FMT_REGION){}
    bool Save(FESurface& surf, FEDataStream& a);
};

