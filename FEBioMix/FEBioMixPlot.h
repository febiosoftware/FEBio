#pragma once
#include "FECore/FEPlotData.h"

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
class FEPlotActualSoluteConcentration_old : public FEPlotDomainData
{
public:
	FEPlotActualSoluteConcentration_old(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);
	bool SetFilter(int nsol);
protected:
	int			m_nsol;
	FEModel*	m_pfem;
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
//! Base class for solute concentration variables
class FEPlotActualSolConcentration_ : public FEPlotDomainData
{
public:
	FEPlotActualSolConcentration_(FEModel* pfem, int nsol) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM), m_nsol(nsol) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	int	m_nsol;
};

//-----------------------------------------------------------------------------
// template class for instantiating solute concentration classes
template<int SOL> class FEPlotActualSolConcentrationT : public FEPlotActualSolConcentration_
{
public:
	FEPlotActualSolConcentrationT(FEModel* pfem) : FEPlotActualSolConcentration_(pfem, SOL) {}
};

//-----------------------------------------------------------------------------
//! Solute flux (for biphasic solute problems)
class FEPlotSoluteFlux_old : public FEPlotDomainData
{
public:
	FEPlotSoluteFlux_old(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);
	bool SetFilter(int nsol);
protected:
	int			m_nsol;
	FEModel*	m_pfem;
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
//! Base class for solute flux variables
class FEPlotSolFlux_ : public FEPlotDomainData
{
public:
	FEPlotSolFlux_(FEModel* pfem, int nsol) : FEPlotDomainData(PLT_VEC3F, FMT_ITEM), m_nsol(nsol) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	int	m_nsol;
};

//-----------------------------------------------------------------------------
// template class for instantiating solute flux classes
template<int SOL> class FEPlotSolFluxT : public FEPlotSolFlux_
{
public:
	FEPlotSolFluxT(FEModel* pfem) : FEPlotSolFlux_(pfem, SOL){}
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
class FEPlotSBMConcentration_old : public FEPlotDomainData
{
public:
	FEPlotSBMConcentration_old(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);
	bool SetFilter(int nsol);
protected:
	int			m_nsbm;
	FEModel*	m_pfem;
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
//! Base class for solid-bound molecule concentration
class FEPlotSBMConcentration_ : public FEPlotDomainData
{
public:
	FEPlotSBMConcentration_(int nsbm) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM), m_nsbm(nsbm) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	int m_nsbm;
};

//-----------------------------------------------------------------------------
// template class for instantiating solid-bound molecule variables
template <int SBM> class FEPlotSBMConcentrationT : public FEPlotSBMConcentration_
{
public:
	FEPlotSBMConcentrationT(FEModel* pfem) : FEPlotSBMConcentration_(SBM) {}
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
class FEPlotEffectiveSoluteConcentration_old : public FEPlotDomainData
{
public:
	FEPlotEffectiveSoluteConcentration_old(FEModel* pfem);
	bool SetFilter(const char* sz);
	bool SetFilter(int nsol);
	bool Save(FEDomain& m, FEDataStream& a);
protected:
	int			m_nsol;
	FEModel*	m_pfem;
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
//! Base class for nodal effective solute concentrations
class FEPlotEffectiveSolConcentration_ : public FEPlotDomainData
{
public:
	FEPlotEffectiveSolConcentration_(FEModel* pfem, int nsol) : FEPlotDomainData(PLT_FLOAT, FMT_NODE), m_nsol(nsol) {}
	bool Save(FEDomain& m, FEDataStream& a);
private:
	int m_nsol;
};

//-----------------------------------------------------------------------------
//! template class for instantiating nodal effective solute concentrations
template<int SOL> class FEPlotEffectiveSolConcentrationT : public FEPlotEffectiveSolConcentration_
{
public:
	FEPlotEffectiveSolConcentrationT(FEModel* pfem) : FEPlotEffectiveSolConcentration_(pfem, SOL) {}
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
//! Base class for nodal effective solute concentrations
class FEPlotEffectiveShellSolConcentration_ : public FEPlotDomainData
{
public:
	FEPlotEffectiveShellSolConcentration_(FEModel* pfem, int nsol) : FEPlotDomainData(PLT_FLOAT, FMT_NODE), m_nsol(nsol) {}
	bool Save(FEDomain& m, FEDataStream& a);
private:
	int m_nsol;
};

//-----------------------------------------------------------------------------
//! template class for instantiating nodal effective solute concentrations
template<int SOL> class FEPlotEffectiveShellSolConcentrationT : public FEPlotEffectiveShellSolConcentration_
{
public:
	FEPlotEffectiveShellSolConcentrationT(FEModel* pfem) : FEPlotEffectiveShellSolConcentration_(pfem, SOL) {}
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
class FEPlotSBMRefAppDensity_old : public FEPlotDomainData
{
public:
	FEPlotSBMRefAppDensity_old(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM), m_nsbm(0) {}
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);
	bool SetFilter(int nsol);
protected:
	int			m_nsbm;
	FEModel*	m_pfem;
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
//! Base class for solid-bound molecule referential apparent density
class FEPlotSBMRefAppDensity_ : public FEPlotDomainData
{
public:
	FEPlotSBMRefAppDensity_(int nsbm) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM), m_nsbm(nsbm) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	int m_nsbm;
};

//-----------------------------------------------------------------------------
// template class for instantiating solid-bound molecule variables
template <int SBM> class FEPlotSBMRefAppDensityT : public FEPlotSBMRefAppDensity_
{
public:
	FEPlotSBMRefAppDensityT(FEModel* pfem) : FEPlotSBMRefAppDensity_(SBM) {}
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

