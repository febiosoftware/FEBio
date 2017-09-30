#pragma once
#include "FECore/FEPlotData.h"

//=============================================================================
//                            N O D E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Nodal effective fluid pressures
class FEPlotEffectiveFluidPressure : public FENodeData
{
public:
    FEPlotEffectiveFluidPressure(FEModel* pfem) : FENodeData(PLT_FLOAT, FMT_NODE){}
    bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal effective downstream fluid pressures
class FEPlotEffectiveShellFluidPressure : public FENodeData
{
public:
    FEPlotEffectiveShellFluidPressure(FEModel* pfem) : FENodeData(PLT_FLOAT, FMT_NODE){}
    bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal effective solute concentrations (for biphasic-solute problems)
class FEPlotEffectiveSoluteConcentration : public FENodeData
{
public:
    FEPlotEffectiveSoluteConcentration(FEModel* pfem);
    bool SetFilter(const char* sz);
    bool SetFilter(int nsol);
    bool Save(FEMesh& m, FEDataStream& a);
protected:
    int            m_nsol;
    FEModel*    m_pfem;
};

//-----------------------------------------------------------------------------
//! Base class for nodal effective solute concentrations
class FEPlotEffectiveSolConcentration_ : public FENodeData
{
public:
    FEPlotEffectiveSolConcentration_(FEModel* pfem, int nsol) : FENodeData(PLT_FLOAT, FMT_NODE), m_nsol(nsol) {}
    bool Save(FEMesh& m, FEDataStream& a);
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
class FEPlotEffectiveShellSoluteConcentration : public FENodeData
{
public:
    FEPlotEffectiveShellSoluteConcentration(FEModel* pfem);
    bool SetFilter(const char* sz);
    bool SetFilter(int nsol);
    bool Save(FEMesh& m, FEDataStream& a);
protected:
    int            m_nsol;
    FEModel*    m_pfem;
};

//-----------------------------------------------------------------------------
//! Base class for nodal effective solute concentrations
class FEPlotEffectiveShellSolConcentration_ : public FENodeData
{
public:
    FEPlotEffectiveShellSolConcentration_(FEModel* pfem, int nsol) : FENodeData(PLT_FLOAT, FMT_NODE), m_nsol(nsol) {}
    bool Save(FEMesh& m, FEDataStream& a);
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

//=============================================================================
//                         S U R F A C E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Fluid flow rate
//!
class FEPlotFluidFlowRate : public FESurfaceData
{
private:
    FEModel*            m_pfem;
    bool                m_binit;
    vector<FEElement*>  m_elem;
    vector<vec3d>       m_area;
    
public:
    FEPlotFluidFlowRate(FEModel* pfem) : FESurfaceData(PLT_FLOAT, FMT_REGION){ m_pfem = pfem; m_binit = true; }
    bool Save(FESurface& surf, FEDataStream& a);
};

//=============================================================================
//							D O M A I N   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Actual fluid pressure
class FEPlotActualFluidPressure : public FEDomainData
{
public:
	FEPlotActualFluidPressure(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid flux
class FEPlotFluidFlux : public FEDomainData
{
public:
	FEPlotFluidFlux(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal Fluid flux
class FEPlotNodalFluidFlux : public FEDomainData
{
public:
	FEPlotNodalFluidFlux(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_MULT){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Actual solute concentration
class FEPlotActualSoluteConcentration : public FEDomainData
{
public:
	FEPlotActualSoluteConcentration(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);
	bool SetFilter(int nsol);
protected:
	int			m_nsol;
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
//! Base class for solute concentration variables
class FEPlotActualSolConcentration_ : public FEDomainData
{
public:
	FEPlotActualSolConcentration_(FEModel* pfem, int nsol) : FEDomainData(PLT_FLOAT, FMT_ITEM), m_nsol(nsol) {}
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
class FEPlotSoluteFlux : public FEDomainData
{
public:
	FEPlotSoluteFlux(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);
	bool SetFilter(int nsol);
protected:
	int			m_nsol;
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
//! Base class for solute flux variables
class FEPlotSolFlux_ : public FEDomainData
{
public:
	FEPlotSolFlux_(FEModel* pfem, int nsol) : FEDomainData(PLT_VEC3F, FMT_ITEM), m_nsol(nsol) {}
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
class FEPlotOsmolarity : public FEDomainData
{
public:
    FEPlotOsmolarity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
class FEPlotSBMConcentration : public FEDomainData
{
public:
	FEPlotSBMConcentration(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);
	bool SetFilter(int nsol);
protected:
	int			m_nsbm;
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
//! Base class for solid-bound molecule concentration
class FEPlotSBMConcentration_ : public FEDomainData
{
public:
	FEPlotSBMConcentration_(int nsbm) : FEDomainData(PLT_FLOAT, FMT_ITEM), m_nsbm(nsbm) {}
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
class FEPlotElectricPotential : public FEDomainData
{
public:
	FEPlotElectricPotential(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Current density
class FEPlotCurrentDensity : public FEDomainData
{
public:
	FEPlotCurrentDensity(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Referential solid volume fraction
class FEPlotReferentialSolidVolumeFraction : public FEDomainData
{
public:
    FEPlotReferentialSolidVolumeFraction(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fixed charge density
class FEPlotFixedChargeDensity : public FEDomainData
{
public:
	FEPlotFixedChargeDensity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Referential fixed charge density
class FEPlotReferentialFixedChargeDensity : public FEDomainData
{
public:
	FEPlotReferentialFixedChargeDensity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Receptor-ligand complex concentration
class FEPlotReceptorLigandConcentration : public FEDomainData
{
public:
	FEPlotReceptorLigandConcentration(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Base class for solid-bound molecule referential apparent density
class FEPlotSBMRefAppDensity : public FEDomainData
{
public:
	FEPlotSBMRefAppDensity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM), m_nsbm(0) {}
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);
	bool SetFilter(int nsol);
protected:
	int			m_nsbm;
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
//! Base class for solid-bound molecule referential apparent density
class FEPlotSBMRefAppDensity_ : public FEDomainData
{
public:
	FEPlotSBMRefAppDensity_(int nsbm) : FEDomainData(PLT_FLOAT, FMT_ITEM), m_nsbm(nsbm) {}
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
class FEPlotEffectiveElasticity : public FEDomainData
{
public:
	FEPlotEffectiveElasticity(FEModel* pfem) : FEDomainData(PLT_TENS4FS, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//=============================================================================
//                         S U R F A C E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Fluid force
//!
class FEPlotFluidForce : public FESurfaceData
{
public:
	FEPlotFluidForce(FEModel* pfem) : FESurfaceData(PLT_VEC3F, FMT_REGION){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid force
//! 
class FEPlotFluidForce2 : public FESurfaceData
{
public:
	FEPlotFluidForce2(FEModel* pfem) : FESurfaceData(PLT_VEC3F, FMT_REGION){}
	bool Save(FESurface& surf, FEDataStream& a);
};


//-----------------------------------------------------------------------------
//! Fluid pressure gap
//!
class FEPlotPressureGap : public FESurfaceData
{
public:
	FEPlotPressureGap(FEModel* pfem) : FESurfaceData(PLT_FLOAT, FMT_MULT){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid load support
//!
class FEPlotFluidLoadSupport : public FESurfaceData
{
public:
    FEPlotFluidLoadSupport(FEModel* pfem) : FESurfaceData(PLT_FLOAT, FMT_REGION){}
    bool Save(FESurface& surf, FEDataStream& a);
};

