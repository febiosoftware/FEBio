/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

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
    bool                m_binit;
    vector<vec3d>       m_area;
    
public:
	FEPlotMixtureFluidFlowRate(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_FLOAT, FMT_REGION){ m_binit = true; }
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
	FEPlotActualFluidPressure(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid flux
class FEPlotFluidFlux : public FEPlotDomainData
{
public:
	FEPlotFluidFlux(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal Fluid flux
class FEPlotNodalFluidFlux : public FEPlotDomainData
{
public:
	FEPlotNodalFluidFlux(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_MULT){}
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
    FEPlotOsmolarity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
class FEPlotSBMConcentration : public FEPlotDomainData
{
public:
	FEPlotSBMConcentration(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
protected:
	vector<int>	m_sbm;
};

//-----------------------------------------------------------------------------
//! Electric potential
class FEPlotElectricPotential : public FEPlotDomainData
{
public:
	FEPlotElectricPotential(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Current density
class FEPlotCurrentDensity : public FEPlotDomainData
{
public:
	FEPlotCurrentDensity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Referential solid volume fraction
class FEPlotReferentialSolidVolumeFraction : public FEPlotDomainData
{
public:
    FEPlotReferentialSolidVolumeFraction(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fixed charge density
class FEPlotFixedChargeDensity : public FEPlotDomainData
{
public:
	FEPlotFixedChargeDensity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Referential fixed charge density
class FEPlotReferentialFixedChargeDensity : public FEPlotDomainData
{
public:
	FEPlotReferentialFixedChargeDensity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal effective fluid pressures
class FEPlotEffectiveFluidPressure : public FEPlotDomainData
{
public:
	FEPlotEffectiveFluidPressure(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_NODE){}
	bool Save(FEDomain& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal effective downstream fluid pressures
class FEPlotEffectiveShellFluidPressure : public FEPlotDomainData
{
public:
	FEPlotEffectiveShellFluidPressure(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_NODE){}
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
};

//-----------------------------------------------------------------------------
//! Receptor-ligand complex concentration
class FEPlotReceptorLigandConcentration : public FEPlotDomainData
{
public:
	FEPlotReceptorLigandConcentration(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
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
};

//-----------------------------------------------------------------------------
//! effective elasticity
class FEPlotEffectiveElasticity : public FEPlotDomainData
{
public:
	FEPlotEffectiveElasticity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_TENS4FS, FMT_ITEM){}
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
	FEPlotFluidForce(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_VEC3F, FMT_REGION){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid force
//! 
class FEPlotFluidForce2 : public FEPlotSurfaceData
{
public:
	FEPlotFluidForce2(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_VEC3F, FMT_REGION){}
	bool Save(FESurface& surf, FEDataStream& a);
};


//-----------------------------------------------------------------------------
//! Fluid pressure gap
//!
class FEPlotPressureGap : public FEPlotSurfaceData
{
public:
	FEPlotPressureGap(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_FLOAT, FMT_MULT){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid load support
//!
class FEPlotFluidLoadSupport : public FEPlotSurfaceData
{
public:
    FEPlotFluidLoadSupport(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_FLOAT, FMT_REGION){}
    bool Save(FESurface& surf, FEDataStream& a);
};

