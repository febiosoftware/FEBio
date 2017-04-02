#pragma once
#include "FECore/FEPlotData.h"
#include <FECore/FEModel.h>

//=============================================================================
//                            N O D E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Nodal displacement
class FEPlotDisplacement : public FENodeData
{
public:
    FEPlotDisplacement(FEModel* pfem) : FENodeData(PLT_VEC3F, FMT_NODE){}
    bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal effective fluid pressures
class FEPlotFluidDilatation : public FENodeData
{
public:
    FEPlotFluidDilatation(FEModel* pfem) : FENodeData(PLT_FLOAT, FMT_NODE){}
    bool Save(FEMesh& m, FEDataStream& a);
};

//=============================================================================
//                         S U R F A C E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Fluid surface force
//!
class FEPlotFluidSurfaceForce : public FESurfaceData
{
private:
    FEModel*            m_pfem;
    bool                m_binit;
    vector<FEElement*>  m_elem;
    vector<vec3d>       m_area;
    
public:
    FEPlotFluidSurfaceForce(FEModel* pfem) : FESurfaceData(PLT_VEC3F, FMT_REGION){ m_pfem = pfem; m_binit = true; }
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid mass flux
//!
class FEPlotFluidMassFlux : public FESurfaceData
{
private:
    FEModel*            m_pfem;
    bool                m_binit;
    vector<FEElement*>  m_elem;
    
public:
    FEPlotFluidMassFlux(FEModel* pfem) : FESurfaceData(PLT_FLOAT, FMT_REGION){ m_pfem = pfem; m_binit = true; }
    bool Save(FESurface& surf, FEDataStream& a);
};

//=============================================================================
//							D O M A I N   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Element elastic fluid pressure
class FEPlotElasticFluidPressure : public FEDomainData
{
public:
	FEPlotElasticFluidPressure(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid volume ratio
class FEPlotFluidVolumeRatio : public FEDomainData
{
public:
    FEPlotFluidVolumeRatio(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid density
class FEPlotFluidDensity : public FEDomainData
{
public:
    FEPlotFluidDensity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid velocity
class FEPlotFluidVelocity : public FEDomainData
{
public:
    FEPlotFluidVelocity(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid acceleration
class FEPlotFluidAcceleration : public FEDomainData
{
public:
    FEPlotFluidAcceleration(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid vorticity
class FEPlotFluidVorticity : public FEDomainData
{
public:
    FEPlotFluidVorticity(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid stresses
class FEPlotElementFluidStress : public FEDomainData
{
public:
    FEPlotElementFluidStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid rate of deformation
class FEPlotElementFluidRateOfDef : public FEDomainData
{
public:
    FEPlotElementFluidRateOfDef(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid stress power
class FEPlotFluidStressPower : public FEDomainData
{
public:
    FEPlotFluidStressPower(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid shear viscosity
class FEPlotFluidShearViscosity : public FEDomainData
{
public:
    FEPlotFluidShearViscosity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

