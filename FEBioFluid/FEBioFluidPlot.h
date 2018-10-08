#pragma once
#include "FECore/FEPlotData.h"
#include <FECore/FEModel.h>

//=============================================================================
//                            N O D E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Nodal displacement
class FEPlotDisplacement : public FEPlotNodeData
{
public:
    FEPlotDisplacement(FEModel* pfem) : FEPlotNodeData(PLT_VEC3F, FMT_NODE){}
    bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal fluid velocity
class FEPlotNodalFluidVelocity : public FEPlotNodeData
{
public:
    FEPlotNodalFluidVelocity(FEModel* pfem) : FEPlotNodeData(PLT_VEC3F, FMT_NODE){}
    bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal relative fluid velocity
class FEPlotNodalRelativeFluidVelocity : public FEPlotNodeData
{
public:
    FEPlotNodalRelativeFluidVelocity(FEModel* pfem) : FEPlotNodeData(PLT_VEC3F, FMT_NODE){}
    bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal effective fluid pressures
class FEPlotFluidDilatation : public FEPlotNodeData
{
public:
    FEPlotFluidDilatation(FEModel* pfem) : FEPlotNodeData(PLT_FLOAT, FMT_NODE){}
    bool Save(FEMesh& m, FEDataStream& a);
};

//=============================================================================
//                         S U R F A C E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Fluid surface force
//!
class FEPlotFluidSurfaceForce : public FEPlotSurfaceData
{
private:
    FEModel*            m_pfem;
    bool                m_binit;
    vector<vec3d>       m_area;
    
public:
    FEPlotFluidSurfaceForce(FEModel* pfem) : FEPlotSurfaceData(PLT_VEC3F, FMT_REGION){ m_pfem = pfem; m_binit = true; }
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid surface traction power
//!
class FEPlotFluidSurfaceTractionPower : public FEPlotSurfaceData
{
private:
    FEModel*            m_pfem;
    bool                m_binit;
    vector<vec3d>       m_area;
    
public:
    FEPlotFluidSurfaceTractionPower(FEModel* pfem) : FEPlotSurfaceData(PLT_FLOAT, FMT_REGION){ m_pfem = pfem; m_binit = true; }
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid surface energy flux
//!
class FEPlotFluidSurfaceEnergyFlux : public FEPlotSurfaceData
{
private:
    FEModel*            m_pfem;
    bool                m_binit;
    vector<vec3d>       m_area;
    
public:
    FEPlotFluidSurfaceEnergyFlux(FEModel* pfem) : FEPlotSurfaceData(PLT_FLOAT, FMT_REGION){ m_pfem = pfem; m_binit = true; }
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid mass flow rate
//!
class FEPlotFluidMassFlowRate : public FEPlotSurfaceData
{
private:
    FEModel*            m_pfem;
    
public:
    FEPlotFluidMassFlowRate(FEModel* pfem) : FEPlotSurfaceData(PLT_FLOAT, FMT_REGION){ m_pfem = pfem; }
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid flow rate
//!
class FEPlotFluidFlowRate : public FEPlotSurfaceData
{
private:
	FEModel*            m_pfem;
	bool                m_binit;
	vector<vec3d>       m_area;

public:
	FEPlotFluidFlowRate(FEModel* pfem) : FEPlotSurfaceData(PLT_FLOAT, FMT_REGION){ m_pfem = pfem; m_binit = true; }
	bool Save(FESurface& surf, FEDataStream& a);
};

//=============================================================================
//							D O M A I N   D A T A
//=============================================================================


//-----------------------------------------------------------------------------
//! Actual fluid pressure
class FEPlotFluidPressure : public FEPlotDomainData
{
public:
	FEPlotFluidPressure(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element elastic fluid pressure
class FEPlotElasticFluidPressure : public FEPlotDomainData
{
public:
	FEPlotElasticFluidPressure(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid temperature
class FEPlotFluidTemperature : public FEPlotDomainData
{
public:
    FEPlotFluidTemperature(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid volume ratio
class FEPlotFluidVolumeRatio : public FEPlotDomainData
{
public:
    FEPlotFluidVolumeRatio(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid density
class FEPlotFluidDensity : public FEPlotDomainData
{
public:
    FEPlotFluidDensity(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid density rate
class FEPlotFluidDensityRate : public FEPlotDomainData
{
public:
    FEPlotFluidDensityRate(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid velocity
class FEPlotFluidVelocity : public FEPlotDomainData
{
public:
    FEPlotFluidVelocity(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element relative fluid velocity
class FEPlotRelativeFluidVelocity : public FEPlotDomainData
{
public:
    FEPlotRelativeFluidVelocity(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid acceleration
class FEPlotFluidAcceleration : public FEPlotDomainData
{
public:
    FEPlotFluidAcceleration(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid vorticity
class FEPlotFluidVorticity : public FEPlotDomainData
{
public:
    FEPlotFluidVorticity(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid stresses
class FEPlotElementFluidStress : public FEPlotDomainData
{
public:
    FEPlotElementFluidStress(FEModel* pfem) : FEPlotDomainData(PLT_MAT3FS, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid rate of deformation
class FEPlotElementFluidRateOfDef : public FEPlotDomainData
{
public:
    FEPlotElementFluidRateOfDef(FEModel* pfem) : FEPlotDomainData(PLT_MAT3FS, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid stress power density
class FEPlotFluidStressPowerDensity : public FEPlotDomainData
{
public:
    FEPlotFluidStressPowerDensity(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid heat supply density
class FEPlotFluidHeatSupplyDensity : public FEPlotDomainData
{
public:
    FEPlotFluidHeatSupplyDensity(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid shear viscosity
class FEPlotFluidShearViscosity : public FEPlotDomainData
{
public:
    FEPlotFluidShearViscosity(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element strain energy density
class FEPlotFluidStrainEnergyDensity : public FEPlotDomainData
{
public:
    FEPlotFluidStrainEnergyDensity(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element kinetic energy density
class FEPlotFluidKineticEnergyDensity : public FEPlotDomainData
{
public:
    FEPlotFluidKineticEnergyDensity(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element energy density
class FEPlotFluidEnergyDensity : public FEPlotDomainData
{
public:
    FEPlotFluidEnergyDensity(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Strain energy
class FEPlotFluidElementStrainEnergy : public FEPlotDomainData
{
public:
    FEPlotFluidElementStrainEnergy(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Kinetic energy
class FEPlotFluidElementKineticEnergy : public FEPlotDomainData
{
public:
    FEPlotFluidElementKineticEnergy(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Center of mass
class FEPlotFluidElementCenterOfMass : public FEPlotDomainData
{
public:
    FEPlotFluidElementCenterOfMass(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Linear momentum
class FEPlotFluidElementLinearMomentum : public FEPlotDomainData
{
public:
    FEPlotFluidElementLinearMomentum(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Angular momentum
class FEPlotFluidElementAngularMomentum : public FEPlotDomainData
{
public:
    FEPlotFluidElementAngularMomentum(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

