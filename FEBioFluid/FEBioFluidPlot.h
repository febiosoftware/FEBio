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
//! Nodal fluid velocity
class FEPlotNodalFluidVelocity : public FENodeData
{
public:
    FEPlotNodalFluidVelocity(FEModel* pfem) : FENodeData(PLT_VEC3F, FMT_NODE){}
    bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal relative fluid velocity
class FEPlotNodalRelativeFluidVelocity : public FENodeData
{
public:
    FEPlotNodalRelativeFluidVelocity(FEModel* pfem) : FENodeData(PLT_VEC3F, FMT_NODE){}
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
//! Fluid surface traction power
//!
class FEPlotFluidSurfaceTractionPower : public FESurfaceData
{
private:
    FEModel*            m_pfem;
    bool                m_binit;
    vector<FEElement*>  m_elem;
    vector<vec3d>       m_area;
    
public:
    FEPlotFluidSurfaceTractionPower(FEModel* pfem) : FESurfaceData(PLT_FLOAT, FMT_REGION){ m_pfem = pfem; m_binit = true; }
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid surface energy flux
//!
class FEPlotFluidSurfaceEnergyFlux : public FESurfaceData
{
private:
    FEModel*            m_pfem;
    bool                m_binit;
    vector<FEElement*>  m_elem;
    vector<vec3d>       m_area;
    
public:
    FEPlotFluidSurfaceEnergyFlux(FEModel* pfem) : FESurfaceData(PLT_FLOAT, FMT_REGION){ m_pfem = pfem; m_binit = true; }
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid mass flow rate
//!
class FEPlotFluidMassFlowRate : public FESurfaceData
{
private:
    FEModel*            m_pfem;
    bool                m_binit;
    vector<FEElement*>  m_elem;
    
public:
    FEPlotFluidMassFlowRate(FEModel* pfem) : FESurfaceData(PLT_FLOAT, FMT_REGION){ m_pfem = pfem; m_binit = true; }
    bool Save(FESurface& surf, FEDataStream& a);
};

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
class FEPlotFluidPressure : public FEDomainData
{
public:
	FEPlotFluidPressure(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

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
//! Element fluid density rate
class FEPlotFluidDensityRate : public FEDomainData
{
public:
    FEPlotFluidDensityRate(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
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
//! Element relative fluid velocity
class FEPlotRelativeFluidVelocity : public FEDomainData
{
public:
    FEPlotRelativeFluidVelocity(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM){}
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
//! Element fluid stress power density
class FEPlotFluidStressPowerDensity : public FEDomainData
{
public:
    FEPlotFluidStressPowerDensity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid heat supply density
class FEPlotFluidHeatSupplyDensity : public FEDomainData
{
public:
    FEPlotFluidHeatSupplyDensity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
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

//-----------------------------------------------------------------------------
//! Element strain energy density
class FEPlotFluidStrainEnergyDensity : public FEDomainData
{
public:
    FEPlotFluidStrainEnergyDensity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element kinetic energy density
class FEPlotFluidKineticEnergyDensity : public FEDomainData
{
public:
    FEPlotFluidKineticEnergyDensity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element energy density
class FEPlotFluidEnergyDensity : public FEDomainData
{
public:
    FEPlotFluidEnergyDensity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Strain energy
class FEPlotFluidElementStrainEnergy : public FEDomainData
{
public:
    FEPlotFluidElementStrainEnergy(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Kinetic energy
class FEPlotFluidElementKineticEnergy : public FEDomainData
{
public:
    FEPlotFluidElementKineticEnergy(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Center of mass
class FEPlotFluidElementCenterOfMass : public FEDomainData
{
public:
    FEPlotFluidElementCenterOfMass(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Linear momentum
class FEPlotFluidElementLinearMomentum : public FEDomainData
{
public:
    FEPlotFluidElementLinearMomentum(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Angular momentum
class FEPlotFluidElementAngularMomentum : public FEDomainData
{
public:
    FEPlotFluidElementAngularMomentum(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

