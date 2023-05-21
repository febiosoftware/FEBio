/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

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
#include "FECore/FEPlotData.h"
#include <FECore/units.h>

//=============================================================================
//                            N O D E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Nodal displacement
class FEPlotDisplacement : public FEPlotNodeData
{
public:
    FEPlotDisplacement(FEModel* pfem) : FEPlotNodeData(pfem, PLT_VEC3F, FMT_NODE) { SetUnits(UNIT_LENGTH); }
    bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal fluid velocity
class FEPlotNodalFluidVelocity : public FEPlotNodeData
{
public:
    FEPlotNodalFluidVelocity(FEModel* pfem) : FEPlotNodeData(pfem, PLT_VEC3F, FMT_NODE) { SetUnits(UNIT_VELOCITY); }
    bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal relative fluid velocity
class FEPlotNodalRelativeFluidVelocity : public FEPlotNodeData
{
public:
    FEPlotNodalRelativeFluidVelocity(FEModel* pfem) : FEPlotNodeData(pfem, PLT_VEC3F, FMT_NODE) { SetUnits(UNIT_VELOCITY); }
    bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal fluid dilatation
class FEPlotFluidDilatation : public FEPlotNodeData
{
public:
    FEPlotFluidDilatation(FEModel* pfem) : FEPlotNodeData(pfem, PLT_FLOAT, FMT_NODE){}
    bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal effective fluid pressures
class FEPlotFluidEffectivePressure : public FEPlotDomainData
{
public:
    FEPlotFluidEffectivePressure(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_NODE) { SetUnits(UNIT_PRESSURE); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal polar fluid angular velocity
class FEPlotNodalPolarFluidAngularVelocity : public FEPlotNodeData
{
public:
    FEPlotNodalPolarFluidAngularVelocity(FEModel* pfem) : FEPlotNodeData(pfem, PLT_VEC3F, FMT_NODE) { SetUnits(UNIT_ANGULAR_VELOCITY); }
    bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal fluid temperature
class FEPlotNodalFluidTemperature : public FEPlotNodeData
{
public:
    FEPlotNodalFluidTemperature(FEModel* pfem) : FEPlotNodeData(pfem, PLT_FLOAT, FMT_NODE) { SetUnits(UNIT_RELATIVE_TEMPERATURE); }
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
    bool                m_binit;
    vector<vec3d>       m_area;
    
public:
    FEPlotFluidSurfaceForce(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_VEC3F, FMT_REGION){ m_binit = true; SetUnits(UNIT_FORCE); }
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid surface moment (polar fluids)
//!
class FEPlotFluidSurfaceMoment : public FEPlotSurfaceData
{
private:
    bool                m_binit;
    vector<vec3d>       m_area;
    
public:
    FEPlotFluidSurfaceMoment(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_VEC3F, FMT_REGION){ m_binit = true; SetUnits(UNIT_MOMENT); }
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid surface pressure
//!
class FEPlotFluidSurfacePressure : public FEPlotSurfaceData
{
public:
    FEPlotFluidSurfacePressure(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_PRESSURE); }
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid surface traction power
//!
class FEPlotFluidSurfaceTractionPower : public FEPlotSurfaceData
{
private:
    bool                m_binit;
    vector<vec3d>       m_area;
    
public:
    FEPlotFluidSurfaceTractionPower(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_FLOAT, FMT_REGION){ m_binit = true; SetUnits(UNIT_POWER); }
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid surface energy flux
//!
class FEPlotFluidSurfaceEnergyFlux : public FEPlotSurfaceData
{
private:
    bool                m_binit;
    vector<vec3d>       m_area;
    
public:
    FEPlotFluidSurfaceEnergyFlux(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_FLOAT, FMT_REGION){ m_binit = true; SetUnits(UNIT_POWER); }
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid mass flow rate
//!
class FEPlotFluidMassFlowRate : public FEPlotSurfaceData
{
public:
    FEPlotFluidMassFlowRate(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_FLOAT, FMT_REGION) { SetUnits(UNIT_MASS_FLOW_RATE); }
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fluid flow rate
//!
class FEPlotFluidFlowRate : public FEPlotSurfaceData
{
private:
	bool                m_binit;
	vector<vec3d>       m_area;

public:
	FEPlotFluidFlowRate(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_FLOAT, FMT_REGION){ m_binit = true; SetUnits(UNIT_FLOW_RATE); }
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
	FEPlotFluidPressure(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_PRESSURE); }
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element elastic fluid pressure
class FEPlotElasticFluidPressure : public FEPlotDomainData
{
public:
	FEPlotElasticFluidPressure(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_PRESSURE); }
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid temperature
class FEPlotFluidTemperature : public FEPlotDomainData
{
public:
    FEPlotFluidTemperature(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_RELATIVE_TEMPERATURE); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid volume ratio
class FEPlotFluidVolumeRatio : public FEPlotDomainData
{
public:
    FEPlotFluidVolumeRatio(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid density
class FEPlotFluidDensity : public FEPlotDomainData
{
public:
    FEPlotFluidDensity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_DENSITY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid density rate
class FEPlotFluidDensityRate : public FEPlotDomainData
{
public:
    FEPlotFluidDensityRate(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_DENSITY_RATE); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid velocity
class FEPlotFluidVelocity : public FEPlotDomainData
{
public:
    FEPlotFluidVelocity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM) { SetUnits(UNIT_VELOCITY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element relative fluid velocity
class FEPlotRelativeFluidVelocity : public FEPlotDomainData
{
public:
    FEPlotRelativeFluidVelocity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM) { SetUnits(UNIT_VELOCITY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element relative fluid flux
class FEPlotFSIFluidFlux : public FEPlotDomainData
{
public:
    FEPlotFSIFluidFlux(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM) { SetUnits(UNIT_VELOCITY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Permeability
class FEPlotPermeability : public FEPlotDomainData
{
public:
    FEPlotPermeability(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3FS, FMT_ITEM) { SetUnits(UNIT_PERMEABILITY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element GradJ
class FEPlotGradJ : public FEPlotDomainData
{
public:
    FEPlotGradJ(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element gradphif
class FEPlotGradPhiF : public FEPlotDomainData
{
public:
    FEPlotGradPhiF(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element porosity
class FEPlotBFSIPorosity : public FEPlotDomainData
{
public:
    FEPlotBFSIPorosity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element solid volume fraction
class FEPlotBFSISolidVolumeFraction : public FEPlotDomainData
{
public:
    FEPlotBFSISolidVolumeFraction(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid acceleration
class FEPlotFluidAcceleration : public FEPlotDomainData
{
public:
    FEPlotFluidAcceleration(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM) { SetUnits(UNIT_ACCELERATION); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid vorticity
class FEPlotFluidVorticity : public FEPlotDomainData
{
public:
    FEPlotFluidVorticity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM) { SetUnits(UNIT_ANGULAR_VELOCITY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element polar fluid angular velocity
class FEPlotPolarFluidAngularVelocity : public FEPlotDomainData
{
public:
    FEPlotPolarFluidAngularVelocity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM) { SetUnits(UNIT_ANGULAR_VELOCITY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element polar fluid relative angular velocity
class FEPlotPolarFluidRelativeAngularVelocity : public FEPlotDomainData
{
public:
    FEPlotPolarFluidRelativeAngularVelocity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM) { SetUnits(UNIT_ANGULAR_VELOCITY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element polar fluid rregional angular velocity
class FEPlotPolarFluidRegionalAngularVelocity : public FEPlotDomainData
{
public:
    FEPlotPolarFluidRegionalAngularVelocity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM) { SetUnits(UNIT_ANGULAR_VELOCITY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid heat flux
class FEPlotFluidHeatFlux : public FEPlotDomainData
{
public:
    FEPlotFluidHeatFlux(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM) { SetUnits(UNIT_ENERGY_FLUX); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid stresses
class FEPlotFluidStress : public FEPlotDomainData
{
public:
    FEPlotFluidStress(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3FS, FMT_ITEM) { SetUnits(UNIT_PRESSURE); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid rate of deformation
class FEPlotElementFluidRateOfDef : public FEPlotDomainData
{
public:
    FEPlotElementFluidRateOfDef(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3FS, FMT_ITEM) { SetUnits(UNIT_RECIPROCAL_TIME); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid stress power density
class FEPlotFluidStressPowerDensity : public FEPlotDomainData
{
public:
    FEPlotFluidStressPowerDensity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_POWER_DENSITY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid heat supply density
class FEPlotFluidHeatSupplyDensity : public FEPlotDomainData
{
public:
    FEPlotFluidHeatSupplyDensity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_POWER_DENSITY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid shear viscosity
class FEPlotFluidShearViscosity : public FEPlotDomainData
{
public:
    FEPlotFluidShearViscosity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_VISCOSITY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element strain energy density
class FEPlotFluidStrainEnergyDensity : public FEPlotDomainData
{
public:
    FEPlotFluidStrainEnergyDensity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_ENERGY_DENSITY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element kinetic energy density
class FEPlotFluidKineticEnergyDensity : public FEPlotDomainData
{
public:
    FEPlotFluidKineticEnergyDensity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_ENERGY_DENSITY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element energy density
class FEPlotFluidEnergyDensity : public FEPlotDomainData
{
public:
    FEPlotFluidEnergyDensity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_ENERGY_DENSITY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element fluid bulk modulus
class FEPlotFluidBulkModulus : public FEPlotDomainData
{
public:
    FEPlotFluidBulkModulus(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_PRESSURE); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Strain energy
class FEPlotFluidElementStrainEnergy : public FEPlotDomainData
{
public:
    FEPlotFluidElementStrainEnergy(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_ENERGY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Kinetic energy
class FEPlotFluidElementKineticEnergy : public FEPlotDomainData
{
public:
    FEPlotFluidElementKineticEnergy(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_ENERGY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Center of mass
class FEPlotFluidElementCenterOfMass : public FEPlotDomainData
{
public:
    FEPlotFluidElementCenterOfMass(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM) { SetUnits(UNIT_LENGTH); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Linear momentum
class FEPlotFluidElementLinearMomentum : public FEPlotDomainData
{
public:
    FEPlotFluidElementLinearMomentum(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM) { SetUnits(UNIT_LINEAR_MOMENTUM); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Angular momentum
class FEPlotFluidElementAngularMomentum : public FEPlotDomainData
{
public:
    FEPlotFluidElementAngularMomentum(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM) { SetUnits(UNIT_ANGULAR_MOMENTUM); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Specific free energy
class FEPlotFluidSpecificFreeEnergy : public FEPlotDomainData
{
public:
    FEPlotFluidSpecificFreeEnergy(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_SPECIFIC_ENERGY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Specific entropy
class FEPlotFluidSpecificEntropy : public FEPlotDomainData
{
public:
    FEPlotFluidSpecificEntropy(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_SPECIFIC_ENTROPY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Specific internal energy
class FEPlotFluidSpecificInternalEnergy : public FEPlotDomainData
{
public:
    FEPlotFluidSpecificInternalEnergy(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_SPECIFIC_ENERGY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Specific gage enthalpy
class FEPlotFluidSpecificGageEnthalpy : public FEPlotDomainData
{
public:
    FEPlotFluidSpecificGageEnthalpy(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_SPECIFIC_ENERGY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Specific free enthalpy
class FEPlotFluidSpecificFreeEnthalpy : public FEPlotDomainData
{
public:
    FEPlotFluidSpecificFreeEnthalpy(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_SPECIFIC_ENERGY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Specific strain energy
class FEPlotFluidSpecificStrainEnergy : public FEPlotDomainData
{
public:
    FEPlotFluidSpecificStrainEnergy(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_SPECIFIC_ENERGY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Specific isochoric heat capacity
class FEPlotFluidIsochoricSpecificHeatCapacity : public FEPlotDomainData
{
public:
    FEPlotFluidIsochoricSpecificHeatCapacity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_SPECIFIC_ENTROPY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Specific isobaric heat capacity
class FEPlotFluidIsobaricSpecificHeatCapacity : public FEPlotDomainData
{
public:
    FEPlotFluidIsobaricSpecificHeatCapacity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_SPECIFIC_ENTROPY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Thermal conductivity
class FEPlotFluidThermalConductivity : public FEPlotDomainData
{
public:
    FEPlotFluidThermalConductivity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_THERMAL_CONDUCTIVITY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element solid stresses
class FEPlotFSISolidStress : public FEPlotDomainData
{
public:
    FEPlotFSISolidStress(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3FS, FMT_ITEM) { SetUnits(UNIT_PRESSURE); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! fluid shear stress error
class FEPlotFluidShearStressError : public FEPlotDomainData
{
public:
	FEPlotFluidShearStressError(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element polar fluid  stresses
class FEPlotPolarFluidStress : public FEPlotDomainData
{
public:
    FEPlotPolarFluidStress(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3F, FMT_ITEM) { SetUnits(UNIT_PRESSURE); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element polar fluid couple stresses
class FEPlotPolarFluidCoupleStress : public FEPlotDomainData
{
public:
    FEPlotPolarFluidCoupleStress(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3F, FMT_ITEM) { SetUnits(UNIT_ENERGY_AREAL_DENSITY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element relative Reynolds number
class FEPlotFluidRelativeReynoldsNumber : public FEPlotDomainData
{
public:
    FEPlotFluidRelativeReynoldsNumber(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_RECIPROCAL_LENGTH); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element relative Peclet number
class FEPlotFluidRelativePecletNumber : public FEPlotDomainData
{
public:
    FEPlotFluidRelativePecletNumber(FEModel* pfem);
    bool Save(FEDomain& dom, FEDataStream& a);
protected:
    vector<int>    m_sol;
};

