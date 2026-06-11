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
#include <FECore/FEPlotData.h>

//! Element fluid temperature
class FEPlotThermoFluidTemperature : public FEPlotDomainData
{
public:
	FEPlotThermoFluidTemperature(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_RELATIVE_TEMPERATURE); }
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal fluid temperature
class FEPlotNodalThermoFluidTemperature : public FEPlotNodeData
{
public:
	FEPlotNodalThermoFluidTemperature(FEModel* pfem) : FEPlotNodeData(pfem, PLT_FLOAT, FMT_NODE) { SetUnits(UNIT_RELATIVE_TEMPERATURE); }
	bool Save(FEMesh& m, FEDataStream& a);
};

//! Fluid pressure tangent temperature
class FEPlotThermoFluidPressureTangentTemperature : public FEPlotDomainData
{
public:
	FEPlotThermoFluidPressureTangentTemperature(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_SPECIFIC_ENTROPY); }
	bool Save(FEDomain& dom, FEDataStream& a);
};

//! Element fluid heat flux
class FEPlotThermoFluidHeatFlux : public FEPlotDomainData
{
public:
	FEPlotThermoFluidHeatFlux(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM) { SetUnits(UNIT_ENERGY_FLUX); }
	bool Save(FEDomain& dom, FEDataStream& a);
};

//! Element relative thermal Peclet number
class FEPlotThermoFluidRelativeThermalPecletNumber : public FEPlotDomainData
{
public:
	FEPlotThermoFluidRelativeThermalPecletNumber(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_RECIPROCAL_LENGTH); }
	bool Save(FEDomain& dom, FEDataStream& a);
};

//! Thermal conductivity
class FEPlotThermoFluidThermalConductivity : public FEPlotDomainData
{
public:
	FEPlotThermoFluidThermalConductivity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_THERMAL_CONDUCTIVITY); }
	bool Save(FEDomain& dom, FEDataStream& a);
};

//! Specific isochoric heat capacity
class FEPlotThermoFluidIsochoricSpecificHeatCapacity : public FEPlotDomainData
{
public:
	FEPlotThermoFluidIsochoricSpecificHeatCapacity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_SPECIFIC_ENTROPY); }
	bool Save(FEDomain& dom, FEDataStream& a);
};

//! Specific isobaric heat capacity
class FEPlotThermoFluidIsobaricSpecificHeatCapacity : public FEPlotDomainData
{
public:
	FEPlotThermoFluidIsobaricSpecificHeatCapacity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_SPECIFIC_ENTROPY); }
	bool Save(FEDomain& dom, FEDataStream& a);
};

//! Specific free enthalpy
class FEPlotThermoFluidSpecificFreeEnthalpy : public FEPlotDomainData
{
public:
	FEPlotThermoFluidSpecificFreeEnthalpy(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_SPECIFIC_ENERGY); }
	bool Save(FEDomain& dom, FEDataStream& a);
};

//! Specific gauge enthalpy
class FEPlotThermoFluidSpecificGaugeEnthalpy : public FEPlotDomainData
{
public:
	FEPlotThermoFluidSpecificGaugeEnthalpy(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_SPECIFIC_ENERGY); }
	bool Save(FEDomain& dom, FEDataStream& a);
};

//! Specific internal energy
class FEPlotThermoFluidSpecificInternalEnergy : public FEPlotDomainData
{
public:
	FEPlotThermoFluidSpecificInternalEnergy(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_SPECIFIC_ENERGY); }
	bool Save(FEDomain& dom, FEDataStream& a);
};

//! Specific entropy
class FEPlotThermoFluidSpecificEntropy : public FEPlotDomainData
{
public:
	FEPlotThermoFluidSpecificEntropy(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_SPECIFIC_ENTROPY); }
	bool Save(FEDomain& dom, FEDataStream& a);
};

//! Specific free energy
class FEPlotThermoFluidSpecificFreeEnergy : public FEPlotDomainData
{
public:
    FEPlotThermoFluidSpecificFreeEnergy(FEModel* fem) : FEPlotDomainData(fem, PLT_FLOAT, FMT_ITEM){
        SetUnits(UNIT_SPECIFIC_ENERGY); }
    bool Save(FEDomain& dom, FEDataStream& a) override;
};

//-----------------------------------------------------------------------------
//! Element ThermoFluid heat supply density
class FEPlotThermoFluidHeatSupplyDensity : public FEPlotDomainData
{
public:
    FEPlotThermoFluidHeatSupplyDensity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_POWER_DENSITY); }
    bool Save(FEDomain& dom, FEDataStream& a);
};

