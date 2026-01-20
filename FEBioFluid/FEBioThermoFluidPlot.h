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
class FEPlotFluidTemperature : public FEPlotDomainData
{
public:
	FEPlotFluidTemperature(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_RELATIVE_TEMPERATURE); }
	bool Save(FEDomain& dom, FEDataStream& a);
};

//! Nodal fluid temperature
class FEPlotNodalFluidTemperature : public FEPlotNodeData
{
public:
	FEPlotNodalFluidTemperature(FEModel* pfem) : FEPlotNodeData(pfem, PLT_FLOAT, FMT_NODE) { SetUnits(UNIT_RELATIVE_TEMPERATURE); }
	bool Save(FEMesh& m, FEDataStream& a);
};

//! Fluid pressure tangent temperature
class FEPlotFluidPressureTangentTemperature : public FEPlotDomainData
{
public:
	FEPlotFluidPressureTangentTemperature(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_SPECIFIC_ENTROPY); }
	bool Save(FEDomain& dom, FEDataStream& a);
};

//! Element fluid heat flux
class FEPlotFluidHeatFlux : public FEPlotDomainData
{
public:
	FEPlotFluidHeatFlux(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM) { SetUnits(UNIT_ENERGY_FLUX); }
	bool Save(FEDomain& dom, FEDataStream& a);
};

//! Element relative thermal Peclet number
class FEPlotFluidRelativeThermalPecletNumber : public FEPlotDomainData
{
public:
	FEPlotFluidRelativeThermalPecletNumber(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_RECIPROCAL_LENGTH); }
	bool Save(FEDomain& dom, FEDataStream& a);
};

//! Thermal conductivity
class FEPlotFluidThermalConductivity : public FEPlotDomainData
{
public:
	FEPlotFluidThermalConductivity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_THERMAL_CONDUCTIVITY); }
	bool Save(FEDomain& dom, FEDataStream& a);
};

//! Specific isochoric heat capacity
class FEPlotFluidIsochoricSpecificHeatCapacity : public FEPlotDomainData
{
public:
	FEPlotFluidIsochoricSpecificHeatCapacity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_SPECIFIC_ENTROPY); }
	bool Save(FEDomain& dom, FEDataStream& a);
};

//! Specific isobaric heat capacity
class FEPlotFluidIsobaricSpecificHeatCapacity : public FEPlotDomainData
{
public:
	FEPlotFluidIsobaricSpecificHeatCapacity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_SPECIFIC_ENTROPY); }
	bool Save(FEDomain& dom, FEDataStream& a);
};
