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
#include <FECore/NodeDataRecord.h>
#include <FECore/ElementDataRecord.h>

//=============================================================================
// N O D E  D A T A
//=============================================================================

//-----------------------------------------------------------------------------
class FENodeThermoFluidTemperature : public FELogNodeData
{
public:
    FENodeThermoFluidTemperature(FEModel* pfem) : FELogNodeData(pfem){}
    double value(const FENode& node) override;
};

//=============================================================================
// E L E M E N T   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
class FELogThermoElasticFluidPressure : public FELogElemData
{
public:
    FELogThermoElasticFluidPressure(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogThermoFluidVolumeRatio : public FELogElemData
{
public:
    FELogThermoFluidVolumeRatio(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogThermoFluidDensity : public FELogElemData
{
public:
    FELogThermoFluidDensity(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogThermoFluidSpecificFreeEnergy : public FELogElemData
{
public:
    FELogThermoFluidSpecificFreeEnergy(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogThermoFluidSpecificEntropy : public FELogElemData
{
public:
    FELogThermoFluidSpecificEntropy(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogThermoFluidSpecificInternalEnergy : public FELogElemData
{
public:
    FELogThermoFluidSpecificInternalEnergy(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogThermoFluidSpecificEnthalpy : public FELogElemData
{
public:
    FELogThermoFluidSpecificEnthalpy(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogThermoFluidSpecificStrainEnergy : public FELogElemData
{
public:
    FELogThermoFluidSpecificStrainEnergy(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogThermoFluidIsochoricSpecificHeatCapacity : public FELogElemData
{
public:
    FELogThermoFluidIsochoricSpecificHeatCapacity(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogThermoFluidIsobaricSpecificHeatCapacity : public FELogElemData
{
public:
    FELogThermoFluidIsobaricSpecificHeatCapacity(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogThermoFluidTemperature : public FELogElemData
{
public:
    FELogThermoFluidTemperature(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};
