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
#include "FEBioThermoFluid.h"
#include <FECore/FECoreKernel.h>
#include "FEThermoFluidSolver.h"
#include "FEThermoFluid.h"
#include "FEThermoFluidDomain3D.h"
#include "FEThermoFluidDomainFactory.h"
#include "FEFixedFluidTemperature.h"
#include "FEInitialFluidTemperature.h"
#include "FEInitialFluidPressureTemperature.h"
#include "FEPrescribedFluidTemperature.h"
#include "FEFluidHeatSupplyConst.h"
#include "FEFluidNormalHeatFlux.h"
#include "FEFluidNaturalHeatFlux.h"
#include "FEFluidConvectiveHeatFlux.h"
#include "FENewtonianThermoFluid.h"
#include "FENewtonianRealVapor.h"
#include "FEIdealGas.h"
#include "FERealGas.h"
#include "FERealVapor.h"
#include "FERealLiquid.h"
#include "FELinearThermoElasticFluid.h"
#include "FEThermoFluidPressureLoad.h"
#include "FETemperatureBackFlowStabilization.h"
#include "FEThermoFluidPressureBC.h"
#include "FEThermoFluidTemperatureBC.h"
#include "FEThermoFluidModule.h"
#include "FEThermoFluidAnalysis.h"
#include "FEBioThermoFluidPlot.h"
#include "FEThermoFluidProperty.h"
#include "FEThermalPropConst.h"
#include "FEThermalPropLnJvirial.h"
#include "FEThermalPropTempDpndnt.h"
#include "FEFluidFourierLaw.h"
#include "FEThermoElasticFluid.h"
#include "FEBioThermoFluidData.h"

const char* FEBioThermoFluid::GetVariableName(FEBioThermoFluid::THERMOFLUID_VARIABLE var)
{
    switch (var)
    {
    case DISPLACEMENT                : return "displacement"               ; break;
    case RELATIVE_FLUID_VELOCITY     : return "relative fluid velocity"    ; break;
    case RELATIVE_FLUID_ACCELERATION : return "relative fluid acceleration"; break;
    case FLUID_DILATATION            : return "fluid dilatation"           ; break;
    case FLUID_DILATATION_TDERIV     : return "fluid dilatation tderiv"    ; break;
    case TEMPERATURE                 : return "temperature"                ; break;
    case TEMPERATURE_TDERIV          : return "temperature tderiv"         ; break;
    }
    assert(false);
    return nullptr;
}

void FEBioThermoFluid::InitModule()
{
    FECoreKernel& febio = FECoreKernel::GetInstance();

    // register domain
    febio.RegisterDomain(new FEThermoFluidDomainFactory);

    // define the thermo-fluid module
    febio.CreateModule(new FEThermoFluidModule, "thermo-fluid",
                       "{"
                       "   \"title\" : \"Thermofluid\","
                       "   \"info\"  : \"Fluid analysis with heat transfer and thermodynamics.\""
                       "}");
    febio.AddModuleDependency("fluid");

    //-----------------------------------------------------------------------------
    // analysis classes (default type must match module name!)
    REGISTER_FECORE_CLASS(FEThermoFluidAnalysis, "thermo-fluid");

    //-----------------------------------------------------------------------------
    REGISTER_FECORE_CLASS(FEThermoFluidSolver, "thermo-fluid");

    REGISTER_FECORE_CLASS(FEThermoFluid, "thermo-fluid");

    REGISTER_FECORE_CLASS(FEThermoFluidDomain3D, "thermo-fluid-3D");

    //-----------------------------------------------------------------------------
    // initial conditions
    REGISTER_FECORE_CLASS(FEInitialFluidTemperature  , "initial fluid temperature");
    REGISTER_FECORE_CLASS(FEInitialFluidPressureTemperature  , "initial fluid pressure and temperature");

    //-----------------------------------------------------------------------------
    // boundary conditions
    REGISTER_FECORE_CLASS(FEFixedFluidTemperature       , "zero fluid temperature"      );
    REGISTER_FECORE_CLASS(FEPrescribedFluidTemperature  , "prescribed fluid temperature");
    REGISTER_FECORE_CLASS(FEThermoFluidPressureBC       , "fluid pressure");
    REGISTER_FECORE_CLASS(FEThermoFluidTemperatureBC    , "natural temperature");

    //-----------------------------------------------------------------------------
    // Surface loads
    REGISTER_FECORE_CLASS(FEFluidNormalHeatFlux, "fluid heat flux");
    REGISTER_FECORE_CLASS(FEFluidNaturalHeatFlux, "fluid natural heat flux");
    REGISTER_FECORE_CLASS(FEFluidConvectiveHeatFlux, "fluid convective heat flux");
    REGISTER_FECORE_CLASS(FETemperatureBackFlowStabilization, "temperature backflow stabilization");

    //-----------------------------------------------------------------------------
    // Body loads
    REGISTER_FECORE_CLASS(FEFluidHeatSupplyConst   , "constant fluid heat supply");
    
    //-----------------------------------------------------------------------------
    // Materials
    
    // viscous thermofluids
    REGISTER_FECORE_CLASS(FENewtonianThermoFluid, "Newtonian thermofluid");
    REGISTER_FECORE_CLASS(FENewtonianRealVapor, "Newtonian real vapor");

    // elastic fluids
    REGISTER_FECORE_CLASS(FEIdealGas   , "ideal gas"   );
    REGISTER_FECORE_CLASS(FERealGas    , "real gas"    );
    REGISTER_FECORE_CLASS(FERealVapor  , "real vapor"  );
    REGISTER_FECORE_CLASS(FERealLiquid , "real liquid" );
    REGISTER_FECORE_CLASS(FELinearThermoElasticFluid , "linear liquid" );

    // heat flux
    REGISTER_FECORE_CLASS(FEFluidFourierLaw, "Fourier's law");

    // thermofluid property
    REGISTER_FECORE_CLASS(FEThermalPropConst, "constant");
    REGISTER_FECORE_CLASS(FEThermalPropLnJvirial, "virial T-lnJ");
    REGISTER_FECORE_CLASS(FEThermalPropTempDpndnt, "temperature dependent");

    //-----------------------------------------------------------------------------
    // loads
    REGISTER_FECORE_CLASS(FEThermoFluidPressureLoad, "fluid pressure constraint");

    //-----------------------------------------------------------------------------
    // plot variables
	REGISTER_FECORE_CLASS(FEPlotThermoFluidTemperature, "fluid temperature");
	REGISTER_FECORE_CLASS(FEPlotNodalThermoFluidTemperature, "nodal fluid temperature");
	REGISTER_FECORE_CLASS(FEPlotThermoFluidPressureTangentTemperature, "fluid pressure tangent temperature");
	REGISTER_FECORE_CLASS(FEPlotThermoFluidRelativeThermalPecletNumber, "fluid relative thermal Peclet number");
	REGISTER_FECORE_CLASS(FEPlotThermoFluidIsochoricSpecificHeatCapacity, "fluid isochoric specific heat capacity");
	REGISTER_FECORE_CLASS(FEPlotThermoFluidIsobaricSpecificHeatCapacity, "fluid isobaric specific heat capacity");
	REGISTER_FECORE_CLASS(FEPlotThermoFluidHeatFlux, "fluid heat flux");
	REGISTER_FECORE_CLASS(FEPlotThermoFluidSpecificEntropy, "fluid specific entropy");
	REGISTER_FECORE_CLASS(FEPlotThermoFluidSpecificInternalEnergy, "fluid specific internal energy");
	REGISTER_FECORE_CLASS(FEPlotThermoFluidSpecificGaugeEnthalpy, "fluid specific gauge enthalpy");
	REGISTER_FECORE_CLASS(FEPlotThermoFluidSpecificFreeEnthalpy, "fluid specific free enthalpy");
    REGISTER_FECORE_CLASS(FEPlotThermoFluidThermalConductivity, "fluid thermal conductivity");
    REGISTER_FECORE_CLASS(FEPlotThermoFluidHeatSupplyDensity, "fluid heat supply density");
    
    //-----------------------------------------------------------------------------
    // log variables
    REGISTER_FECORE_CLASS(FENodeThermoFluidTemperature, "nodal fluid temperature");
    REGISTER_FECORE_CLASS(FELogThermoElasticFluidPressure, "fp");
    REGISTER_FECORE_CLASS(FELogThermoFluidVolumeRatio, "fJ");
    REGISTER_FECORE_CLASS(FELogThermoFluidDensity, "fd");
    REGISTER_FECORE_CLASS(FELogThermoFluidSpecificFreeEnergy, "af");
    REGISTER_FECORE_CLASS(FELogThermoFluidSpecificEntropy,    "sf");
    REGISTER_FECORE_CLASS(FELogThermoFluidSpecificInternalEnergy, "uf");
    REGISTER_FECORE_CLASS(FELogThermoFluidSpecificStrainEnergy, "wf");
    REGISTER_FECORE_CLASS(FELogThermoFluidSpecificEnthalpy, "hf");
    REGISTER_FECORE_CLASS(FELogThermoFluidIsochoricSpecificHeatCapacity, "cvf");
    REGISTER_FECORE_CLASS(FELogThermoFluidIsobaricSpecificHeatCapacity, "cpf");
    REGISTER_FECORE_CLASS(FELogThermoFluidTemperature, "Tf");

	febio.SetActiveModule(0);
}
