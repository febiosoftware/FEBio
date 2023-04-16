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



#include "stdafx.h"
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
#include "FEIdealGas.h"
#include "FERealGas.h"
#include "FERealLiquid.h"
#include "FEFluidConstantConductivity.h"
#include "FETempDependentConductivity.h"
#include "FEThermoFluidPressureLoad.h"
#include "FETemperatureBackFlowStabilization.h"
#include "FEThermoFluidPressureBC.h"
#include "FEThermoFluidTemperatureBC.h"
#include "FEFluidModule.h"
#include "FEThermoFluidAnalysis.h"
#include <FECore/FEModelUpdate.h>
#include <FECore/FETimeStepController.h>

//-----------------------------------------------------------------------------
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
    // analyis classes (default type must match module name!)
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
    REGISTER_FECORE_CLASS(FETemperatureBackFlowStabilization, "temperature backflow stabilization");

    //-----------------------------------------------------------------------------
    // Body loads
    REGISTER_FECORE_CLASS(FEFluidHeatSupplyConst   , "constant fluid heat supply");
    
    //-----------------------------------------------------------------------------
    // Materials
    REGISTER_FECORE_CLASS(FEIdealGas   , "ideal gas"   );
    REGISTER_FECORE_CLASS(FERealGas    , "real gas"    );
    REGISTER_FECORE_CLASS(FERealLiquid , "real liquid" );
    REGISTER_FECORE_CLASS(FEFluidConstantConductivity, "constant thermal conductivity");
    REGISTER_FECORE_CLASS(FETempDependentConductivity, "temp-dependent thermal conductivity");
    REGISTER_FECORE_CLASS(FEThermoFluidPressureLoad, "fluid pressure");

    //-----------------------------------------------------------------------------
    // Reset solver parameters to preferred default settings
    febio.OnCreateEvent(CallWhenCreating<FENewtonStrategy>([](FENewtonStrategy* pc) {
        pc->m_maxups = 50;
    }));
    
    febio.OnCreateEvent(CallWhenCreating<FETimeStepController>([](FETimeStepController* pc) {
        pc->m_iteopt = 50;
    }));
    
    febio.OnCreateEvent(CallWhenCreating<FEThermoFluidAnalysis>([](FEThermoFluidAnalysis* pc) {
        pc->m_nanalysis = FEThermoFluidAnalysis::DYNAMIC;
    }));
    
    febio.OnCreateEvent(CallWhenCreating<FENewtonSolver>([](FENewtonSolver* pc) {
        pc->m_maxref = 5;
        pc->m_Rmax = 1.0e+20;
        // turn off reform on each time step and diverge reform
        pc->m_breformtimestep = false;
        pc->m_bdivreform = false;
    }));
    
    febio.SetActiveModule(0);
}
