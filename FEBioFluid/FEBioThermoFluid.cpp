/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in
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
#include "FEFluidNormalHeatFlux.h"
#include "FEIdealGas.h"
#include "FEIdealLiquid.h"
#include "FEFluidConstantConductivity.h"
#include "FEThermoFluidPressureLoad.h"
#include "FETemperatureBackFlowStabilization.h"

//-----------------------------------------------------------------------------
const char* FEBioThermoFluid::GetVariableName(FEBioThermoFluid::THERMOFLUID_VARIABLE var)
{
    switch (var)
    {
    case DISPLACEMENT                : return "displacement"               ; break;
    case RELATIVE_FLUID_VELOCITY     : return "relative fluid velocity"    ; break;
    case RELATIVE_FLUID_ACCELERATION : return "relative fluid acceleration"; break;
    case FLUID_DILATATION            : return "fluid dilation"             ; break;
    case FLUID_DILATATION_TDERIV     : return "fluid dilation tderiv"      ; break;
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
    febio.CreateModule("thermo-fluid");
    febio.SetModuleDependency("fluid");

    REGISTER_FECORE_CLASS(FEThermoFluidSolver, "thermo-fluid");

    REGISTER_FECORE_CLASS(FEThermoFluid, "thermo-fluid");

    REGISTER_FECORE_CLASS(FEThermoFluidDomain3D, "thermo-fluid-3D");

    REGISTER_FECORE_CLASS(FEFluidNormalHeatFlux, "fluid heat flux");
    REGISTER_FECORE_CLASS(FETemperatureBackFlowStabilization, "temperature backflow stabilization");

    REGISTER_FECORE_CLASS(FEIdealGas   , "ideal gas"   );
    REGISTER_FECORE_CLASS(FEIdealLiquid, "ideal liquid");
    REGISTER_FECORE_CLASS(FEFluidConstantConductivity, "constant thermal conductivity");
    REGISTER_FECORE_CLASS(FEThermoFluidPressureLoad, "thermo-fluid pressure load");

    febio.SetActiveModule(0);
}
