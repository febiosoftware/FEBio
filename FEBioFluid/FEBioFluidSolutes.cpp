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
#include "FEFluidModule.h"
#include "FEBioFluidSolutes.h"
#include "FEFluidSolutesSolver.h"
#include "FEFluidSolutes.h"
#include "FEFluidSolutesDomain3D.h"
#include "FEFluidSolutesDomainFactory.h"
#include "FESoluteBackflowStabilization.h"
#include "FEInitialFluidSolutesPressure.h"
#include "FEFluidSolutesFlux.h"
#include "FEFluidSolutesNaturalFlux.h"
#include "FEFluidSolutesPressure.h"
#include "FEFluidSolutesPressureBC.h"
#include "FEFluidSolutesResistanceBC.h"
#include "FEFluidSolutesRCRBC.h"
#include "FESoluteConvectiveFlow.h"
#include "FEFluidSolutesPressureLC.h"
#include "FEFluidSolutesGradientLC.h"
#include "FEFluidSolutesDomainFactory.h"
#include "FESolutesSolver.h"
#include "FESolutesMaterial.h"
#include "FESolutesDomain.h"
#include "FESolutesDomainFactory.h"
#include "FEBioFluidPlot.h"
#include <FEBioMix/FESoluteFlux.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FEModelUpdate.h>
#include <FECore/FETimeStepController.h>
#include "FEFluidSolutesAnalysis.h"

//-----------------------------------------------------------------------------
const char* FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_SOLUTES_VARIABLE var)
{
    switch (var)
    {
        case DISPLACEMENT                : return "displacement"               ; break;
        case RELATIVE_FLUID_VELOCITY     : return "relative fluid velocity"    ; break;
        case RELATIVE_FLUID_ACCELERATION : return "relative fluid acceleration"; break;
        case FLUID_DILATATION            : return "fluid dilatation"             ; break;
        case FLUID_DILATATION_TDERIV     : return "fluid dilatation tderiv"      ; break;
        case FLUID_CONCENTRATION         : return "concentration"              ; break;
        case FLUID_CONCENTRATION_TDERIV  : return "concentration tderiv"       ; break;
    }
    assert(false);
    return nullptr;
}

void FEBioFluidSolutes::InitModule()
{
    FECoreKernel& febio = FECoreKernel::GetInstance();
    
    // register domain
    febio.RegisterDomain(new FEFluidSolutesDomainFactory);
    
    // define the fsi module
    febio.CreateModule(new FEFluidSolutesModule, "fluid-solutes",
                       "{"
                       "   \"title\" : \"Fluid-Solutes\","
                       "   \"info\"  : \"Fluid analysis with solute transport and reactive processes.\""
                       "}");
	febio.AddModuleDependency("fluid");
    febio.AddModuleDependency("multiphasic"); // also pulls in solid, biphasic, solutes
    
    //-----------------------------------------------------------------------------
    // analyis classes (default type must match module name!)
    REGISTER_FECORE_CLASS(FEFluidSolutesAnalysis, "fluid-solutes");

	// monolithic fluid-solutes solver
    REGISTER_FECORE_CLASS(FEFluidSolutesSolver, "fluid-solutes");
    REGISTER_FECORE_CLASS(FEFluidSolutes, "fluid-solutes");
    REGISTER_FECORE_CLASS(FEFluidSolutesDomain3D, "fluid-solutes-3D");
    
    // loads
    REGISTER_FECORE_CLASS(FEFluidSolutesFlux           , "solute flux"                  );
    REGISTER_FECORE_CLASS(FESoluteBackflowStabilization, "solute backflow stabilization");
    REGISTER_FECORE_CLASS(FEFluidSolutesNaturalFlux    , "solute natural flux");
    REGISTER_FECORE_CLASS(FEFluidSolutesPressure       , "fluid pressure"               , 0x0300); // deprecated, use BC version
    REGISTER_FECORE_CLASS(FESoluteConvectiveFlow       , "solute convective flow"       , FECORE_EXPERIMENTAL);

    // bcs
    REGISTER_FECORE_CLASS(FEFluidSolutesPressureBC     , "fluid pressure"  );
    REGISTER_FECORE_CLASS(FEFluidSolutesResistanceBC   , "fluid resistance");
    REGISTER_FECORE_CLASS(FEFluidSolutesRCRBC          , "fluid RCR"       );

    // ics
    REGISTER_FECORE_CLASS(FEInitialFluidSolutesPressure, "initial fluid pressure");

    // constraints
    REGISTER_FECORE_CLASS(FEFluidSolutesPressureLC     , "fluid pressure constraint", FECORE_EXPERIMENTAL);
    REGISTER_FECORE_CLASS(FEFluidSolutesGradientLC     , "zero concentration gradient", FECORE_EXPERIMENTAL);
    
    //-----------------------------------------------------------------------------
    // classes derived from FEPlotData
    REGISTER_FECORE_CLASS(FEPlotFluidRelativePecletNumber, "solute relative Peclet number");

	// solutes solver classes
	febio.RegisterDomain(new FESolutesDomainFactory);
	REGISTER_FECORE_CLASS(FESolutesSolver, "solutes");
	REGISTER_FECORE_CLASS(FESolutesMaterial, "solutes");
	REGISTER_FECORE_CLASS(FESolutesDomain, "solutes-3D");

    //-----------------------------------------------------------------------------
    // Reset solver parameters to preferred default settings
    febio.OnCreateEvent(CallWhenCreating<FENewtonStrategy>([](FENewtonStrategy* pc) {
        pc->m_maxups = 20;
    }));
    
    febio.OnCreateEvent(CallWhenCreating<FETimeStepController>([](FETimeStepController* pc) {
        pc->m_iteopt = 100;
    }));
    
    febio.OnCreateEvent(CallWhenCreating<FEFluidSolutesAnalysis>([](FEFluidSolutesAnalysis* pc) {
        pc->m_nanalysis = FEFluidSolutesAnalysis::DYNAMIC;
    }));
    
    febio.OnCreateEvent(CallWhenCreating<FENewtonSolver>([](FENewtonSolver* pc) {
        pc->m_maxref = 5;
        pc->m_Rmax = 1.0e+20;
        // turn off reform on each time step and diverge reform
        pc->m_breformtimestep = true;
        pc->m_bdivreform = false;
    }));
    
    febio.SetActiveModule(0);
}
