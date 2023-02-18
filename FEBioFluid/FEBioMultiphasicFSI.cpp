/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
#include "FEBioMultiphasicFSI.h"
#include <FECore/FECoreKernel.h>
#include "FEMultiphasicFSISolver.h"
#include "FEMultiphasicFSI.h"
#include "FEMultiphasicFSIDomain3D.h"
#include "FEBiphasicFSITraction.h"
#include "FEMultiphasicFSIDomainFactory.h"
#include "FEBackFlowFSIStabilization.h"
#include "FETangentialFlowFSIStabilization.h"
#include "FEMultiphasicFSISoluteFlux.h"
#include "FEMultiphasicFSIPressure.h"
#include "FEMultiphasicFSIPressureBC.h"
#include "FESoluteConvectiveFlow.h"
#include "FEMultiphasicFSISoluteBackflowStabilization.h"
#include "FEFluidModule.h"
#include "FEMultiphasicFSIAnalysis.h"
#include <FECore/FEModelUpdate.h>
#include <FECore/FETimeStepController.h>

//-----------------------------------------------------------------------------
const char* FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::MULTIPHASIC_FSI_VARIABLE var)
{
    switch (var)
    {
        case DISPLACEMENT                : return "displacement"               ; break;
        case VELOCITY                    : return "velocity"                   ; break;
        case SHELL_ROTATION              : return "shell rotation"             ; break;
        case SHELL_DISPLACEMENT          : return "shell displacement"         ; break;
        case SHELL_VELOCITY              : return "shell velocity"             ; break;
        case SHELL_ACCELERATION          : return "shell acceleration"         ; break;
        case RIGID_ROTATION              : return "rigid rotation"             ; break;
        case RELATIVE_FLUID_VELOCITY     : return "relative fluid velocity"    ; break;
        case RELATIVE_FLUID_ACCELERATION : return "relative fluid acceleration"; break;
        case FLUID_VELOCITY              : return "fluid velocity"             ; break;
        case FLUID_ACCELERATION          : return "fluid acceleration"         ; break;
        case FLUID_DILATATION            : return "fluid dilatation"           ; break;
        case FLUID_DILATATION_TDERIV     : return "fluid dilatation tderiv"    ; break;
        case FLUID_CONCENTRATION         : return "concentration"              ; break;
        case FLUID_CONCENTRATION_TDERIV  : return "concentration tderiv"       ; break;
    }
    assert(false);
    return nullptr;
}

void FEBioMultiphasicFSI::InitModule()
{
    FECoreKernel& febio = FECoreKernel::GetInstance();
    
    // register domain
    febio.RegisterDomain(new FEMultiphasicFSIDomainFactory);
    
    // define the fsi module
    febio.CreateModule(new FEMultiphasicFSIModule, "multiphasic-FSI");
    febio.AddModuleDependency("fluid");
    febio.AddModuleDependency("multiphasic");    // also pulls in solid, biphasic, solutes
    
    //-----------------------------------------------------------------------------
    // analyis classes (default type must match module name!)
    REGISTER_FECORE_CLASS(FEMultiphasicFSIAnalysis, "multiphasic-FSI");

    //-----------------------------------------------------------------------------
    REGISTER_FECORE_CLASS(FEMultiphasicFSISolver, "multiphasic-FSI");
    
    REGISTER_FECORE_CLASS(FEMultiphasicFSI, "multiphasic-FSI");
    
    REGISTER_FECORE_CLASS(FEMultiphasicFSIDomain3D, "multiphasic-FSI-3D");
    
    REGISTER_FECORE_CLASS(FEBiphasicFSITraction, "multiphasic-FSI traction");
    
    REGISTER_FECORE_CLASS(FEBackFlowFSIStabilization, "fluid backflow stabilization");
    
    REGISTER_FECORE_CLASS(FETangentialFlowFSIStabilization, "fluid tangential stabilization");
    
    // loads
    REGISTER_FECORE_CLASS(FEMultiphasicFSISoluteFlux, "solute flux");
    REGISTER_FECORE_CLASS(FEMultiphasicFSISoluteBackflowStabilization, "solute backflow stabilization");
    REGISTER_FECORE_CLASS(FEMultiphasicFSIPressure, "fluid pressure", 0x0300); // deprecated, use BC version

    // bcs
    REGISTER_FECORE_CLASS(FEMultiphasicFSIPressureBC, "fluid pressure");

    //-----------------------------------------------------------------------------
    // Reset solver parameters to preferred default settings
    febio.OnCreateEvent(CallWhenCreating<FENewtonStrategy>([](FENewtonStrategy* pc) {
        pc->m_maxups = 50;
    }));
    
    febio.OnCreateEvent(CallWhenCreating<FETimeStepController>([](FETimeStepController* pc) {
        pc->m_iteopt = 50;
    }));
    
    febio.OnCreateEvent(CallWhenCreating<FEMultiphasicFSIAnalysis>([](FEMultiphasicFSIAnalysis* pc) {
        pc->m_nanalysis = FEMultiphasicFSIAnalysis::DYNAMIC;
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
