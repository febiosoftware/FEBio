#include "stdafx.h"
#include "FEBioOpt.h"
#include "FEOptimize.h"
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
//! Initialization of the FEBioOpt module. This function registers all the classes
//! in this module with the FEBio framework.
void FEBioOpt::InitModule()
{
REGISTER_FECORE_CLASS(FEOptimize, "optimize");
}
