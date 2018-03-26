#include "FEBioFSI.h"
#include <FECore/FECoreKernel.h>
#include "FEFluidFSISolver.h"
#include "FEFluidFSI.h"
#include "FEFluidFSIDomain3D.h"
#include "FEFluidFSITraction.h"
#include "FEFluidFSIDomainFactory.h"

void FEBioFSI::InitModule()
{
	FECoreKernel& febio = FECoreKernel::GetInstance();

	// register domain
	febio.RegisterDomain(new FEFluidFSIDomainFactory);

	// define the fsi module
	febio.CreateModule("fluid-FSI");
	febio.SetModuleDependency("fluid");

	REGISTER_FECORE_CLASS(FEFluidFSISolver, FESOLVER_ID, "fluid-FSI");

	REGISTER_FECORE_CLASS(FEFluidFSI, FEMATERIAL_ID, "fluid-FSI");

	REGISTER_FECORE_CLASS(FEFluidFSIDomain3D, FEDOMAIN_ID, "fluid-FSI-3D");

	REGISTER_FECORE_CLASS(FEFluidFSITraction, FESURFACELOAD_ID, "fluid-FSI traction");

	febio.SetActiveModule(0);
}
