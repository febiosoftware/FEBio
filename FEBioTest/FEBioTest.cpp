#include "FEBioTest.h"
#include <FECore/FECoreKernel.h>
#include "FEBioDiagnostic.h"
#include "FETangentDiagnostic.h"
#include "FERestartDiagnostics.h"

namespace FEBioTest
{

void InitModule()
{
	REGISTER_FECORE_CLASS(FEBioDiagnostic, "diagnose");
	REGISTER_FECORE_CLASS(FERestartDiagnostic, "restart_test");
}
}
