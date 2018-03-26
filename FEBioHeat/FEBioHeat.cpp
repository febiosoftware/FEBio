#include "FEBioHeat.h"
#include "FEIsotropicFourier.h"
#include "FEPlotHeatFlux.h"
#include "FEHeatFlux.h"
#include "FEConvectiveHeatFlux.h"
#include "FEHeatSolver.h"
#include "FEHeatSolidDomain.h"
#include "FEHeatDomainFactory.h"
#include "FEHeatSource.h"
#include "FEThermoElasticSolver.h"
#include "FEThermoElasticMaterial.h"
#include "FEThermoNeoHookean.h"
#include "FEThermalConductivity.h"
#include "FEHeatSolidDomain.h"
#include "FEThermoElasticSolidDomain.h"

namespace FEBioHeat {

//============================================================================
void InitModule()
{
	// Domain factory
	FECoreKernel& febio = FECoreKernel::GetInstance();
	febio.RegisterDomain(new FEHeatDomainFactory);

	// create the heat-transfer module
	febio.CreateModule("heat");

	// Solvers
	REGISTER_FECORE_CLASS(FEHeatSolver         , FESOLVER_ID, "heat" );
//	REGISTER_FECORE_CLASS(FEThermoElasticSolver, FESOLVER_ID, "thermo-elastic");

	// Materials
	REGISTER_FECORE_CLASS(FEIsotropicFourier     , FEMATERIAL_ID, "isotropic Fourier");
	REGISTER_FECORE_CLASS(FEThermoElasticMaterial, FEMATERIAL_ID, "thermo-elastic"   );
	REGISTER_FECORE_CLASS(FEThermoNeoHookean     , FEMATERIAL_ID, "thermo-neo-Hookean");
	REGISTER_FECORE_CLASS(FEConstReferenceThermalConductivity, FEMATERIAL_ID, "cond-ref-iso");
	REGISTER_FECORE_CLASS(FEConstThermalConductivity         , FEMATERIAL_ID, "cond-const"  );

	// Domain classes
	REGISTER_FECORE_CLASS(FEHeatSolidDomain         , FEDOMAIN_ID, "heat-solid"          );
	REGISTER_FECORE_CLASS(FEThermoElasticSolidDomain, FEDOMAIN_ID, "thermo-elastic-solid");

	// Body loads
	REGISTER_FECORE_CLASS(FEHeatSource, FEBODYLOAD_ID, "heat_source");

	// Surface loads
	REGISTER_FECORE_CLASS(FEHeatFlux          , FESURFACELOAD_ID, "heatflux"           );
	REGISTER_FECORE_CLASS(FEConvectiveHeatFlux, FESURFACELOAD_ID, "convective_heatflux");

	// Plot data fields
	REGISTER_FECORE_CLASS(FEPlotHeatFlux		, FEPLOTDATA_ID, "heat flux"	);

	febio.SetActiveModule(0);
}

}
