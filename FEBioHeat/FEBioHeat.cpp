#include "FEBioHeat.h"
#include "FEIsotropicFourier.h"
#include "FEPlotNodeTemperature.h"
#include "FEPlotHeatFlux.h"
#include "FEHeatTransferAnalysis.h"
#include "FEHeatFlux.h"
#include "FEConvectiveHeatFlux.h"
#include "FEHeatSolver.h"
#include "FEHeatSolidDomain.h"
#include "FEHeatDomainFactory.h"
#include "FEHeatSource.h"
#include "FEBioHeatData.h"

namespace FEBioHeat {

//============================================================================
void InitModule()
{
	// Domain factory
	FECoreKernel& febio = FECoreKernel::GetInstance();
	febio.RegisterDomain(new FEHeatDomainFactory);

	// Analysis
	REGISTER_FECORE_CLASS(FEHeatTransferAnalysis, FEANALYSIS_ID, "heat transfer");

	// Solvers
	REGISTER_FECORE_CLASS(FEHeatSolver, FESOLVER_ID, "heat transfer");

	// Materials
	REGISTER_FECORE_CLASS(FEIsotropicFourier, FEMATERIAL_ID, "isotropic Fourier");

	// Body loads
	REGISTER_FECORE_CLASS(FEHeatSource, FEBODYLOAD_ID, "heat_source");

	// Surface loads
	REGISTER_FECORE_CLASS(FEHeatFlux          , FESURFACELOAD_ID, "heatflux"           );
	REGISTER_FECORE_CLASS(FEConvectiveHeatFlux, FESURFACELOAD_ID, "convective_heatflux");

	// Plot data fields
	REGISTER_FECORE_CLASS(FEPlotNodeTemperature	, FEPLOTDATA_ID, "temperature"	);
	REGISTER_FECORE_CLASS(FEPlotHeatFlux		, FEPLOTDATA_ID, "heat flux"	);

	// log data fields
	REGISTER_FECORE_CLASS(FENodeTemp, FENODELOGDATA_ID, "T");
}

}
