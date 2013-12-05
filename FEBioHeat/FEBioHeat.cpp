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

namespace FEBioHeat {

//============================================================================
void InitModule()
{
	// Domain factory
	FEBioKernel& febio = FEBioKernel::GetInstance();
	febio.RegisterDomain(new FEHeatDomainFactory);

	// Analysis
	REGISTER_FEBIO_CLASS(FEHeatTransferAnalysis, FEANALYSIS_ID, "heat transfer");

	// Solvers
	REGISTER_FEBIO_CLASS(FEHeatSolver, FESOLVER_ID, "heat transfer");

	// Materials
	REGISTER_FEBIO_CLASS(FEIsotropicFourier, FEMATERIAL_ID, "isotropic Fourier");

	// Body loads
	REGISTER_FEBIO_CLASS(FEHeatSource, FEBODYLOAD_ID, "heat_source");

	// Surface loads
	REGISTER_FEBIO_CLASS(FEHeatFlux          , FESURFACELOAD_ID, "heatflux"           );
	REGISTER_FEBIO_CLASS(FEConvectiveHeatFlux, FESURFACELOAD_ID, "convective_heatflux");

	// Plot data fields
	REGISTER_FEBIO_CLASS(FEPlotNodeTemperature	, FEPLOTDATA_ID, "temperature"	);
	REGISTER_FEBIO_CLASS(FEPlotHeatFlux			, FEPLOTDATA_ID, "heat flux"	);
}

}
