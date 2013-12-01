#include "FEBioHeat.h"
#include "FEIsotropicFourier.h"
#include "FEPlotNodeTemperature.h"
#include "FEPlotHeatFlux.h"
#include "FEHeatTransferAnalysis.h"
#include "FEHeatFlux.h"
#include "FEConvectiveHeatFlux.h"

namespace FEBioHeat {

void InitModule()
{
// Analysis
REGISTER_FEBIO_CLASS(FEHeatTransferAnalysis, FEAnalysis, "heat transfer");

// Materials
REGISTER_MATERIAL(FEIsotropicFourier, "isotropic Fourier");

// Surface loads
REGISTER_FEBIO_CLASS(FEHeatFlux          , FESurfaceLoad, "heatflux"           );
REGISTER_FEBIO_CLASS(FEConvectiveHeatFlux, FESurfaceLoad, "convective_heatflux");

// Plot data fields
REGISTER_FEBIO_CLASS(FEPlotNodeTemperature	, FEPlotData, "temperature"	);
REGISTER_FEBIO_CLASS(FEPlotHeatFlux			, FEPlotData, "heat flux"	);
}

}
