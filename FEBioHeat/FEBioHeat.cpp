#include "FEBioHeat.h"
#include "FEIsotropicFourier.h"
#include "FEPlotNodeTemperature.h"
#include "FEPlotHeatFlux.h"
#include "FEHeatTransferAnalysis.h"

namespace FEBioHeat {

void InitModule()
{
// Analysis
REGISTER_FEBIO_CLASS(FEHeatTransferAnalysis, FEAnalysis, "heat transfer");

// Materials
REGISTER_MATERIAL(FEIsotropicFourier, "isotropic Fourier");

// Plot data fields
REGISTER_FEBIO_CLASS(FEPlotNodeTemperature	, FEPlotData, "temperature"	);
REGISTER_FEBIO_CLASS(FEPlotHeatFlux			, FEPlotData, "heat flux"	);

}

}
