#include "FEHeatTransferAnalysis.h"
#include "FECore/FEModel.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"

//-----------------------------------------------------------------------------
FEHeatTransferAnalysis::FEHeatTransferAnalysis(FEModel* pfem) : FEAnalysis(pfem, FE_HEAT)
{
}
