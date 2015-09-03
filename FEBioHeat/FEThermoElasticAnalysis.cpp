#include "FEThermoElasticAnalysis.h"
#include "FECore/FEModel.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"

//-----------------------------------------------------------------------------
FEThermoElasticAnalysis::FEThermoElasticAnalysis(FEModel* pfem) : FEAnalysis(pfem, FE_THERMO_ELASTIC)
{
}
