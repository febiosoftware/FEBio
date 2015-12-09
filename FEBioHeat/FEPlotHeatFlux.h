#include "FECore/FEPlotData.h"

//-----------------------------------------------------------------------------
//! Class the outputs the heat flux
class FEPlotHeatFlux : public FEDomainData
{
public:
	FEPlotHeatFlux(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};
