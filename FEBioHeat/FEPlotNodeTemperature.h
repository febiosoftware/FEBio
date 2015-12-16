#pragma once
#include "FECore/FEPlotData.h"

//-----------------------------------------------------------------------------
//! Nodal temperatures
class FEPlotNodeTemperature : public FENodeData
{
public:
	FEPlotNodeTemperature(FEModel* pfem) : FENodeData(PLT_FLOAT, FMT_NODE) {}
	bool Save(FEMesh& m, FEDataStream& a);
};
