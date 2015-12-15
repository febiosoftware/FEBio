#pragma once
#include "FECore/FEPlotData.h"

//-----------------------------------------------------------------------------
//! Nodal temperatures
class FEPlotNodeTemperature : public FENodeData
{
public:
	FEPlotNodeTemperature(FEModel* pfem) : FENodeData(PLT_FLOAT, FMT_NODE), m_pfem(pfem) {}
	bool Save(FEMesh& m, FEDataStream& a);
private:
	FEModel*	m_pfem;
};
