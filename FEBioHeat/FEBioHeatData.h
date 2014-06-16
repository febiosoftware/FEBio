#pragma once
#include "FECore/DataStore.h"
#include "FECore/NodeDataRecord.h"

//-----------------------------------------------------------------------------
class FENodeTemp : public FENodeLogData
{ 
public: 
	FENodeTemp(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};
