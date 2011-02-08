#pragma once
#include "FEPlotData.h"

//=============================================================================
//                            N O D E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Nodal displacements
//!
class FEPlotNodeDisplacement : public FENodeData
{
public:
	FEPlotNodeDisplacement() : FENodeData(VEC3F, FMT_NODE){}
	bool Save(FEMesh& m, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Nodal velocities
//!
class FEPlotNodeVelocity : public FENodeData
{
public:
	FEPlotNodeVelocity() : FENodeData(VEC3F, FMT_NODE){}
	bool Save(FEMesh& m, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Nodal accelerations
//!
class FEPlotNodeAcceleration : public FENodeData
{
public:
	FEPlotNodeAcceleration() : FENodeData(VEC3F, FMT_NODE){}
	bool Save(FEMesh& m, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Nodal temperatures
class FEPlotNodeTemperature : public FENodeData
{
public:
	FEPlotNodeTemperature() : FENodeData(FLOAT, FMT_NODE){}
	bool Save(FEMesh& m, vector<float>& a);
};
