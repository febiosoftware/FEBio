#pragma once
#include "FECore/FEPlotData.h"

//=============================================================================
//                            N O D E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Nodal displacements
//!
class FEPlotNodeDisplacement : public FENodeData
{
public:
	FEPlotNodeDisplacement(FEModel* pfem) : FENodeData(PLT_VEC3F, FMT_NODE){}
	bool Save(FEMesh& m, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Nodal velocities
//!
class FEPlotNodeVelocity : public FENodeData
{
public:
	FEPlotNodeVelocity(FEModel* pfem) : FENodeData(PLT_VEC3F, FMT_NODE){}
	bool Save(FEMesh& m, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Nodal accelerations
//!
class FEPlotNodeAcceleration : public FENodeData
{
public:
	FEPlotNodeAcceleration(FEModel* pfem) : FENodeData(PLT_VEC3F, FMT_NODE){}
	bool Save(FEMesh& m, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Nodal temperatures
class FEPlotNodeTemperature : public FENodeData
{
public:
	FEPlotNodeTemperature(FEModel* pfem) : FENodeData(PLT_FLOAT, FMT_NODE){}
	bool Save(FEMesh& m, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Nodal reaction forces
class FEPlotNodeReactionForces : public FENodeData
{
public:
	FEPlotNodeReactionForces(FEModel* pfem) : FENodeData(PLT_VEC3F, FMT_NODE){}
	bool Save(FEMesh& m, vector<float>& a);
};
