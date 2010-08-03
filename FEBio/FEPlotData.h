#pragma once

#include "Archive.h"
class FEM;

//-----------------------------------------------------------------------------
//! This is the base class for all classes that wisht to store data in the 
//! plot file. 
//!
class FEPlotData
{
public:
	virtual void Save(FEM& fem, FILE* fp) = 0;
};

//-----------------------------------------------------------------------------
//! Nodal displacements
//!
class FEPlotNodeDisplacement : public FEPlotData
{
public:
	void Save(FEM& fem, FILE* fp);
};

//-----------------------------------------------------------------------------
//! Nodal velocities
//!
class FEPlotNodeVelocity : public FEPlotData
{
public:
	void Save(FEM& fem, FILE* fp);
};

//-----------------------------------------------------------------------------
//! Nodal accelerations
//!
class FEPlotNodeAcceleration : public FEPlotData
{
public:
	void Save(FEM& fem, FILE* fp);
};

//-----------------------------------------------------------------------------
//! Element stresses
//!
class FEPlotElementStress : public FEPlotData
{
public:
	void Save(FEM& fem, FILE* fp);
};

//-----------------------------------------------------------------------------
//! Contact gap
//!
class FEPlotContactGap : public FEPlotData
{
public:
	void Save(FEM& fem, FILE* fp);
};

//-----------------------------------------------------------------------------
//! Contact traction
//!
class FEPlotContactTraction : public FEPlotData
{
public:
	void Save(FEM& fem, FILE* fp);
};

//-----------------------------------------------------------------------------
class FEPlotFluidPressure : public FEPlotData
{
public:
	void Save(FEM& fem, FILE* fp);
};

//-----------------------------------------------------------------------------
class FEPlotFluidFlux : public FEPlotData
{
public:
	void Save(FEM& fem, FILE* fp);
};

//-----------------------------------------------------------------------------
class FEPlotFiberVector : public FEPlotData
{
public:
	void Save(FEM& fem, FILE* fp);
};
