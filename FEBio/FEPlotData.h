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
	virtual void Save(FEM& fem, Archive& ar) = 0;
};

//-----------------------------------------------------------------------------
//! Nodal displacements
//!
class FEPlotNodeDisplacement : public FEPlotData
{
public:
	void Save(FEM& fem, Archive& ar);
};

//-----------------------------------------------------------------------------
//! Nodal velocities
//!
class FEPlotNodeVelocity : public FEPlotData
{
public:
	void Save(FEM& fem, Archive& ar);
};

//-----------------------------------------------------------------------------
//! Nodal accelerations
//!
class FEPlotNodeAcceleration : public FEPlotData
{
public:
	void Save(FEM& fem, Archive& ar);
};

//-----------------------------------------------------------------------------
//! Element stresses
//!
class FEPlotElementStress : public FEPlotData
{
public:
	void Save(FEM& fem, Archive& ar);
};

//-----------------------------------------------------------------------------
//! Contact gap
//!
class FEPlotContactGap : public FEPlotData
{
public:
	void Save(FEM& fem, Archive& ar);
};

//-----------------------------------------------------------------------------
//! Contact traction
//!
class FEPlotContactTraction : public FEPlotData
{
public:
	void Save(FEM& fem, Archive& ar);
};

//-----------------------------------------------------------------------------
class FEPlotFluidPressure : public FEPlotData
{
public:
	void Save(FEM& fem, Archive& ar);
};

//-----------------------------------------------------------------------------
class FEPlotFluidFlux : public FEPlotData
{
public:
	void Save(FEM& fem, Archive& ar);
};
