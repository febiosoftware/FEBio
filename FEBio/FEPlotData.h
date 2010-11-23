#pragma once

#include "DumpFile.h"
#include "FEMesh.h"
#include "Archive.h"
class FEM;

//-----------------------------------------------------------------------------
//! This is the base class for all classes that wish to store data to the 
//! plot file. However, classes will not use this class directly as their
//! base class. Instead they'll use one of the more specialized classes
//! defined below.
//!
class FEPlotData
{
public:
	virtual void Save(FEM& fem, Archive& ar) = 0;
};

//-----------------------------------------------------------------------------
//! This is the base class for node data. Classes that wish to store data
//! associated with each node of the mesh, will use this base class.
class FENodeData : public FEPlotData
{
public:
	void Save(FEM& fem, Archive& ar);
	virtual void Save(FEMesh& m, Archive& ar) = 0;
};

//-----------------------------------------------------------------------------
//! This is the base class for element data. Classes that wish to store data
//! associated with each element of a domain, will use this base class.
class FEElementData : public FEPlotData
{
public:
	void Save(FEM& fem, Archive& ar);
	virtual void Save(FEDomain& D, Archive& ar) = 0;
};

//-----------------------------------------------------------------------------
//! This is the base class for facet data. Classes that wish to store data
//! associated with each facet of a surface, will use this base class.
class FEFaceData : public FEPlotData
{
public:
	void Save(FEM& fem, Archive& ar);
	virtual void Save(FESurface& S, Archive& ar) = 0;
};

//=============================================================================
//                            N O D E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Nodal displacements
//!
class FEPlotNodeDisplacement : public FENodeData
{
public:
	void Save(FEMesh& m, Archive& ar);
};

//-----------------------------------------------------------------------------
//! Nodal velocities
//!
class FEPlotNodeVelocity : public FENodeData
{
public:
	void Save(FEMesh& m, Archive& ar);
};

//-----------------------------------------------------------------------------
//! Nodal accelerations
//!
class FEPlotNodeAcceleration : public FENodeData
{
public:
	void Save(FEMesh& m, Archive& ar);
};

//-----------------------------------------------------------------------------
class FEPlotFluidPressure : public FENodeData
{
public:
	void Save(FEMesh& m, Archive& ar);
};

//=============================================================================
//                         E L E M E N T   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Element stresses
class FEPlotElementStress : public FEElementData
{
public:
	void Save(FEDomain& dom, Archive& ar);
};

//-----------------------------------------------------------------------------
//! Fluid flux
class FEPlotFluidFlux : public FEElementData
{
public:
	void Save(FEDomain& dom, Archive& ar);
};

//-----------------------------------------------------------------------------
//! Material fibers
class FEPlotFiberVector : public FEElementData
{
public:
	void Save(FEDomain& dom, Archive& ar);
};

//=============================================================================
//                         S U R F A C E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Contact gap
//!
class FEPlotContactGap : public FEFaceData
{
public:
	void Save(FESurface& surf, Archive& ar);
};

//-----------------------------------------------------------------------------
//! Contact traction
//!
class FEPlotContactTraction : public FEFaceData
{
public:
	void Save(FESurface& surf, Archive& ar);
};
