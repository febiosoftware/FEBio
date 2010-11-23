#pragma once

#include "Archive.h"
#include "FEMesh.h"
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
	virtual void Save(FEM& fem, FILE* fp) = 0;
};

//-----------------------------------------------------------------------------
//! This is the base class for node data. Classes that wish to store data
//! associated with each node of the mesh, will use this base class.
class FENodeData : public FEPlotData
{
public:
	void Save(FEM& fem, FILE* fp);
	virtual void Save(FEMesh& m, FILE* fp) = 0;
};

//-----------------------------------------------------------------------------
//! This is the base class for element data. Classes that wish to store data
//! associated with each element of a domain, will use this base class.
class FEElementData : public FEPlotData
{
public:
	void Save(FEM& fem, FILE* fp);
	virtual void Save(FEDomain& D, FILE* fp) = 0;
};

//-----------------------------------------------------------------------------
//! This is the base class for facet data. Classes that wish to store data
//! associated with each facet of a surface, will use this base class.
class FEFaceData : public FEPlotData
{
public:
	void Save(FEM& fem, FILE* fp);
	virtual void Save(FESurface& S, FILE* fp) = 0;
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
	void Save(FEMesh& m, FILE* fp);
};

//-----------------------------------------------------------------------------
//! Nodal velocities
//!
class FEPlotNodeVelocity : public FENodeData
{
public:
	void Save(FEMesh& m, FILE* fp);
};

//-----------------------------------------------------------------------------
//! Nodal accelerations
//!
class FEPlotNodeAcceleration : public FENodeData
{
public:
	void Save(FEMesh& m, FILE* fp);
};

//-----------------------------------------------------------------------------
class FEPlotFluidPressure : public FENodeData
{
public:
	void Save(FEMesh& m, FILE* fp);
};

//=============================================================================
//                         E L E M E N T   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Element stresses
class FEPlotElementStress : public FEElementData
{
public:
	void Save(FEDomain& dom, FILE* fp);
};

//-----------------------------------------------------------------------------
//! Fluid flux
class FEPlotFluidFlux : public FEElementData
{
public:
	void Save(FEDomain& dom, FILE* fp);
};

//-----------------------------------------------------------------------------
//! Material fibers
class FEPlotFiberVector : public FEElementData
{
public:
	void Save(FEDomain& dom, FILE* fp);
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
	void Save(FESurface& surf, FILE* fp);
};

//-----------------------------------------------------------------------------
//! Contact traction
//!
class FEPlotContactTraction : public FEFaceData
{
public:
	void Save(FESurface& surf, FILE* fp);
};
