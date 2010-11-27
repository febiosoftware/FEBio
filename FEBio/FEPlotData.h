#pragma once

#include "DumpFile.h"
#include "FEMesh.h"
#include "Archive.h"
class FEM;

// --- data types ---
enum Var_Type { FLOAT, VEC3F, MAT3FS };

// --- storage format ---
// ITEM_DATA : one data stored for each item
// NODE_DATA : data stored for each node of each item
enum Storage_Fmt { ITEM_DATA, NODE_DATA };

//-----------------------------------------------------------------------------
//! This is the base class for all classes that wish to store data to the 
//! plot file. However, classes will not use this class directly as their
//! base class. Instead they'll use one of the more specialized classes
//! defined below.
//!
class FEPlotData
{
public:
	FEPlotData(Var_Type t, Storage_Fmt s) { m_ntype = t; m_sfmt = s; }
	virtual void Save(FEM& fem, Archive& ar) = 0;

	Var_Type DataType() { return m_ntype; }
	Storage_Fmt StorageFormat() { return m_sfmt; }

	int VarSize(Var_Type t);

protected:
	Var_Type		m_ntype;
	Storage_Fmt		m_sfmt;
};

//-----------------------------------------------------------------------------
//! This is the base class for node data. Classes that wish to store data
//! associated with each node of the mesh, will use this base class.
class FENodeData : public FEPlotData
{
public:
	FENodeData(Var_Type t, Storage_Fmt s) : FEPlotData(t, s) {}
	void Save(FEM& fem, Archive& ar);
	virtual void Save(FEMesh& m, vector<float>& a) = 0;
};

//-----------------------------------------------------------------------------
//! This is the base class for element data. Classes that wish to store data
//! associated with each element of a domain, will use this base class.
class FEElementData : public FEPlotData
{
public:
	FEElementData(Var_Type t, Storage_Fmt s) : FEPlotData(t, s) {}
	void Save(FEM& fem, Archive& ar);
	virtual void Save(FEDomain& D, vector<float>& a) = 0;
};

//-----------------------------------------------------------------------------
//! This is the base class for facet data. Classes that wish to store data
//! associated with each facet of a surface, will use this base class.
class FEFaceData : public FEPlotData
{
public:
	FEFaceData(Var_Type t, Storage_Fmt s) : FEPlotData(t, s) {}
	void Save(FEM& fem, Archive& ar);
	virtual void Save(FESurface& S, vector<float>& a) = 0;
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
	FEPlotNodeDisplacement() : FENodeData(VEC3F, ITEM_DATA){}
	void Save(FEMesh& m, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Nodal velocities
//!
class FEPlotNodeVelocity : public FENodeData
{
public:
	FEPlotNodeVelocity() : FENodeData(VEC3F, ITEM_DATA){}
	void Save(FEMesh& m, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Nodal accelerations
//!
class FEPlotNodeAcceleration : public FENodeData
{
public:
	FEPlotNodeAcceleration() : FENodeData(VEC3F, ITEM_DATA){}
	void Save(FEMesh& m, vector<float>& a);
};

//-----------------------------------------------------------------------------
class FEPlotFluidPressure : public FENodeData
{
public:
	FEPlotFluidPressure() : FENodeData(FLOAT, ITEM_DATA){}
	void Save(FEMesh& m, vector<float>& a);
};

//=============================================================================
//                         E L E M E N T   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Element stresses
class FEPlotElementStress : public FEElementData
{
public:
	FEPlotElementStress() : FEElementData(MAT3FS, ITEM_DATA){}
	void Save(FEDomain& dom, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Fluid flux
class FEPlotFluidFlux : public FEElementData
{
public:
	FEPlotFluidFlux() : FEElementData(VEC3F, ITEM_DATA){}
	void Save(FEDomain& dom, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Material fibers
class FEPlotFiberVector : public FEElementData
{
public:
	FEPlotFiberVector() : FEElementData(VEC3F, ITEM_DATA){}
	void Save(FEDomain& dom, vector<float>& a);
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
	FEPlotContactGap() : FEFaceData(FLOAT, ITEM_DATA){}
	void Save(FESurface& surf, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Contact traction
//!
class FEPlotContactTraction : public FEFaceData
{
public:
	FEPlotContactTraction() : FEFaceData(VEC3F, ITEM_DATA){}
	void Save(FESurface& surf, vector<float>& a);
};
