#pragma once

#include "FEMesh.h"
#include "FESurface.h"
#include "FECoreBase.h"
#include "Archive.h"

//-----------------------------------------------------------------------------
// forward declaration of model class
class FEModel;

// --- data types ---
enum Var_Type { PLT_FLOAT, PLT_VEC3F, PLT_MAT3FS, PLT_MAT3FD, PLT_TENS4FS, PLT_MAT3F };

// --- storage format ---
// FMT_NODE : one value stored for each node of a region
// FMT_ITEM : one value stored for each item (e.g. element) of a region
// FMT_MULT : one value for each node of each item of a region
enum Storage_Fmt { FMT_NODE, FMT_ITEM, FMT_MULT };

//-----------------------------------------------------------------------------
//! This is the base class for all classes that wish to store data to the 
//! plot file. However, classes will not use this class directly as their
//! base class. Instead they'll use one of the more specialized classes
//! defined below.
//!
class FEPlotData : public FECoreBase
{
public:
	FEPlotData(Var_Type t, Storage_Fmt s) : FECoreBase(FEPLOTDATA_ID) { m_ntype = t; m_sfmt = s; }
	virtual void Save(FEModel& fem, Archive& ar) = 0;

	Var_Type DataType() { return m_ntype; }
	Storage_Fmt StorageFormat() { return m_sfmt; }

	int VarSize(Var_Type t);

	void SetItemList(vector<int>& item) { m_item = item; }

protected:
	Var_Type		m_ntype;	//!< data type
	Storage_Fmt		m_sfmt;		//!< data storage format
	vector<int>		m_item;		//!< Data will only be stored for the item's in this list
};

//-----------------------------------------------------------------------------
//! This is the base class for node data. Classes that wish to store data
//! associated with each node of the mesh, will use this base class.
class FENodeData : public FEPlotData
{
public:
	FENodeData(Var_Type t, Storage_Fmt s) : FEPlotData(t, s) {}
	void Save(FEModel& fem, Archive& ar);
	virtual bool Save(FEMesh& m, vector<float>& a) = 0;
};

//-----------------------------------------------------------------------------
//! This is the base class for domain data. Classes that wish to store data
//! associated with each element or node of a domain, will use this base class.
class FEDomainData : public FEPlotData
{
public:
	FEDomainData(Var_Type t, Storage_Fmt s) : FEPlotData(t, s) {}
	void Save(FEModel& fem, Archive& ar);
	virtual bool Save(FEDomain& D, vector<float>& a) = 0;
};

//-----------------------------------------------------------------------------
//! This is the base class for surface data. Classes that wish to store data
//! associated with each node or facet of a surface, will use this base class.
class FESurfaceData : public FEPlotData
{
public:
	FESurfaceData(Var_Type t, Storage_Fmt s) : FEPlotData(t, s) {}
	void Save(FEModel& fem, Archive& ar);
	virtual bool Save(FESurface& S, vector<float>& a) = 0;
};
