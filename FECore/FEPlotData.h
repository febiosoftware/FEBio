#pragma once

#include "FEMesh.h"
#include "FESurface.h"
#include "FECoreBase.h"
#include "Archive.h"
#include "FE_enum.h"

//-----------------------------------------------------------------------------
// forward declaration of model class
class FEModel;

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

	// The filter can be used to pass additional information to the plot field.
	// The interpretation of this filter is up to the derived class, but could
	// be used e.g. for disambiguation. There are currently two flavors of this
	// function: one that takes a const char* and one that takes an int. Derived
	// classes can overload both or just the one that makes sense for that plot field.
	// Note that these functions return false by default. This implies that trying
	// to specify a filter on a variable that doesn't support it will automatically cause 
	// an error.
	virtual bool SetFilter(const char* sz) { return false; }
	virtual bool SetFilter(int n) { return false;}

	Var_Type DataType() { return m_ntype; }
	Storage_Fmt StorageFormat() { return m_sfmt; }

	int VarSize(Var_Type t);

	void SetItemList(vector<int>& item) { m_item = item; }
    
    void SetDomainName(const char* szdom) { strcpy(m_szdom, szdom); }

protected:
	Var_Type		m_ntype;	//!< data type
	Storage_Fmt		m_sfmt;		//!< data storage format
	vector<int>		m_item;		//!< Data will only be stored for the item's in this list
    char			m_szdom[64];//!< Data will only be stored for the domain with this name
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
