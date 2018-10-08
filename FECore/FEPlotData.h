#pragma once

#include "FECoreBase.h"
#include "FE_enum.h"
#include "FEDataStream.h"
#include <functional>

//-----------------------------------------------------------------------------
// Region types
enum Region_Type {
	FE_REGION_NODE,
	FE_REGION_DOMAIN,
	FE_REGION_SURFACE
};

//-----------------------------------------------------------------------------
// forward declarations
class FEModel;
class FENode;
class FEMesh;
class FESurface;
class FEDomain;
class FESolidDomain;
class FEMaterialPoint;
class FENodeSet;
class FEElement;

//-----------------------------------------------------------------------------
//! This is the base class for all classes that wish to store data to the 
//! plot file. However, classes will not use this class directly as their
//! base class. Instead they'll use one of the more specialized classes
//! defined below.
//!
class FECORE_API FEPlotData : public FECoreBase
{
public:
	FEPlotData(FEModel* fem);
	FEPlotData(FEModel* fem, Region_Type R, Var_Type t, Storage_Fmt s);

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

	Region_Type RegionType() { return m_nregion; }
	Var_Type DataType() { return m_ntype; }
	Storage_Fmt StorageFormat() { return m_sfmt; }

	int VarSize(Var_Type t);

	void SetItemList(vector<int>& item) { m_item = item; }

	vector<int> GetItemList() { return m_item; }
    
    void SetDomainName(const char* szdom);
	const char* GetDomainName() { return m_szdom;  }

	FEModel* GetFEModel() { return m_fem; }

protected:
	void SetRegionType(Region_Type rt) { m_nregion = rt; }
	void SetVarType(Var_Type vt) { m_ntype = vt; }
	void SetStorageFormat(Storage_Fmt sf) { m_sfmt = sf; }

public: // override one of these functions depending on the Region_Type
	virtual bool Save(FEMesh&    m, FEDataStream& a) { return false; }		// for FE_REGION_NODE
	virtual bool Save(FEDomain&  D, FEDataStream& a) { return false; }		// for FE_REGION_DOMAIN
	virtual bool Save(FESurface& S, FEDataStream& a) { return false; }		// for FE_REGION_SURFACE

public: // used by array variables
	void SetArraySize(int n);
	int GetArraysize() const;

	void SetArrayNames(vector<string>& s) { m_arrayNames = s; }
	vector<string>& GetArrayNames() { return m_arrayNames; }

private:
	Region_Type		m_nregion;		//!< region type
	Var_Type		m_ntype;		//!< data type
	Storage_Fmt		m_sfmt;			//!< data storage format
	vector<int>		m_item;			//!< Data will only be stored for the item's in this list
    char			m_szdom[64];	//!< Data will only be stored for the domain with this name
	FEModel*		m_fem;

	int				m_arraySize;	//!< size of arrays (used by arrays)
	vector<string>	m_arrayNames;	//!< optional names of array components (used by arrays)
};

//-----------------------------------------------------------------------------
//! This is the base class for node data. Classes that wish to store data
//! associated with each node of the mesh, will use this base class.
class FECORE_API FEPlotNodeData : public FEPlotData
{
public:
	FEPlotNodeData(FEModel* fem, Var_Type t, Storage_Fmt s) : FEPlotData(fem, FE_REGION_NODE, t, s) {}
};

//-----------------------------------------------------------------------------
//! This is the base class for domain data. Classes that wish to store data
//! associated with each element or node of a domain, will use this base class.
class FECORE_API FEPlotDomainData : public FEPlotData
{
public:
	FEPlotDomainData(FEModel* fem, Var_Type t, Storage_Fmt s) : FEPlotData(fem, FE_REGION_DOMAIN, t, s) {}
};

//-----------------------------------------------------------------------------
//! This is the base class for surface data. Classes that wish to store data
//! associated with each node or facet of a surface, will use this base class.
class FECORE_API FEPlotSurfaceData : public FEPlotData
{
public:
	FEPlotSurfaceData(FEModel* fem, Var_Type t, Storage_Fmt s) : FEPlotData(fem, FE_REGION_SURFACE, t, s) {}
};

//-----------------------------------------------------------------------------
// helper functions for writing nodal values
void writeNodalValues(FEDomain& dom, FEDataStream& ar, std::function<double(int)> f);
void writeNodalValues(FEDomain& dom, FEDataStream& ar, std::function<mat3ds(int)> f);

void writeNodalValues(FEMesh& mesh, FEDataStream& ar, std::function<double(const FENode& node)> f);
void writeNodalValues(FEMesh& mesh, FEDataStream& ar, std::function<vec3d (const FENode& node)> f);

//-----------------------------------------------------------------------------
// helper functions for writing a single value per element
void writeElementValue(FEDomain& dom, FEDataStream& ar, std::function<double(int nface)> f);
void writeElementValue(FEDomain& dom, FEDataStream& ar, std::function<vec3d (int nface)> f);

//-----------------------------------------------------------------------------
// helper function for summing element values
void writeSummedElementValue(FEDomain& dom, FEDataStream& ar, std::function<double(const FEMaterialPoint& mp)> fnc);
void writeSummedElementValue(FEDomain& dom, FEDataStream& ar, std::function<vec3d (const FEMaterialPoint& mp)> fnc);

//-----------------------------------------------------------------------------
// helper functions for writing averaged element values
void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<double  (const FEMaterialPoint& mp)> fnc);
void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<vec3d   (const FEMaterialPoint& mp)> fnc);
void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<mat3ds  (const FEMaterialPoint& mp)> fnc);
void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<tens4ds (const FEMaterialPoint& mp)> fnc);

//-----------------------------------------------------------------------------
// helper functions for writing averaged element values
void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<mat3ds(FEElement& el, int ip)> fnc);

//-----------------------------------------------------------------------------
// helper functions for writing averaged element values using a filter.
void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<vec3d   (const FEMaterialPoint& mp)> fnc, std::function<double (const vec3d&    m)> flt);
void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<vec3d   (const FEMaterialPoint& mp)> fnc, std::function<vec3d  (const vec3d&    m)> flt);
void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<mat3d   (const FEMaterialPoint& mp)> fnc, std::function<double (const mat3d&    m)> flt);
void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<mat3ds  (const FEMaterialPoint& mp)> fnc, std::function<double (const mat3ds&   m)> flt);
void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<tens3drs(const FEMaterialPoint& mp)> fnc, std::function<double (const tens3drs& m)> flt);

//-----------------------------------------------------------------------------
// helper functions for writing averaged element values using a filter
void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<mat3ds(FEElement& el, int ip)> fnc, std::function<double (const mat3ds& m)> flt);

//-----------------------------------------------------------------------------
// helper function for writing integrals over elements
void writeIntegratedElementValue(FESolidDomain& dom, FEDataStream& ar, std::function<double (const FEMaterialPoint& mp)> fnc);
void writeIntegratedElementValue(FESolidDomain& dom, FEDataStream& ar, std::function<vec3d  (const FEMaterialPoint& mp)> fnc);

//-----------------------------------------------------------------------------
// helper functions for writing SPR projected element values
// TODO: I needed to give these functions a different name because of the implicit conversion between mat3ds and mat3dd
void writeSPRElementValueMat3dd(FESolidDomain& dom, FEDataStream& ar, std::function<mat3dd (const FEMaterialPoint&)> fnc, int interpolOrder = -1);
void writeSPRElementValueMat3ds(FESolidDomain& dom, FEDataStream& ar, std::function<mat3ds (const FEMaterialPoint&)> fnc, int interpolOrder = -1);

//-----------------------------------------------------------------------------
// helper functions for writing nodal projected element values
void writeNodalProjectedElementValues(FEDomain& dom, FEDataStream& ar, std::function<double (const FEMaterialPoint&)> var);
void writeNodalProjectedElementValues(FEDomain& dom, FEDataStream& ar, std::function<vec3d  (const FEMaterialPoint&)> var);
void writeNodalProjectedElementValues(FEDomain& dom, FEDataStream& ar, std::function<mat3ds (const FEMaterialPoint&)> var);

void writeNodalProjectedElementValues(FESurface& dom, FEDataStream& ar, std::function<double(const FEMaterialPoint&)> var);
void writeNodalProjectedElementValues(FESurface& dom, FEDataStream& ar, std::function<vec3d (const FEMaterialPoint&)> var);

void writeNodalValues(FENodeSet& set, FEDataStream& ar, std::function<double(const FEMaterialPoint&)> var);
