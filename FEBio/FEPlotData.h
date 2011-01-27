#pragma once

#include "DumpFile.h"
#include "FEMesh.h"
#include "Archive.h"
#include "FESlidingInterface.h"
#include "FESlidingInterface2.h"
#include "FEFacet2FacetSliding.h"
#include "FETiedInterface.h"
class FEM;

// --- data types ---
enum Var_Type { FLOAT, VEC3F, MAT3FS };

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
	virtual bool Save(FEMesh& m, vector<float>& a) = 0;
};

//-----------------------------------------------------------------------------
//! This is the base class for domain data. Classes that wish to store data
//! associated with each element or node of a domain, will use this base class.
class FEDomainData : public FEPlotData
{
public:
	FEDomainData(Var_Type t, Storage_Fmt s) : FEPlotData(t, s) {}
	void Save(FEM& fem, Archive& ar);
	virtual bool Save(FEDomain& D, vector<float>& a) = 0;
};

//-----------------------------------------------------------------------------
//! This is the base class for surface data. Classes that wish to store data
//! associated with each node or facet of a surface, will use this base class.
class FESurfaceData : public FEPlotData
{
public:
	FESurfaceData(Var_Type t, Storage_Fmt s) : FEPlotData(t, s) {}
	void Save(FEM& fem, Archive& ar);
	virtual bool Save(FESurface& S, vector<float>& a) = 0;
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
//! Nodal effective fluid pressures
class FEPlotEffectiveFluidPressure : public FEDomainData
	{
	public:
		FEPlotEffectiveFluidPressure() : FEDomainData(FLOAT, FMT_NODE){}
		bool Save(FEDomain& m, vector<float>& a);
	};

//-----------------------------------------------------------------------------
//! Nodal effective solute concentrations
class FEPlotEffectiveSoluteConcentration : public FEDomainData
	{
	public:
		FEPlotEffectiveSoluteConcentration() : FEDomainData(FLOAT, FMT_NODE){}
		bool Save(FEDomain& m, vector<float>& a);
	};

//=============================================================================
//                        D O M A I N   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Element stresses
class FEPlotElementStress : public FEDomainData
{
public:
	FEPlotElementStress() : FEDomainData(MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& dom, vector<float>& a);

protected:
	bool WriteSolidStress(FEElasticSolidDomain& d, vector<float>& a);
	bool WriteShellStress(FEElasticShellDomain& d, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Relative volume
class FEPlotRelativeVolume : public FEDomainData
	{
	public:
		FEPlotRelativeVolume() : FEDomainData(FLOAT, FMT_ITEM){}
		bool Save(FEDomain& dom, vector<float>& a);
	};

//-----------------------------------------------------------------------------
//! Actual fluid pressure
class FEPlotActualFluidPressure : public FEDomainData
	{
	public:
		FEPlotActualFluidPressure() : FEDomainData(FLOAT, FMT_ITEM){}
		bool Save(FEDomain& dom, vector<float>& a);
	};

//-----------------------------------------------------------------------------
//! Fluid flux
class FEPlotFluidFlux : public FEDomainData
{
public:
	FEPlotFluidFlux() : FEDomainData(VEC3F, FMT_ITEM){}
	bool Save(FEDomain& dom, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Actual solute concentration
class FEPlotActualSoluteConcentration : public FEDomainData
	{
	public:
		FEPlotActualSoluteConcentration() : FEDomainData(FLOAT, FMT_ITEM){}
		bool Save(FEDomain& dom, vector<float>& a);
	};

//-----------------------------------------------------------------------------
//! Solute flux
class FEPlotSoluteFlux : public FEDomainData
	{
	public:
		FEPlotSoluteFlux() : FEDomainData(VEC3F, FMT_ITEM){}
		bool Save(FEDomain& dom, vector<float>& a);
	};

//-----------------------------------------------------------------------------
//! Material fibers
class FEPlotFiberVector : public FEDomainData
{
public:
	FEPlotFiberVector() : FEDomainData(VEC3F, FMT_ITEM){}
	bool Save(FEDomain& dom, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Shell thicknesses
class FEPlotShellThickness : public FEDomainData
{
public:
	FEPlotShellThickness() : FEDomainData(FLOAT, FMT_MULT){}
	bool Save(FEDomain& dom, vector<float>& a);
};
