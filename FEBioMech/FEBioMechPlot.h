#pragma once
#include "FECore/FEPlotData.h"
#include "FEBioMech/FEElasticShellDomain.h"
#include "FEBioMech/FEElasticSolidDomain.h"
#include "FEBioMech/FELinearSolidDomain.h"

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
//! Nodal reaction forces
class FEPlotNodeReactionForces : public FENodeData
{
public:
	FEPlotNodeReactionForces(FEModel* pfem) : FENodeData(PLT_VEC3F, FMT_NODE){}
	bool Save(FEMesh& m, vector<float>& a);
};


//=============================================================================
//							D O M A I N   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Element stresses
class FEPlotElementStress : public FEDomainData
{
public:
	FEPlotElementStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& dom, vector<float>& a);

protected:
	bool WriteSolidStress(FEElasticSolidDomain& d, vector<float>& a);
	bool WriteShellStress(FEElasticShellDomain& d, vector<float>& a);

	bool WriteLinearSolidStress(FELinearSolidDomain& d, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Strain energy density
class FEPlotStrainEnergyDensity : public FEDomainData
{
public:
	FEPlotStrainEnergyDensity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Specific strain energy
class FEPlotSpecificStrainEnergy : public FEDomainData
{
public:
	FEPlotSpecificStrainEnergy(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Mass density
class FEPlotDensity : public FEDomainData
{
public:
	FEPlotDensity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Relative volume
class FEPlotRelativeVolume : public FEDomainData
{
public:
	FEPlotRelativeVolume(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Material fibers
class FEPlotFiberVector : public FEDomainData
{
public:
	FEPlotFiberVector(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM){}
	bool Save(FEDomain& dom, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Shell thicknesses
class FEPlotShellThickness : public FEDomainData
{
public:
	FEPlotShellThickness(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_MULT){}
	bool Save(FEDomain& dom, vector<float>& a);
};


//-----------------------------------------------------------------------------
//! Damage reduction factor
class FEPlotDamage : public FEDomainData
{
public:
	FEPlotDamage(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& m, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Mixture volume fraction
class FEPlotMixtureVolumeFraction : public FEDomainData
{
public:
	FEPlotMixtureVolumeFraction(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& m, vector<float>& a);
};

//-----------------------------------------------------------------------------
//! Class that outputs the element nodal stresses for UT4 domains
class FEPlotUT4NodalStresses : public FEDomainData
{
public:
	FEPlotUT4NodalStresses(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_NODE) {}
	bool Save(FEDomain& dom, vector<float>& a);
};
