#pragma once
#include "FECore/FEPlotData.h"
#include "FEElasticShellDomain.h"
#include "FEElasticSolidDomain.h"
#include "FELinearSolidDomain.h"

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
	bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal velocities
//!
class FEPlotNodeVelocity : public FENodeData
{
public:
	FEPlotNodeVelocity(FEModel* pfem) : FENodeData(PLT_VEC3F, FMT_NODE){}
	bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal accelerations
//!
class FEPlotNodeAcceleration : public FENodeData
{
public:
	FEPlotNodeAcceleration(FEModel* pfem) : FENodeData(PLT_VEC3F, FMT_NODE){}
	bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal reaction forces
class FEPlotNodeReactionForces : public FENodeData
{
public:
	FEPlotNodeReactionForces(FEModel* pfem) : FENodeData(PLT_VEC3F, FMT_NODE){}
	bool Save(FEMesh& m, FEDataStream& a);
};

//=============================================================================
//                         S U R F A C E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Contact gap
//!
class FEPlotContactGap : public FESurfaceData
{
public:
	FEPlotContactGap(FEModel* pfem) : FESurfaceData(PLT_FLOAT, FMT_MULT){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Contact pressure
//!
class FEPlotContactPressure : public FESurfaceData
{
public:
	FEPlotContactPressure(FEModel* pfem) : FESurfaceData(PLT_FLOAT, FMT_MULT){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Contact traction
//!
class FEPlotContactTraction : public FESurfaceData
{
public:
	FEPlotContactTraction(FEModel* pfem) : FESurfaceData(PLT_VEC3F, FMT_MULT){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Contact force
//!
class FEPlotContactForce : public FESurfaceData
{
public:
	FEPlotContactForce(FEModel* pfem) : FESurfaceData(PLT_VEC3F, FMT_REGION){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Contact area
//!
class FEPlotContactArea : public FESurfaceData
{
public:
	FEPlotContactArea(FEModel* pfem) : FESurfaceData(PLT_FLOAT, FMT_MULT){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Contact penalty parameter
class FEPlotContactPenalty : public FESurfaceData
{
public:
	FEPlotContactPenalty(FEModel* pfem) : FESurfaceData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Mortar gap
class FEPlotMortarContactGap : public FESurfaceData
{
public:
	FEPlotMortarContactGap(FEModel* pfem) : FESurfaceData(PLT_FLOAT, FMT_NODE){}
	bool Save(FESurface& S, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Mortar gap
class FEPlotMortarContactGapVector : public FESurfaceData
{
public:
	FEPlotMortarContactGapVector(FEModel* pfem) : FESurfaceData(PLT_VEC3F, FMT_NODE){}
	bool Save(FESurface& S, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Mortar gap
class FEPlotMortarContactNormal : public FESurfaceData
{
public:
	FEPlotMortarContactNormal(FEModel* pfem) : FESurfaceData(PLT_VEC3F, FMT_NODE){}
	bool Save(FESurface& S, FEDataStream& a);
};

//=============================================================================
//							D O M A I N   D A T A
//=============================================================================
//-----------------------------------------------------------------------------
//! Element norm for G
class FEPlotElementGnorm : public FEDomainData
{
public:
	FEPlotElementGnorm(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element stresses
class FEPlotElementStress : public FEDomainData
{
public:
	FEPlotElementStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element norm for Cauchy stress
class FEPlotElementsnorm : public FEDomainData
{
public:
	FEPlotElementsnorm(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element norm for Cauchy stress moment
class FEPlotElementtaunorm : public FEDomainData
{
public:
	FEPlotElementtaunorm(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element norm for PK1 stress
class FEPlotElementPK1norm : public FEDomainData
{
public:
	FEPlotElementPK1norm(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element norm for PK1 stress moment
class FEPlotElementQK1norm : public FEDomainData
{
public:
	FEPlotElementQK1norm(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element norm for PK2 stress
class FEPlotElementSnorm : public FEDomainData
{
public:
	FEPlotElementSnorm(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element norm for PK2 stress moment
class FEPlotElementTnorm : public FEDomainData
{
public:
	FEPlotElementTnorm(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element infinitesimal strain gradiet norm
class FEPlotElementinfstrnorm : public FEDomainData
{
public:
	FEPlotElementinfstrnorm(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element Green-Lagrange strain gradient norm
class FEPlotElementGLstrnorm : public FEDomainData
{
public:
	FEPlotElementGLstrnorm(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element Euler-Almansi strain gradient norm
class FEPlotElementEAstrnorm : public FEDomainData
{
public:
	FEPlotElementEAstrnorm(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element macro energy
class FEPlotElementMacroEnergy : public FEDomainData
{
public:
	FEPlotElementMacroEnergy(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element micro energy
class FEPlotElementMicroEnergy : public FEDomainData
{
public:
	FEPlotElementMicroEnergy(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};


//-----------------------------------------------------------------------------
//! Element difference between macro and micro energy
class FEPlotElementenergydiff : public FEDomainData
{
public:
	FEPlotElementenergydiff(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Strain energy density
class FEPlotStrainEnergyDensity : public FEDomainData
{
public:
	FEPlotStrainEnergyDensity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Deviatoric strain energy density
class FEPlotDevStrainEnergyDensity : public FEDomainData
{
public:
	FEPlotDevStrainEnergyDensity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Specific strain energy
class FEPlotSpecificStrainEnergy : public FEDomainData
{
public:
	FEPlotSpecificStrainEnergy(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Mass density
class FEPlotDensity : public FEDomainData
{
public:
	FEPlotDensity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Relative volume
class FEPlotRelativeVolume : public FEDomainData
{
public:
	FEPlotRelativeVolume(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Material fibers
class FEPlotFiberVector : public FEDomainData
{
public:
	FEPlotFiberVector(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! fiber stretch
class FEPlotFiberStretch : public FEDomainData
{
public:
	FEPlotFiberStretch(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! deviatoric fiber stretch
class FEPlotDevFiberStretch : public FEDomainData
{
public:
	FEPlotDevFiberStretch(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Shell thicknesses
class FEPlotShellThickness : public FEDomainData
{
public:
	FEPlotShellThickness(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_MULT){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element elasticity tensor
class FEPlotElementElasticity : public FEDomainData
{
public:
	FEPlotElementElasticity(FEModel* pfem) : FEDomainData(PLT_TENS4FS, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Damage reduction factor
class FEPlotDamage : public FEDomainData
{
public:
	FEPlotDamage(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Mixture volume fraction
class FEPlotMixtureVolumeFraction : public FEDomainData
{
public:
	FEPlotMixtureVolumeFraction(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Class that outputs the element nodal stresses for UT4 domains
class FEPlotUT4NodalStresses : public FEDomainData
{
public:
	FEPlotUT4NodalStresses(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_NODE) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Store shell strains
class FEPlotShellStrain : public FEDomainData
{
public:
	FEPlotShellStrain(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! class the projects stresses from integration points to nodes using
//! SPR (superconvergergent patch recovery)
class FEPlotSPRStresses : public FEDomainData
{
public:
	FEPlotSPRStresses(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_NODE){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! class the projects stresses from integration points to nodes using
//! SPR (superconvergergent patch recovery)
class FEPlotSPRLinearStresses : public FEDomainData
{
public:
	FEPlotSPRLinearStresses(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_NODE){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! class that projects the principal stresses from integration points to nodes
//! using the SPR projection method
class FEPlotSPRPrincStresses : public FEDomainData
{
public:
	FEPlotSPRPrincStresses(FEModel* pfem) : FEDomainData(PLT_MAT3FD, FMT_NODE){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! class the projects stresses from integration points to nodes using
//! SPR (superconvergergent patch recovery)
class FEPlotSPRTestLinear: public FEDomainData
{
public:
	FEPlotSPRTestLinear(FEModel* pfem) : FEDomainData(PLT_MAT3FD, FMT_NODE){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! class the projects stresses from integration points to nodes using
//! SPR (superconvergergent patch recovery)
class FEPlotSPRTestQuadratic: public FEDomainData
{
public:
	FEPlotSPRTestQuadratic(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_NODE){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Rigid body displacement
class FEPlotRigidDisplacement : public FEDomainData
{
public:
	FEPlotRigidDisplacement(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_REGION),  m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rigid body velocity
class FEPlotRigidVelocity : public FEDomainData
{
public:
	FEPlotRigidVelocity(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_REGION),  m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rigid body acceleration
class FEPlotRigidAcceleration : public FEDomainData
{
public:
	FEPlotRigidAcceleration(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_REGION),  m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rigid body rotation
class FEPlotRigidRotation : public FEDomainData
{
public:
	FEPlotRigidRotation(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_REGION),  m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rigid body angular velocity
class FEPlotRigidAngularVelocity : public FEDomainData
{
public:
	FEPlotRigidAngularVelocity(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_REGION),  m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rigid body angular acceleration
class FEPlotRigidAngularAcceleration : public FEDomainData
{
public:
	FEPlotRigidAngularAcceleration(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_REGION),  m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rigid body kinetic energy
class FEPlotRigidKineticEnergy : public FEDomainData
{
public:
	FEPlotRigidKineticEnergy(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_REGION),  m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rigid Euler angles
class FEPlotRigidEuler : public FEDomainData
{
public:
	FEPlotRigidEuler(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_REGION), m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rotation vector
class FEPlotRigidRotationVector : public FEDomainData
{
public:
	FEPlotRigidRotationVector(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_REGION), m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Class that projects stresses from integration points to the nodes
//! TODO: This only works with tet10 and hex8 -domains
class FEPlotNodalStresses : public FEDomainData
{
public:
	FEPlotNodalStresses(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_MULT){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Lagrange strains
class FEPlotLagrangeStrain : public FEDomainData
{
public:
	FEPlotLagrangeStrain(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Rigid body reaction force
class FEPlotRigidReactionForce : public FEDomainData
{
public:
	FEPlotRigidReactionForce(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_REGION),  m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rigid body reaction torque
class FEPlotRigidReactionTorque : public FEDomainData
{
public:
    FEPlotRigidReactionTorque(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_REGION), m_pfem(pfem) {}
    bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};
