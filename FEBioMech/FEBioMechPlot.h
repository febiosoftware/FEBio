#pragma once
#include <FECore/FEPlotData.h>

//=============================================================================
//                            N O D E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Nodal velocities
//!
class FEPlotNodeVelocity : public FEPlotNodeData
{
public:
	FEPlotNodeVelocity(FEModel* pfem) : FEPlotNodeData(PLT_VEC3F, FMT_NODE){}
	bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal accelerations
//!
class FEPlotNodeAcceleration : public FEPlotNodeData
{
public:
	FEPlotNodeAcceleration(FEModel* pfem) : FEPlotNodeData(PLT_VEC3F, FMT_NODE){}
	bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal reaction forces
class FEPlotNodeReactionForces : public FEPlotNodeData
{
public:
	FEPlotNodeReactionForces(FEModel* pfem) : FEPlotNodeData(PLT_VEC3F, FMT_NODE){}
	bool Save(FEMesh& m, FEDataStream& a);
};

//=============================================================================
//                         S U R F A C E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Contact gap
//!
class FEPlotContactGap : public FEPlotSurfaceData
{
public:
    FEPlotContactGap(FEModel* pfem) : FEPlotSurfaceData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Vector gap
//!
class FEPlotVectorGap : public FEPlotSurfaceData
{
public:
    FEPlotVectorGap(FEModel* pfem) : FEPlotSurfaceData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Contact pressure
//!
class FEPlotContactPressure : public FEPlotSurfaceData
{
public:
    FEPlotContactPressure(FEModel* pfem) : FEPlotSurfaceData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Contact traction
//!
class FEPlotContactTraction : public FEPlotSurfaceData
{
public:
    FEPlotContactTraction(FEModel* pfem) : FEPlotSurfaceData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal contact gap
//!
class FEPlotNodalContactGap : public FEPlotSurfaceData
{
public:
	FEPlotNodalContactGap(FEModel* pfem) : FEPlotSurfaceData(PLT_FLOAT, FMT_MULT){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal vector gap
//!
class FEPlotNodalVectorGap : public FEPlotSurfaceData
{
public:
    FEPlotNodalVectorGap(FEModel* pfem) : FEPlotSurfaceData(PLT_VEC3F, FMT_MULT){}
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal contact pressure
//!
class FEPlotNodalContactPressure : public FEPlotSurfaceData
{
public:
    FEPlotNodalContactPressure(FEModel* pfem) : FEPlotSurfaceData(PLT_FLOAT, FMT_MULT){}
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal contact traction
//!
class FEPlotNodalContactTraction : public FEPlotSurfaceData
{
public:
	FEPlotNodalContactTraction(FEModel* pfem) : FEPlotSurfaceData(PLT_VEC3F, FMT_MULT){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Surface traction
//!
class FEPlotSurfaceTraction : public FEPlotSurfaceData
{
public:
    FEPlotSurfaceTraction(FEModel* pfem) : FEPlotSurfaceData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal surface traction
//!
class FEPlotNodalSurfaceTraction : public FEPlotSurfaceData
{
public:
    FEPlotNodalSurfaceTraction(FEModel* pfem) : FEPlotSurfaceData(PLT_VEC3F, FMT_MULT){}
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Stick status
//!
class FEPlotStickStatus : public FEPlotSurfaceData
{
public:
    FEPlotStickStatus(FEModel* pfem) : FEPlotSurfaceData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Contact force
//!
class FEPlotContactForce : public FEPlotSurfaceData
{
public:
	FEPlotContactForce(FEModel* pfem) : FEPlotSurfaceData(PLT_VEC3F, FMT_REGION){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Contact area
//!
class FEPlotContactArea : public FEPlotSurfaceData
{
public:
	FEPlotContactArea(FEModel* pfem) : FEPlotSurfaceData(PLT_FLOAT, FMT_MULT){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Contact penalty parameter
class FEPlotContactPenalty : public FEPlotSurfaceData
{
public:
	FEPlotContactPenalty(FEModel* pfem) : FEPlotSurfaceData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Mortar gap
class FEPlotMortarContactGap : public FEPlotSurfaceData
{
public:
	FEPlotMortarContactGap(FEModel* pfem) : FEPlotSurfaceData(PLT_FLOAT, FMT_NODE){}
	bool Save(FESurface& S, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Enclosed volume
//!
class FEPlotEnclosedVolume : public FEPlotSurfaceData
{
private:
    FEModel*            m_pfem;
    bool                m_binit;
    vector<FEElement*>  m_elem;
    vector<vec3d>       m_area;
    
public:
    FEPlotEnclosedVolume(FEModel* pfem) : FEPlotSurfaceData(PLT_FLOAT, FMT_REGION){ m_pfem = pfem; m_binit = true; }
    bool Save(FESurface& surf, FEDataStream& a);
};

//=============================================================================
//							D O M A I N   D A T A
//=============================================================================
//-----------------------------------------------------------------------------
//! Velocity
class FEPlotElementVelocity : public FEPlotDomainData
{
public:
    FEPlotElementVelocity(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Acceleration
class FEPlotElementAcceleration : public FEPlotDomainData
{
public:
    FEPlotElementAcceleration(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element norm for G
class FEPlotElementGnorm : public FEPlotDomainData
{
public:
	FEPlotElementGnorm(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element stresses
class FEPlotElementStress : public FEPlotDomainData
{
public:
	FEPlotElementStress(FEModel* pfem) : FEPlotDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element uncoupled pressure
class FEPlotElementUncoupledPressure : public FEPlotDomainData
{
public:
    FEPlotElementUncoupledPressure(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element norm for Cauchy stress
class FEPlotElementsnorm : public FEPlotDomainData
{
public:
	FEPlotElementsnorm(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element norm for Cauchy stress moment
class FEPlotElementtaunorm : public FEPlotDomainData
{
public:
	FEPlotElementtaunorm(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element norm for PK1 stress
class FEPlotElementPK1norm : public FEPlotDomainData
{
public:
	FEPlotElementPK1norm(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element norm for PK1 stress moment
class FEPlotElementQK1norm : public FEPlotDomainData
{
public:
	FEPlotElementQK1norm(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element norm for PK2 stress
class FEPlotElementSnorm : public FEPlotDomainData
{
public:
	FEPlotElementSnorm(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element norm for PK2 stress moment
class FEPlotElementTnorm : public FEPlotDomainData
{
public:
	FEPlotElementTnorm(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element infinitesimal strain gradiet norm
class FEPlotElementinfstrnorm : public FEPlotDomainData
{
public:
	FEPlotElementinfstrnorm(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element Green-Lagrange strain gradient norm
class FEPlotElementGLstrnorm : public FEPlotDomainData
{
public:
	FEPlotElementGLstrnorm(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element Euler-Almansi strain gradient norm
class FEPlotElementEAstrnorm : public FEPlotDomainData
{
public:
	FEPlotElementEAstrnorm(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element macro energy
class FEPlotElementMacroEnergy : public FEPlotDomainData
{
public:
	FEPlotElementMacroEnergy(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element micro energy
class FEPlotElementMicroEnergy : public FEPlotDomainData
{
public:
	FEPlotElementMicroEnergy(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};


//-----------------------------------------------------------------------------
//! Element difference between macro and micro energy
class FEPlotElementenergydiff : public FEPlotDomainData
{
public:
	FEPlotElementenergydiff(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Strain energy density
class FEPlotStrainEnergyDensity : public FEPlotDomainData
{
public:
	FEPlotStrainEnergyDensity(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Deviatoric strain energy density
class FEPlotDevStrainEnergyDensity : public FEPlotDomainData
{
public:
	FEPlotDevStrainEnergyDensity(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Specific strain energy
class FEPlotSpecificStrainEnergy : public FEPlotDomainData
{
public:
	FEPlotSpecificStrainEnergy(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Kinetic energy density
class FEPlotKineticEnergyDensity : public FEPlotDomainData
{
public:
    FEPlotKineticEnergyDensity(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Mass density
class FEPlotDensity : public FEPlotDomainData
{
public:
	FEPlotDensity(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Strain energy
class FEPlotElementStrainEnergy : public FEPlotDomainData
{
public:
    FEPlotElementStrainEnergy(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Kinetic energy
class FEPlotElementKineticEnergy : public FEPlotDomainData
{
public:
    FEPlotElementKineticEnergy(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Center of mass
class FEPlotElementCenterOfMass : public FEPlotDomainData
{
public:
    FEPlotElementCenterOfMass(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Linear momentum
class FEPlotElementLinearMomentum : public FEPlotDomainData
{
public:
    FEPlotElementLinearMomentum(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Angular momentum
class FEPlotElementAngularMomentum : public FEPlotDomainData
{
public:
    FEPlotElementAngularMomentum(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Stress power
class FEPlotElementStressPower : public FEPlotDomainData
{
public:
    FEPlotElementStressPower(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Strain energy at current time
class FEPlotCurrentElementStrainEnergy : public FEPlotDomainData
{
public:
    FEPlotCurrentElementStrainEnergy(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Kinetic energy at current time
class FEPlotCurrentElementKineticEnergy : public FEPlotDomainData
{
public:
    FEPlotCurrentElementKineticEnergy(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Center of mass at current time
class FEPlotCurrentElementCenterOfMass : public FEPlotDomainData
{
public:
    FEPlotCurrentElementCenterOfMass(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Linear momentum at current time
class FEPlotCurrentElementLinearMomentum : public FEPlotDomainData
{
public:
    FEPlotCurrentElementLinearMomentum(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Angular momentum at current time
class FEPlotCurrentElementAngularMomentum : public FEPlotDomainData
{
public:
    FEPlotCurrentElementAngularMomentum(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Relative volume
class FEPlotRelativeVolume : public FEPlotDomainData
{
public:
	FEPlotRelativeVolume(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Material fibers
class FEPlotFiberVector : public FEPlotDomainData
{
public:
	FEPlotFiberVector(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Material axes
class FEPlotMaterialAxes : public FEPlotDomainData
{
public:
	FEPlotMaterialAxes(FEModel* pfem) : FEPlotDomainData(PLT_MAT3F, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! fiber stretch
class FEPlotFiberStretch : public FEPlotDomainData
{
public:
	FEPlotFiberStretch(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! deviatoric fiber stretch
class FEPlotDevFiberStretch : public FEPlotDomainData
{
public:
	FEPlotDevFiberStretch(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Shell thicknesses
class FEPlotShellThickness : public FEPlotDomainData
{
public:
	FEPlotShellThickness(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_MULT){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Shell directors
class FEPlotShellDirector : public FEPlotDomainData
{
public:
	FEPlotShellDirector(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_MULT){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element elasticity tensor
class FEPlotElementElasticity : public FEPlotDomainData
{
public:
	FEPlotElementElasticity(FEModel* pfem) : FEPlotDomainData(PLT_TENS4FS, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Damage reduction factor
class FEPlotDamage : public FEPlotDomainData
{
public:
	FEPlotDamage(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Damage reduction factor
class FEPlotNestedDamage : public FEPlotDomainData
{
public:
    FEPlotNestedDamage(FEModel* pfem);
    bool Save(FEDomain& m, FEDataStream& a);
    bool SetFilter(int nsol);
protected:
    int			m_nmat;
    FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
//! Mixture volume fraction
class FEPlotMixtureVolumeFraction : public FEPlotDomainData
{
public:
	FEPlotMixtureVolumeFraction(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Class that outputs the element nodal stresses for UT4 domains
class FEPlotUT4NodalStresses : public FEPlotDomainData
{
public:
	FEPlotUT4NodalStresses(FEModel* pfem) : FEPlotDomainData(PLT_MAT3FS, FMT_NODE) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Store shell strains
class FEPlotShellStrain : public FEPlotDomainData
{
public:
	FEPlotShellStrain(FEModel* pfem) : FEPlotDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Shell relative volume
class FEPlotShellRelativeVolume : public FEPlotDomainData
{
public:
    FEPlotShellRelativeVolume(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! class the projects stresses from integration points to nodes using
//! SPR (superconvergergent patch recovery)
class FEPlotSPRStresses : public FEPlotDomainData
{
public:
	FEPlotSPRStresses(FEModel* pfem) : FEPlotDomainData(PLT_MAT3FS, FMT_NODE){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! class the projects stresses from integration points to nodes using
//! SPR (superconvergergent patch recovery)
class FEPlotSPRLinearStresses : public FEPlotDomainData
{
public:
	FEPlotSPRLinearStresses(FEModel* pfem) : FEPlotDomainData(PLT_MAT3FS, FMT_NODE){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! class that projects the principal stresses from integration points to nodes
//! using the SPR projection method
class FEPlotSPRPrincStresses : public FEPlotDomainData
{
public:
	FEPlotSPRPrincStresses(FEModel* pfem) : FEPlotDomainData(PLT_MAT3FD, FMT_NODE){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! class the projects stresses from integration points to nodes using
//! SPR (superconvergergent patch recovery)
class FEPlotSPRTestLinear: public FEPlotDomainData
{
public:
	FEPlotSPRTestLinear(FEModel* pfem) : FEPlotDomainData(PLT_MAT3FD, FMT_NODE){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! class the projects stresses from integration points to nodes using
//! SPR (superconvergergent patch recovery)
class FEPlotSPRTestQuadratic: public FEPlotDomainData
{
public:
	FEPlotSPRTestQuadratic(FEModel* pfem) : FEPlotDomainData(PLT_MAT3FS, FMT_NODE){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Rigid body displacement
class FEPlotRigidDisplacement : public FEPlotDomainData
{
public:
	FEPlotRigidDisplacement(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_REGION),  m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rigid body velocity
class FEPlotRigidVelocity : public FEPlotDomainData
{
public:
	FEPlotRigidVelocity(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_REGION),  m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rigid body acceleration
class FEPlotRigidAcceleration : public FEPlotDomainData
{
public:
	FEPlotRigidAcceleration(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_REGION),  m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rigid body rotation
class FEPlotRigidRotation : public FEPlotDomainData
{
public:
	FEPlotRigidRotation(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_REGION),  m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rigid body angular velocity
class FEPlotRigidAngularVelocity : public FEPlotDomainData
{
public:
	FEPlotRigidAngularVelocity(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_REGION),  m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rigid body angular acceleration
class FEPlotRigidAngularAcceleration : public FEPlotDomainData
{
public:
	FEPlotRigidAngularAcceleration(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_REGION),  m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rigid body kinetic energy
class FEPlotRigidKineticEnergy : public FEPlotDomainData
{
public:
	FEPlotRigidKineticEnergy(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_REGION),  m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rigid body linear momentum
class FEPlotRigidLinearMomentum : public FEPlotDomainData
{
public:
    FEPlotRigidLinearMomentum(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_REGION),  m_pfem(pfem) {}
    bool Save(FEDomain& dom, FEDataStream& a);
private:
    FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rigid body angular momentum
class FEPlotRigidAngularMomentum : public FEPlotDomainData
{
public:
    FEPlotRigidAngularMomentum(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_REGION),  m_pfem(pfem) {}
    bool Save(FEDomain& dom, FEDataStream& a);
private:
    FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rigid Euler angles
class FEPlotRigidEuler : public FEPlotDomainData
{
public:
	FEPlotRigidEuler(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_REGION), m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rotation vector
class FEPlotRigidRotationVector : public FEPlotDomainData
{
public:
	FEPlotRigidRotationVector(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_REGION), m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Class that projects stresses from integration points to the nodes
//! TODO: This only works with tet10 and hex8 -domains
class FEPlotNodalStresses : public FEPlotDomainData
{
public:
	FEPlotNodalStresses(FEModel* pfem) : FEPlotDomainData(PLT_MAT3FS, FMT_MULT){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Lagrange strains
class FEPlotLagrangeStrain : public FEPlotDomainData
{
public:
	FEPlotLagrangeStrain(FEModel* pfem) : FEPlotDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Lagrange strains
class FEPlotSPRLagrangeStrain : public FEPlotDomainData
{
public:
	FEPlotSPRLagrangeStrain(FEModel* pfem) : FEPlotDomainData(PLT_MAT3FS, FMT_NODE){}
	bool Save(FEDomain& dom, FEDataStream& a);
};


//-----------------------------------------------------------------------------
//! Rigid body reaction force
class FEPlotRigidReactionForce : public FEPlotDomainData
{
public:
	FEPlotRigidReactionForce(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_REGION),  m_pfem(pfem) {}
	bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};

//-----------------------------------------------------------------------------
//! Rigid body reaction torque
class FEPlotRigidReactionTorque : public FEPlotDomainData
{
public:
    FEPlotRigidReactionTorque(FEModel* pfem) : FEPlotDomainData(PLT_VEC3F, FMT_REGION), m_pfem(pfem) {}
    bool Save(FEDomain& dom, FEDataStream& a);
private:
	FEModel* m_pfem;
};
