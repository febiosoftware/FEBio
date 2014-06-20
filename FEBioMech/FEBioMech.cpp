#include "stdafx.h"

#include "FEBioMech.h"
#include "FE2DFiberNeoHookean.h"
#include "FE2DTransIsoMooneyRivlin.h"
#include "FE2DTransIsoVerondaWestmann.h"
#include "FEArrudaBoyce.h"
#include "FECarterHayesOld.h"
#include "FECellGrowth.h"
#include "FEDamageMooneyRivlin.h"
#include "FEDamageNeoHookean.h"
#include "FEDamageTransIsoMooneyRivlin.h"
#include "FESpringMaterial.h"
#include "FEDonnanEquilibrium.h"
#include "FEEFD.h"
#include "FEEFDDonnanEquilibrium.h"
#include "FEEFDMooneyRivlin.h"
#include "FEEFDNeoHookean.h"
#include "FEEFDUncoupled.h"
#include "FEEFDVerondaWestmann.h"
#include "FEElasticMixture.h"
#include "FEElasticMultigeneration.h"
#include "FEEllipsoidalFiberDistribution.h"
#include "FEFiberExpPow.h"
#include "FEFiberExpPowUncoupled.h"
#include "FEFiberNeoHookean.h"
#include "FEFungOrthoCompressible.h"
#include "FEFungOrthotropic.h"
#include "FEGasserOgdenHolzapfel.h"
#include "FEHolmesMow.h"
#include "FEHuiskesSupply.h"
#include "FEIncompNeoHookean.h"
#include "FEIsotropicElastic.h"
#include "FELinearElastic.h"
#include "FELinearOrthotropic.h"
#include "FELinearTransIso.h"
#include "FEMooneyRivlin.h"
#include "FEMRVonMisesFibers.h"
#include "FEMuscleMaterial.h"
#include "FENeoHookean.h"
#include "FENeoHookeanTransIso.h"
#include "FEOgdenMaterial.h"
#include "FEOgdenUnconstrained.h"
#include "FEOrthoElastic.h"
#include "FEPerfectOsmometer.h"
#include "FERemodelingElasticMaterial.h"
#include "FERigidMaterial.h"
#include "FESphericalFiberDistribution.h"
#include "FEStVenantKirchhoff.h"
#include "FETCNonlinearOrthotropic.h"
#include "FETendonMaterial.h"
#include "FETransIsoMooneyRivlin.h"
#include "FETransIsoVerondaWestmann.h"
#include "FETrussMaterial.h"
#include "FEUncoupledActiveContraction.h"
#include "FEUncoupledElasticMixture.h"
#include "FEUncoupledViscoElasticMaterial.h"
#include "FEVerondaWestmann.h"
#include "FEViscoElasticMaterial.h"
#include "FEVonMisesPlasticity.h"
#include "FEPreStrainTransIsoMR.h"
#include "FEPreStrainCoupledTransIsoMR.h"
#include "FEPreStrainElastic.h"
#include "FEElasticFiberMaterial.h"
#include "FEFiberDensityDistribution.h"
#include "FEContinuousFiberDistribution.h"
#include "FEFiberIntegrationGauss.h"
#include "FEFiberIntegrationTrapezoidal.h"
#include "FEFiberIntegrationGeodesic.h"
#include "FECoupledTransIsoMooneyRivlin.h"
#include "FECoupledTransIsoVerondaWestmann.h"
#include "FEMicroMaterial.h"
#include "FESpringMaterial.h"
#include "FEElasticMultigeneration.h"
#include "FEPRLig.h"
#include "FECoupledMooneyRivlin.h"
#include "FECoupledVerondaWestmann.h"

#include "FEPressureLoad.h"
#include "FETractionLoad.h"
#include "FEConstBodyForce.h"
#include "FEPointBodyForce.h"

#include "FEFacet2FacetSliding.h"
#include "FEPeriodicBoundary.h"
#include "FERigidWallInterface.h"
#include "FESlidingInterface.h"
#include "FESlidingInterfaceBW.h"
#include "FESurfaceConstraint.h"
#include "FETiedInterface.h"
#include "FEStickyInterface.h"
#include "FEInSituStretch.h"
#include "FEPointConstraint.h"
#include "FEFacet2FacetTied.h"

#include "FEAugLagLinearConstraint.h"
#include "FERigidJoint.h"
#include "FERigidSphericalJoint.h"
#include "FERigidPinJoint.h"

#include "FESolidAnalysis.h"
#include "FESolidSolver.h"
#include "FELinearSolidSolver.h"
#include "FEExplicitSolidSolver.h"
#include "FECGSolidSolver.h"

#include "FEBioMechPlot.h"
#include "FEBioMechData.h"

#include "FESolidDomainFactory.h"

//-----------------------------------------------------------------------------
//! Register all the classes of the FEBioMech module with the FEBio framework.
void FEBioMech::InitModule()
{
//-----------------------------------------------------------------------------
// Domain factory
	FECoreKernel& febio = FECoreKernel::GetInstance();
	febio.RegisterDomain(new FESolidDomainFactory);

//-----------------------------------------------------------------------------
// Analysis classes
REGISTER_FECORE_CLASS(FESolidAnalysis        , FEANALYSIS_ID, "solid"         );
REGISTER_FECORE_CLASS(FEExplicitSolidAnalysis, FEANALYSIS_ID, "explicit-solid");
REGISTER_FECORE_CLASS(FELinearSolidAnalysis  , FEANALYSIS_ID, "linear-solid"  );

//-----------------------------------------------------------------------------
// Solver classes
REGISTER_FECORE_CLASS(FESolidSolver        , FESOLVER_ID, "solid"         );
REGISTER_FECORE_CLASS(FEExplicitSolidSolver, FESOLVER_ID, "explicit-solid");
REGISTER_FECORE_CLASS(FELinearSolidSolver  , FESOLVER_ID, "linear-solid"  );
REGISTER_FECORE_CLASS(FECGSolidSolver      , FESOLVER_ID, "CG-solid"      );

//-----------------------------------------------------------------------------
// material classes
REGISTER_FECORE_CLASS(FE2DFiberNeoHookean            ,FEMATERIAL_ID, "2D fiber neo-Hookean"          );
REGISTER_FECORE_CLASS(FE2DTransIsoMooneyRivlin       ,FEMATERIAL_ID, "2D trans iso Mooney-Rivlin"    );
REGISTER_FECORE_CLASS(FE2DTransIsoVerondaWestmann    ,FEMATERIAL_ID, "2D trans iso Veronda-Westmann" );
REGISTER_FECORE_CLASS(FEArrudaBoyce                  ,FEMATERIAL_ID, "Arruda-Boyce"                  );
REGISTER_FECORE_CLASS(FECellGrowth                   ,FEMATERIAL_ID, "cell growth"                   );
REGISTER_FECORE_CLASS(FEDamageMooneyRivlin           ,FEMATERIAL_ID, "damage Mooney-Rivlin"          );
REGISTER_FECORE_CLASS(FEDamageNeoHookean             ,FEMATERIAL_ID, "damage neo-Hookean"            );
REGISTER_FECORE_CLASS(FEDamageTransIsoMooneyRivlin   ,FEMATERIAL_ID, "damage trans iso Mooney-Rivlin");
REGISTER_FECORE_CLASS(FEDonnanEquilibrium            ,FEMATERIAL_ID, "Donnan equilibrium"            );
REGISTER_FECORE_CLASS(FEEFD                          ,FEMATERIAL_ID, "EFD"                           );
REGISTER_FECORE_CLASS(FEEFDDonnanEquilibrium         ,FEMATERIAL_ID, "EFD Donnan equilibrium"        );
REGISTER_FECORE_CLASS(FEEFDMooneyRivlin              ,FEMATERIAL_ID, "EFD Mooney-Rivlin"             );
//REGISTER_FECORE_CLASS(FEEFDNeoHookean                ,FEMATERIAL_ID, "EFD neo-Hookean"               );
REGISTER_FECORE_CLASS(FEEFDNeoHookeanOld             ,FEMATERIAL_ID, "EFD neo-Hookean"               );
REGISTER_FECORE_CLASS(FEEFDUncoupled                 ,FEMATERIAL_ID, "EFD uncoupled"                 );
REGISTER_FECORE_CLASS(FEEFDVerondaWestmann           ,FEMATERIAL_ID, "EFD Veronda-Westmann"          );
REGISTER_FECORE_CLASS(FEElasticMixture               ,FEMATERIAL_ID, "solid mixture"                 );
REGISTER_FECORE_CLASS(FEEllipsoidalFiberDistribution ,FEMATERIAL_ID, "ellipsoidal fiber distribution");
REGISTER_FECORE_CLASS(FEEllipsoidalFiberDistributionOld,FEMATERIAL_ID, "ellipsoidal fiber distribution (old)");
REGISTER_FECORE_CLASS(FEFiberExpPow                  ,FEMATERIAL_ID, "fiber-exp-pow"                 );
REGISTER_FECORE_CLASS(FEFiberExpPowUncoupled         ,FEMATERIAL_ID, "fiber-exp-pow-uncoupled"       );
REGISTER_FECORE_CLASS(FEFiberNeoHookean              ,FEMATERIAL_ID, "fiber neo-Hookean"             );
REGISTER_FECORE_CLASS(FEFungOrthoCompressible        ,FEMATERIAL_ID, "Fung-ortho-compressible"       );
REGISTER_FECORE_CLASS(FEFungOrthotropic              ,FEMATERIAL_ID, "Fung orthotropic"              );
REGISTER_FECORE_CLASS(FEGasserOgdenHolzapfel         ,FEMATERIAL_ID, "Gasser-Ogden-Holzapfel"        );
REGISTER_FECORE_CLASS(FEHolmesMow                    ,FEMATERIAL_ID, "Holmes-Mow"                    );
REGISTER_FECORE_CLASS(FEIncompNeoHookean             ,FEMATERIAL_ID, "incomp neo-Hookean"            );
REGISTER_FECORE_CLASS(FEIsotropicElastic             ,FEMATERIAL_ID, "isotropic elastic"             );
REGISTER_FECORE_CLASS(FELinearElastic                ,FEMATERIAL_ID, "linear elastic"                );
REGISTER_FECORE_CLASS(FELinearOrthotropic            ,FEMATERIAL_ID, "linear orthotropic"            );
REGISTER_FECORE_CLASS(FELinearTransIso               ,FEMATERIAL_ID, "linear trans iso"              );
REGISTER_FECORE_CLASS(FEMooneyRivlin                 ,FEMATERIAL_ID, "Mooney-Rivlin"                 );
REGISTER_FECORE_CLASS(FECoupledMooneyRivlin          ,FEMATERIAL_ID, "coupled Mooney-Rivlin"         );
REGISTER_FECORE_CLASS(FECoupledVerondaWestmann       ,FEMATERIAL_ID, "coupled Veronda-Westmann"      );
REGISTER_FECORE_CLASS(FEMuscleMaterial               ,FEMATERIAL_ID, "muscle material"               );
REGISTER_FECORE_CLASS(FENeoHookean                   ,FEMATERIAL_ID, "neo-Hookean"                   );
REGISTER_FECORE_CLASS(FENeoHookeanTransIso           ,FEMATERIAL_ID, "neo-Hookean transiso"          );
REGISTER_FECORE_CLASS(FEOgdenMaterial                ,FEMATERIAL_ID, "Ogden"                         );
REGISTER_FECORE_CLASS(FEOgdenUnconstrained           ,FEMATERIAL_ID, "Ogden unconstrained"           );
REGISTER_FECORE_CLASS(FEOrthoElastic                 ,FEMATERIAL_ID, "orthotropic elastic"           );
REGISTER_FECORE_CLASS(FEPerfectOsmometer             ,FEMATERIAL_ID, "perfect osmometer"             );
REGISTER_FECORE_CLASS(FERigidMaterial                ,FEMATERIAL_ID, "rigid body"                    );
REGISTER_FECORE_CLASS(FESphericalFiberDistribution   ,FEMATERIAL_ID, "spherical fiber distribution"  );
REGISTER_FECORE_CLASS(FEStVenantKirchhoff            ,FEMATERIAL_ID, "St.Venant-Kirchhoff"           );
REGISTER_FECORE_CLASS(FETCNonlinearOrthotropic       ,FEMATERIAL_ID, "TC nonlinear orthotropic"      );
REGISTER_FECORE_CLASS(FETendonMaterial               ,FEMATERIAL_ID, "tendon material"               );
REGISTER_FECORE_CLASS(FETransIsoMooneyRivlin         ,FEMATERIAL_ID, "trans iso Mooney-Rivlin"       );
REGISTER_FECORE_CLASS(FETransIsoVerondaWestmann      ,FEMATERIAL_ID, "trans iso Veronda-Westmann"    );
REGISTER_FECORE_CLASS(FETrussMaterial                ,FEMATERIAL_ID, "linear truss"                  );
REGISTER_FECORE_CLASS(FEUncoupledElasticMixture      ,FEMATERIAL_ID, "uncoupled solid mixture"       );
REGISTER_FECORE_CLASS(FEVerondaWestmann              ,FEMATERIAL_ID, "Veronda-Westmann"              );
REGISTER_FECORE_CLASS(FEViscoElasticMaterial         ,FEMATERIAL_ID, "viscoelastic"                  );
REGISTER_FECORE_CLASS(FEUncoupledViscoElasticMaterial,FEMATERIAL_ID, "uncoupled viscoelastic"        );
REGISTER_FECORE_CLASS(FEVonMisesPlasticity           ,FEMATERIAL_ID, "von-Mises plasticity"          );
REGISTER_FECORE_CLASS(FEElasticMultigeneration       ,FEMATERIAL_ID, "multigeneration"               );
REGISTER_FECORE_CLASS(FEMRVonMisesFibers             ,FEMATERIAL_ID, "Mooney-Rivlin von Mises Fibers");
REGISTER_FECORE_CLASS(FEUncoupledActiveContraction   ,FEMATERIAL_ID, "uncoupled active contraction"  );
REGISTER_FECORE_CLASS(FEHuiskesSupply                ,FEMATERIAL_ID, "Huiskes-supply"                );
REGISTER_FECORE_CLASS(FERemodelingElasticMaterial    ,FEMATERIAL_ID, "remodeling solid"              );
REGISTER_FECORE_CLASS(FECarterHayesOld               ,FEMATERIAL_ID, "Carter-Hayes (old)"            );
REGISTER_FECORE_CLASS(FEActiveFiberContraction       ,FEMATERIAL_ID, "active_contraction"            );
REGISTER_FECORE_CLASS(FEPreStrainTransIsoMR          ,FEMATERIAL_ID, "pre-strain trans iso Mooney-Rivlin");
REGISTER_FECORE_CLASS(FEPreStrainCoupledTransIsoMR   ,FEMATERIAL_ID, "pre-strain coupled trans iso Mooney-Rivlin");
REGISTER_FECORE_CLASS(FEPreStrainElastic             ,FEMATERIAL_ID, "pre-strain elastic"            );
REGISTER_FECORE_CLASS(FEFiberExponentialPower        ,FEMATERIAL_ID, "fiber-exponential-power-law"   );
REGISTER_FECORE_CLASS(FEFiberNH                      ,FEMATERIAL_ID, "fiber-NH"                      );
REGISTER_FECORE_CLASS(FESphericalFiberDensityDistribution  , FEMATERIAL_ID, "spherical"   );
REGISTER_FECORE_CLASS(FEEllipsodialFiberDensityDistribution, FEMATERIAL_ID, "ellipsoidal" );
REGISTER_FECORE_CLASS(FEVonMises3DFiberDensityDistribution , FEMATERIAL_ID, "von-Mises-3d");
REGISTER_FECORE_CLASS(FECircularFiberDensityDistribution   , FEMATERIAL_ID, "circular"    );
REGISTER_FECORE_CLASS(FEEllipticalFiberDensityDistribution , FEMATERIAL_ID, "elliptical"  );
REGISTER_FECORE_CLASS(FEVonMises2DFiberDensityDistribution , FEMATERIAL_ID, "von-Mises-2d");
REGISTER_FECORE_CLASS(FEContinuousFiberDistribution        , FEMATERIAL_ID, "continuous fiber distribution");
REGISTER_FECORE_CLASS(FEFiberIntegrationGauss              , FEMATERIAL_ID, "fibers-3d-gauss");
REGISTER_FECORE_CLASS(FEFiberIntegrationTrapezoidal        , FEMATERIAL_ID, "fibers-2d-trapezoidal");
REGISTER_FECORE_CLASS(FEFiberIntegrationGeodesic           , FEMATERIAL_ID, "fibers-3d-geodesic");
REGISTER_FECORE_CLASS(FECoupledTransIsoVerondaWestmann     , FEMATERIAL_ID, "coupled trans-iso Veronda-Westmann");
REGISTER_FECORE_CLASS(FECoupledTransIsoMooneyRivlin        , FEMATERIAL_ID, "coupled trans-iso Mooney-Rivlin");
REGISTER_FECORE_CLASS(FEMicroMaterial                      , FEMATERIAL_ID, "micro-material");
REGISTER_FECORE_CLASS(FEGenerationMaterial                 , FEMATERIAL_ID, "generation");
REGISTER_FECORE_CLASS(FEPRLig					           , FEMATERIAL_ID, "PRLig");

REGISTER_FECORE_CLASS(FELinearSpring           , FEMATERIAL_ID, "linear");
REGISTER_FECORE_CLASS(FETensionOnlyLinearSpring, FEMATERIAL_ID, "tension-only linear");
REGISTER_FECORE_CLASS(FENonLinearSpring        , FEMATERIAL_ID, "nonlinear");

//-----------------------------------------------------------------------------
// classes derived from FESurfaceLoad
REGISTER_FECORE_CLASS(FEPressureLoad, FESURFACELOAD_ID, "pressure");
REGISTER_FECORE_CLASS(FETractionLoad, FESURFACELOAD_ID, "traction");

//-----------------------------------------------------------------------------
// classes derived from FEBodyForce
REGISTER_FECORE_CLASS(FEConstBodyForce      , FEBODYLOAD_ID, "const"      );
REGISTER_FECORE_CLASS(FENonConstBodyForce   , FEBODYLOAD_ID, "non-const"  );
REGISTER_FECORE_CLASS(FECentrifugalBodyForce, FEBODYLOAD_ID, "centrifugal");
REGISTER_FECORE_CLASS(FEPointBodyForce      , FEBODYLOAD_ID, "point"      );

//-----------------------------------------------------------------------------
// constraint classes
REGISTER_FECORE_CLASS(FEPointConstraint    , FENLCONSTRAINT_ID, "point"                );
REGISTER_FECORE_CLASS(FEInSituStretch      , FENLCONSTRAINT_ID, "in-situ stretch"      );
REGISTER_FECORE_CLASS(FELinearConstraintSet, FENLCONSTRAINT_ID, "linear constraint"    );
REGISTER_FECORE_CLASS(FERigidJoint         , FENLCONSTRAINT_ID, "rigid joint"          );
REGISTER_FECORE_CLASS(FERigidSphericalJoint, FENLCONSTRAINT_ID, "rigid spherical joint");
REGISTER_FECORE_CLASS(FERigidPinJoint      , FENLCONSTRAINT_ID, "rigid pin joint"      );

//-----------------------------------------------------------------------------
// classes derived from FEContactInterface
REGISTER_FECORE_CLASS(FEFacet2FacetSliding   , FESURFACEPAIRINTERACTION_ID, "facet-to-facet sliding"     );
REGISTER_FECORE_CLASS(FEPeriodicBoundary     , FESURFACEPAIRINTERACTION_ID, "periodic boundary"          );
REGISTER_FECORE_CLASS(FERigidWallInterface   , FESURFACEPAIRINTERACTION_ID, "rigid_wall"                 );
REGISTER_FECORE_CLASS(FESlidingInterface     , FESURFACEPAIRINTERACTION_ID, "sliding_with_gaps"          );
REGISTER_FECORE_CLASS(FESlidingInterfaceBW   , FESURFACEPAIRINTERACTION_ID, "sliding-tension-compression");
REGISTER_FECORE_CLASS(FESurfaceConstraint    , FESURFACEPAIRINTERACTION_ID, "surface constraint"         );
REGISTER_FECORE_CLASS(FETiedInterface        , FESURFACEPAIRINTERACTION_ID, "tied"                       );
REGISTER_FECORE_CLASS(FEStickyInterface      , FESURFACEPAIRINTERACTION_ID, "sticky"                     );
REGISTER_FECORE_CLASS(FEFacet2FacetTied      , FESURFACEPAIRINTERACTION_ID, "facet-to-facet tied"        );

//-----------------------------------------------------------------------------
// classes derived from FEPlotData
REGISTER_FECORE_CLASS(FEPlotStrainEnergyDensity  , FEPLOTDATA_ID, "strain energy density" );
REGISTER_FECORE_CLASS(FEPlotSpecificStrainEnergy , FEPLOTDATA_ID, "specific strain energy");
REGISTER_FECORE_CLASS(FEPlotDensity              , FEPLOTDATA_ID, "density"               );
REGISTER_FECORE_CLASS(FEPlotElementStress        , FEPLOTDATA_ID, "stress"                );
REGISTER_FECORE_CLASS(FEPlotElementElasticity    , FEPLOTDATA_ID, "elasticity"            );
REGISTER_FECORE_CLASS(FEPlotRelativeVolume       , FEPLOTDATA_ID, "relative volume"       );
REGISTER_FECORE_CLASS(FEPlotFiberVector          , FEPLOTDATA_ID, "fiber vector"          );
REGISTER_FECORE_CLASS(FEPlotShellThickness       , FEPLOTDATA_ID, "shell thickness"       );
REGISTER_FECORE_CLASS(FEPlotDamage               , FEPLOTDATA_ID, "damage"                );
REGISTER_FECORE_CLASS(FEPlotMixtureVolumeFraction, FEPLOTDATA_ID, "volume fraction"       );
REGISTER_FECORE_CLASS(FEPlotUT4NodalStresses     , FEPLOTDATA_ID, "ut4 nodal stress"      );
REGISTER_FECORE_CLASS(FEPlotFiberPreStretch		 , FEPLOTDATA_ID, "in-situ fiber stretch" );
REGISTER_FECORE_CLASS(FEPlotShellStrain          , FEPLOTDATA_ID, "shell strain"          );
REGISTER_FECORE_CLASS(FEPlotContactGap			 , FEPLOTDATA_ID, "contact gap"           );
REGISTER_FECORE_CLASS(FEPlotContactPressure		 , FEPLOTDATA_ID, "contact pressure"      );
REGISTER_FECORE_CLASS(FEPlotContactTraction		 , FEPLOTDATA_ID, "contact traction"      );
REGISTER_FECORE_CLASS(FEPlotContactForce 		 , FEPLOTDATA_ID, "contact force"         );
REGISTER_FECORE_CLASS(FEPlotContactArea 		 , FEPLOTDATA_ID, "contact area"          );
REGISTER_FECORE_CLASS(FEPlotSPRStresses          , FEPLOTDATA_ID, "SPR stress"            );
REGISTER_FECORE_CLASS(FEPlotSPRPrincStresses     , FEPLOTDATA_ID, "SPR principal stress"  );
REGISTER_FECORE_CLASS(FEPlotSPRTestLinear		 , FEPLOTDATA_ID, "SPR test linear"       );
REGISTER_FECORE_CLASS(FEPlotSPRTestQuadratic	 , FEPLOTDATA_ID, "SPR test quadratic"    );

//-----------------------------------------------------------------------------
REGISTER_FECORE_CLASS(FEPlotNodeDisplacement   , FEPLOTDATA_ID, "displacement"   );
REGISTER_FECORE_CLASS(FEPlotNodeVelocity       , FEPLOTDATA_ID, "velocity"       );
REGISTER_FECORE_CLASS(FEPlotNodeAcceleration   , FEPLOTDATA_ID, "acceleration"   );
REGISTER_FECORE_CLASS(FEPlotNodeReactionForces , FEPLOTDATA_ID, "reaction forces");
REGISTER_FECORE_CLASS(FEPlotRigidReactionTorque, FEPLOTDATA_ID, "rigid torque"   );
REGISTER_FECORE_CLASS(FEPlotRigidDisplacement       , FEPLOTDATA_ID, "rigid position"            );
REGISTER_FECORE_CLASS(FEPlotRigidVelocity           , FEPLOTDATA_ID, "rigid velocity"            );
REGISTER_FECORE_CLASS(FEPlotRigidAcceleration       , FEPLOTDATA_ID, "rigid acceleration"        );
REGISTER_FECORE_CLASS(FEPlotRigidRotation           , FEPLOTDATA_ID, "rigid angular position"    );
REGISTER_FECORE_CLASS(FEPlotRigidAngularVelocity    , FEPLOTDATA_ID, "rigid angular velocity"    );
REGISTER_FECORE_CLASS(FEPlotRigidAngularAcceleration, FEPLOTDATA_ID, "rigid angular acceleration");
REGISTER_FECORE_CLASS(FEPlotRigidKineticEnergy      , FEPLOTDATA_ID, "rigid kinetic energy"      );

//-----------------------------------------------------------------------------
// Derived from FENodeLogData
REGISTER_FECORE_CLASS(FENodeXPos  , FENODELOGDATA_ID, "x");
REGISTER_FECORE_CLASS(FENodeYPos  , FENODELOGDATA_ID, "y");
REGISTER_FECORE_CLASS(FENodeZPos  , FENODELOGDATA_ID, "z");
REGISTER_FECORE_CLASS(FENodeXDisp , FENODELOGDATA_ID, "ux");
REGISTER_FECORE_CLASS(FENodeYDisp , FENODELOGDATA_ID, "uy");
REGISTER_FECORE_CLASS(FENodeZDisp , FENODELOGDATA_ID, "uz");
REGISTER_FECORE_CLASS(FENodeXVel  , FENODELOGDATA_ID, "vx");
REGISTER_FECORE_CLASS(FENodeYVel  , FENODELOGDATA_ID, "vy");
REGISTER_FECORE_CLASS(FENodeZVel  , FENODELOGDATA_ID, "vz");
REGISTER_FECORE_CLASS(FENodeXAcc  , FENODELOGDATA_ID, "ax");
REGISTER_FECORE_CLASS(FENodeYAcc  , FENODELOGDATA_ID, "ay");
REGISTER_FECORE_CLASS(FENodeZAcc  , FENODELOGDATA_ID, "az");
REGISTER_FECORE_CLASS(FENodeForceX, FENODELOGDATA_ID, "Rx");
REGISTER_FECORE_CLASS(FENodeForceY, FENODELOGDATA_ID, "Ry");
REGISTER_FECORE_CLASS(FENodeForceZ, FENODELOGDATA_ID, "Rz");

//-----------------------------------------------------------------------------
// Derived from FELogElemData
REGISTER_FECORE_CLASS(FELogElemPosX    , FEELEMLOGDATA_ID, "x");
REGISTER_FECORE_CLASS(FELogElemPosY    , FEELEMLOGDATA_ID, "y");
REGISTER_FECORE_CLASS(FELogElemPosZ    , FEELEMLOGDATA_ID, "z");
REGISTER_FECORE_CLASS(FELogElemJacobian, FEELEMLOGDATA_ID, "J");
REGISTER_FECORE_CLASS(FELogElemStrainX , FEELEMLOGDATA_ID, "Ex");
REGISTER_FECORE_CLASS(FELogElemStrainY , FEELEMLOGDATA_ID, "Ey");
REGISTER_FECORE_CLASS(FELogElemStrainZ , FEELEMLOGDATA_ID, "Ez");
REGISTER_FECORE_CLASS(FELogElemStrainXY, FEELEMLOGDATA_ID, "Exy");
REGISTER_FECORE_CLASS(FELogElemStrainYZ, FEELEMLOGDATA_ID, "Eyz");
REGISTER_FECORE_CLASS(FELogElemStrainXZ, FEELEMLOGDATA_ID, "Exz");
REGISTER_FECORE_CLASS(FELogElemStrain1 , FEELEMLOGDATA_ID, "E1");
REGISTER_FECORE_CLASS(FELogElemStrain2 , FEELEMLOGDATA_ID, "E2");
REGISTER_FECORE_CLASS(FELogElemStrain3 , FEELEMLOGDATA_ID, "E3");
REGISTER_FECORE_CLASS(FELogElemStressX , FEELEMLOGDATA_ID, "sx");
REGISTER_FECORE_CLASS(FELogElemStressY , FEELEMLOGDATA_ID, "sy");
REGISTER_FECORE_CLASS(FELogElemStressZ , FEELEMLOGDATA_ID, "sz");
REGISTER_FECORE_CLASS(FELogElemStressXY, FEELEMLOGDATA_ID, "sxy");
REGISTER_FECORE_CLASS(FELogElemStressYZ, FEELEMLOGDATA_ID, "syz");
REGISTER_FECORE_CLASS(FELogElemStressXZ, FEELEMLOGDATA_ID, "sxz");
REGISTER_FECORE_CLASS(FELogElemStress1 , FEELEMLOGDATA_ID, "s1");
REGISTER_FECORE_CLASS(FELogElemStress2 , FEELEMLOGDATA_ID, "s2");
REGISTER_FECORE_CLASS(FELogElemStress3 , FEELEMLOGDATA_ID, "s3");
REGISTER_FECORE_CLASS(FELogElemDeformationGradientXX, FEELEMLOGDATA_ID, "Fxx");
REGISTER_FECORE_CLASS(FELogElemDeformationGradientXY, FEELEMLOGDATA_ID, "Fxy");
REGISTER_FECORE_CLASS(FELogElemDeformationGradientXZ, FEELEMLOGDATA_ID, "Fxz");
REGISTER_FECORE_CLASS(FELogElemDeformationGradientYX, FEELEMLOGDATA_ID, "Fyx");
REGISTER_FECORE_CLASS(FELogElemDeformationGradientYY, FEELEMLOGDATA_ID, "Fyy");
REGISTER_FECORE_CLASS(FELogElemDeformationGradientYZ, FEELEMLOGDATA_ID, "Fyz");
REGISTER_FECORE_CLASS(FELogElemDeformationGradientZX, FEELEMLOGDATA_ID, "Fzx");
REGISTER_FECORE_CLASS(FELogElemDeformationGradientZY, FEELEMLOGDATA_ID, "Fzy");
REGISTER_FECORE_CLASS(FELogElemDeformationGradientZZ, FEELEMLOGDATA_ID, "Fzz");

//-----------------------------------------------------------------------------
// Derived from FELogObjectData
REGISTER_FECORE_CLASS(FELogRigidBodyPosX   , FEOBJLOGDATA_ID, "x");
REGISTER_FECORE_CLASS(FELogRigidBodyPosY   , FEOBJLOGDATA_ID, "y");
REGISTER_FECORE_CLASS(FELogRigidBodyPosZ   , FEOBJLOGDATA_ID, "z");
REGISTER_FECORE_CLASS(FELogRigidBodyVelX   , FEOBJLOGDATA_ID, "vx");
REGISTER_FECORE_CLASS(FELogRigidBodyVelY   , FEOBJLOGDATA_ID, "vy");
REGISTER_FECORE_CLASS(FELogRigidBodyVelZ   , FEOBJLOGDATA_ID, "vz");
REGISTER_FECORE_CLASS(FELogRigidBodyAccX   , FEOBJLOGDATA_ID, "ax");
REGISTER_FECORE_CLASS(FELogRigidBodyAccY   , FEOBJLOGDATA_ID, "ay");
REGISTER_FECORE_CLASS(FELogRigidBodyAccZ   , FEOBJLOGDATA_ID, "az");
REGISTER_FECORE_CLASS(FELogRigidBodyAngPosX, FEOBJLOGDATA_ID, "thx");
REGISTER_FECORE_CLASS(FELogRigidBodyAngPosY, FEOBJLOGDATA_ID, "thy");
REGISTER_FECORE_CLASS(FELogRigidBodyAngPosZ, FEOBJLOGDATA_ID, "thz");
REGISTER_FECORE_CLASS(FELogRigidBodyAngVelX, FEOBJLOGDATA_ID, "omx");
REGISTER_FECORE_CLASS(FELogRigidBodyAngVelY, FEOBJLOGDATA_ID, "omy");
REGISTER_FECORE_CLASS(FELogRigidBodyAngVelZ, FEOBJLOGDATA_ID, "omz");
REGISTER_FECORE_CLASS(FELogRigidBodyAngAccX, FEOBJLOGDATA_ID, "alx");
REGISTER_FECORE_CLASS(FELogRigidBodyAngAccY, FEOBJLOGDATA_ID, "aly");
REGISTER_FECORE_CLASS(FELogRigidBodyAngAccZ, FEOBJLOGDATA_ID, "alz");
REGISTER_FECORE_CLASS(FELogRigidBodyQuatX  , FEOBJLOGDATA_ID, "qx");
REGISTER_FECORE_CLASS(FELogRigidBodyQuatY  , FEOBJLOGDATA_ID, "qy");
REGISTER_FECORE_CLASS(FELogRigidBodyQuatZ  , FEOBJLOGDATA_ID, "qz");
REGISTER_FECORE_CLASS(FELogRigidBodyQuatW  , FEOBJLOGDATA_ID, "qw");
REGISTER_FECORE_CLASS(FELogRigidBodyR11    , FEOBJLOGDATA_ID, "R11");
REGISTER_FECORE_CLASS(FELogRigidBodyR12    , FEOBJLOGDATA_ID, "R12");
REGISTER_FECORE_CLASS(FELogRigidBodyR13    , FEOBJLOGDATA_ID, "R13");
REGISTER_FECORE_CLASS(FELogRigidBodyR21    , FEOBJLOGDATA_ID, "R21");
REGISTER_FECORE_CLASS(FELogRigidBodyR22    , FEOBJLOGDATA_ID, "R22");
REGISTER_FECORE_CLASS(FELogRigidBodyR23    , FEOBJLOGDATA_ID, "R23");
REGISTER_FECORE_CLASS(FELogRigidBodyR31    , FEOBJLOGDATA_ID, "R31");
REGISTER_FECORE_CLASS(FELogRigidBodyR32    , FEOBJLOGDATA_ID, "R32");
REGISTER_FECORE_CLASS(FELogRigidBodyR33    , FEOBJLOGDATA_ID, "R33");
REGISTER_FECORE_CLASS(FELogRigidBodyForceX , FEOBJLOGDATA_ID, "Fx");
REGISTER_FECORE_CLASS(FELogRigidBodyForceY , FEOBJLOGDATA_ID, "Fy");
REGISTER_FECORE_CLASS(FELogRigidBodyForceZ , FEOBJLOGDATA_ID, "Fz");
REGISTER_FECORE_CLASS(FELogRigidBodyTorqueX, FEOBJLOGDATA_ID, "Mx");
REGISTER_FECORE_CLASS(FELogRigidBodyTorqueY, FEOBJLOGDATA_ID, "My");
REGISTER_FECORE_CLASS(FELogRigidBodyTorqueZ, FEOBJLOGDATA_ID, "Mz");
}
