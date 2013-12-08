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
#include "FEDiscreteMaterial.h"
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
#include "FERigid.h"
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
#include "FEElasticFiberMaterial.h"
#include "FEFiberDensityDistribution.h"
#include "FEContinuousFiberDistribution.h"
#include "FEFiberIntegrationGauss.h"
#include "FEFiberIntegrationTrapezoidal.h"
#include "FEFiberIntegrationGeodesic.h"

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
#include "FEInSituStretch.h"
#include "FEPointConstraint.h"
#include "FEFacet2FacetTied.h"

#include "FESolidAnalysis.h"
#include "FESolidSolver.h"
#include "FELinearSolidSolver.h"
#include "FEExplicitSolidSolver.h"

#include "FEBioMechPlot.h"
#include "FEBioMechData.h"

#include "FESolidDomainFactory.h"

//-----------------------------------------------------------------------------
//! Register all the classes of the FEBioMech module with the FEBio framework.
void FEBioMech::InitModule()
{
//-----------------------------------------------------------------------------
// Domain factory
	FEBioKernel& febio = FEBioKernel::GetInstance();
	febio.RegisterDomain(new FESolidDomainFactory);

//-----------------------------------------------------------------------------
// Analysis classes
REGISTER_FEBIO_CLASS(FESolidAnalysis        , FEANALYSIS_ID, "solid"         );
REGISTER_FEBIO_CLASS(FEExplicitSolidAnalysis, FEANALYSIS_ID, "explicit-solid");
REGISTER_FEBIO_CLASS(FELinearSolidAnalysis  , FEANALYSIS_ID, "linear-solid"  );

//-----------------------------------------------------------------------------
// Solver classes
REGISTER_FEBIO_CLASS(FESolidSolver        , FESOLVER_ID, "solid"         );
REGISTER_FEBIO_CLASS(FEExplicitSolidSolver, FESOLVER_ID, "explicit-solid");
REGISTER_FEBIO_CLASS(FELinearSolidSolver  , FESOLVER_ID, "linear-solid"  );

//-----------------------------------------------------------------------------
// material classes
REGISTER_FEBIO_CLASS(FE2DFiberNeoHookean            ,FEMATERIAL_ID, "2D fiber neo-Hookean"          );
REGISTER_FEBIO_CLASS(FE2DTransIsoMooneyRivlin       ,FEMATERIAL_ID, "2D trans iso Mooney-Rivlin"    );
REGISTER_FEBIO_CLASS(FE2DTransIsoVerondaWestmann    ,FEMATERIAL_ID, "2D trans iso Veronda-Westmann" );
REGISTER_FEBIO_CLASS(FEArrudaBoyce                  ,FEMATERIAL_ID, "Arruda-Boyce"                  );
REGISTER_FEBIO_CLASS(FECellGrowth                   ,FEMATERIAL_ID, "cell growth"                   );
REGISTER_FEBIO_CLASS(FEDamageMooneyRivlin           ,FEMATERIAL_ID, "damage Mooney-Rivlin"          );
REGISTER_FEBIO_CLASS(FEDamageNeoHookean             ,FEMATERIAL_ID, "damage neo-Hookean"            );
REGISTER_FEBIO_CLASS(FEDamageTransIsoMooneyRivlin   ,FEMATERIAL_ID, "damage trans iso Mooney-Rivlin");
REGISTER_FEBIO_CLASS(FEDonnanEquilibrium            ,FEMATERIAL_ID, "Donnan equilibrium"            );
REGISTER_FEBIO_CLASS(FEEFD                          ,FEMATERIAL_ID, "EFD"                           );
REGISTER_FEBIO_CLASS(FEEFDDonnanEquilibrium         ,FEMATERIAL_ID, "EFD Donnan equilibrium"        );
REGISTER_FEBIO_CLASS(FEEFDMooneyRivlin              ,FEMATERIAL_ID, "EFD Mooney-Rivlin"             );
//REGISTER_FEBIO_CLASS(FEEFDNeoHookean                ,FEMATERIAL_ID, "EFD neo-Hookean"               );
REGISTER_FEBIO_CLASS(FEEFDNeoHookeanOld             ,FEMATERIAL_ID, "EFD neo-Hookean"               );
REGISTER_FEBIO_CLASS(FEEFDUncoupled                 ,FEMATERIAL_ID, "EFD uncoupled"                 );
REGISTER_FEBIO_CLASS(FEEFDVerondaWestmann           ,FEMATERIAL_ID, "EFD Veronda-Westmann"          );
REGISTER_FEBIO_CLASS(FEElasticMixture               ,FEMATERIAL_ID, "solid mixture"                 );
REGISTER_FEBIO_CLASS(FEEllipsoidalFiberDistribution ,FEMATERIAL_ID, "ellipsoidal fiber distribution");
REGISTER_FEBIO_CLASS(FEEllipsoidalFiberDistributionOld,FEMATERIAL_ID, "ellipsoidal fiber distribution (old)");
REGISTER_FEBIO_CLASS(FEFiberExpPow                  ,FEMATERIAL_ID, "fiber-exp-pow"                 );
REGISTER_FEBIO_CLASS(FEFiberExpPowUncoupled         ,FEMATERIAL_ID, "fiber-exp-pow-uncoupled"       );
REGISTER_FEBIO_CLASS(FEFiberNeoHookean              ,FEMATERIAL_ID, "fiber neo-Hookean"             );
REGISTER_FEBIO_CLASS(FEFungOrthoCompressible        ,FEMATERIAL_ID, "Fung-ortho-compressible"       );
REGISTER_FEBIO_CLASS(FEFungOrthotropic              ,FEMATERIAL_ID, "Fung orthotropic"              );
REGISTER_FEBIO_CLASS(FEGasserOgdenHolzapfel         ,FEMATERIAL_ID, "Gasser-Ogden-Holzapfel"        );
REGISTER_FEBIO_CLASS(FEHolmesMow                    ,FEMATERIAL_ID, "Holmes-Mow"                    );
REGISTER_FEBIO_CLASS(FEIncompNeoHookean             ,FEMATERIAL_ID, "incomp neo-Hookean"            );
REGISTER_FEBIO_CLASS(FEIsotropicElastic             ,FEMATERIAL_ID, "isotropic elastic"             );
REGISTER_FEBIO_CLASS(FELinearElastic                ,FEMATERIAL_ID, "linear elastic"                );
REGISTER_FEBIO_CLASS(FELinearOrthotropic            ,FEMATERIAL_ID, "linear orthotropic"            );
REGISTER_FEBIO_CLASS(FELinearSpring                 ,FEMATERIAL_ID, "linear spring"                 );
REGISTER_FEBIO_CLASS(FELinearTransIso               ,FEMATERIAL_ID, "linear trans iso"              );
REGISTER_FEBIO_CLASS(FEMooneyRivlin                 ,FEMATERIAL_ID, "Mooney-Rivlin"                 );
REGISTER_FEBIO_CLASS(FEMuscleMaterial               ,FEMATERIAL_ID, "muscle material"               );
REGISTER_FEBIO_CLASS(FENeoHookean                   ,FEMATERIAL_ID, "neo-Hookean"                   );
REGISTER_FEBIO_CLASS(FENeoHookeanTransIso           ,FEMATERIAL_ID, "neo-Hookean transiso"          );
REGISTER_FEBIO_CLASS(FENonLinearSpring              ,FEMATERIAL_ID, "nonlinear spring"              );
REGISTER_FEBIO_CLASS(FEOgdenMaterial                ,FEMATERIAL_ID, "Ogden"                         );
REGISTER_FEBIO_CLASS(FEOgdenUnconstrained           ,FEMATERIAL_ID, "Ogden unconstrained"           );
REGISTER_FEBIO_CLASS(FEOrthoElastic                 ,FEMATERIAL_ID, "orthotropic elastic"           );
REGISTER_FEBIO_CLASS(FEPerfectOsmometer             ,FEMATERIAL_ID, "perfect osmometer"             );
REGISTER_FEBIO_CLASS(FERigidMaterial                ,FEMATERIAL_ID, "rigid body"                    );
REGISTER_FEBIO_CLASS(FESphericalFiberDistribution   ,FEMATERIAL_ID, "spherical fiber distribution"  );
REGISTER_FEBIO_CLASS(FEStVenantKirchhoff            ,FEMATERIAL_ID, "St.Venant-Kirchhoff"           );
REGISTER_FEBIO_CLASS(FETCNonlinearOrthotropic       ,FEMATERIAL_ID, "TC nonlinear orthotropic"      );
REGISTER_FEBIO_CLASS(FETendonMaterial               ,FEMATERIAL_ID, "tendon material"               );
REGISTER_FEBIO_CLASS(FETensionOnlyLinearSpring      ,FEMATERIAL_ID, "tension only linear spring"    );
REGISTER_FEBIO_CLASS(FETransIsoMooneyRivlin         ,FEMATERIAL_ID, "trans iso Mooney-Rivlin"       );
REGISTER_FEBIO_CLASS(FETransIsoVerondaWestmann      ,FEMATERIAL_ID, "trans iso Veronda-Westmann"    );
REGISTER_FEBIO_CLASS(FETrussMaterial                ,FEMATERIAL_ID, "linear truss"                  );
REGISTER_FEBIO_CLASS(FEUncoupledElasticMixture      ,FEMATERIAL_ID, "uncoupled solid mixture"       );
REGISTER_FEBIO_CLASS(FEVerondaWestmann              ,FEMATERIAL_ID, "Veronda-Westmann"              );
REGISTER_FEBIO_CLASS(FEViscoElasticMaterial         ,FEMATERIAL_ID, "viscoelastic"                  );
REGISTER_FEBIO_CLASS(FEUncoupledViscoElasticMaterial,FEMATERIAL_ID, "uncoupled viscoelastic"        );
REGISTER_FEBIO_CLASS(FEVonMisesPlasticity           ,FEMATERIAL_ID, "von-Mises plasticity"          );
REGISTER_FEBIO_CLASS(FEElasticMultigeneration       ,FEMATERIAL_ID, "multigeneration"               );
REGISTER_FEBIO_CLASS(FEMRVonMisesFibers             ,FEMATERIAL_ID, "Mooney-Rivlin von Mises Fibers");
REGISTER_FEBIO_CLASS(FEUncoupledActiveContraction   ,FEMATERIAL_ID, "uncoupled active contraction"  );
REGISTER_FEBIO_CLASS(FEHuiskesSupply                ,FEMATERIAL_ID, "Huiskes-supply"                );
REGISTER_FEBIO_CLASS(FERemodelingElasticMaterial    ,FEMATERIAL_ID, "remodeling solid"              );
REGISTER_FEBIO_CLASS(FECarterHayesOld               ,FEMATERIAL_ID, "Carter-Hayes (old)"            );
REGISTER_FEBIO_CLASS(FEActiveFiberContraction       ,FEMATERIAL_ID, "active_contraction"            );
REGISTER_FEBIO_CLASS(FEPreStrainTransIsoMR          ,FEMATERIAL_ID, "pre-strain trans iso Mooney-Rivlin");
REGISTER_FEBIO_CLASS(FEFiberExponentialPower        ,FEMATERIAL_ID, "fiber-exponential-power-law"   );
REGISTER_FEBIO_CLASS(FEFiberNH                      ,FEMATERIAL_ID, "fiber-NH"                      );
REGISTER_FEBIO_CLASS(FESphericalFiberDensityDistribution  , FEMATERIAL_ID, "spherical"   );
REGISTER_FEBIO_CLASS(FEEllipsodialFiberDensityDistribution, FEMATERIAL_ID, "ellipsoidal" );
REGISTER_FEBIO_CLASS(FEVonMises3DFiberDensityDistribution , FEMATERIAL_ID, "von-Mises-3d");
REGISTER_FEBIO_CLASS(FECircularFiberDensityDistribution   , FEMATERIAL_ID, "circular"    );
REGISTER_FEBIO_CLASS(FEEllipticalFiberDensityDistribution , FEMATERIAL_ID, "elliptical"  );
REGISTER_FEBIO_CLASS(FEVonMises2DFiberDensityDistribution , FEMATERIAL_ID, "von-Mises-2d");
REGISTER_FEBIO_CLASS(FEContinuousFiberDistribution        , FEMATERIAL_ID, "continuous fiber distribution");
REGISTER_FEBIO_CLASS(FEFiberIntegrationGauss              , FEMATERIAL_ID, "fibers-3d-gauss");
REGISTER_FEBIO_CLASS(FEFiberIntegrationTrapezoidal        , FEMATERIAL_ID, "fibers-2d-trapezoidal");
REGISTER_FEBIO_CLASS(FEFiberIntegrationGeodesic           , FEMATERIAL_ID, "fibers-3d-geodesic");

//-----------------------------------------------------------------------------
// classes derived from FESurfaceLoad
REGISTER_FEBIO_CLASS(FEPressureLoad, FESURFACELOAD_ID, "pressure");
REGISTER_FEBIO_CLASS(FETractionLoad, FESURFACELOAD_ID, "traction");

//-----------------------------------------------------------------------------
// classes derived from FEBodyForce
REGISTER_FEBIO_CLASS(FEConstBodyForce      , FEBODYLOAD_ID, "const"      );
REGISTER_FEBIO_CLASS(FENonConstBodyForce   , FEBODYLOAD_ID, "non-const"  );
REGISTER_FEBIO_CLASS(FECentrifugalBodyForce, FEBODYLOAD_ID, "centrifugal");
REGISTER_FEBIO_CLASS(FEPointBodyForce      , FEBODYLOAD_ID, "point"      );

//-----------------------------------------------------------------------------
// constraint classes
REGISTER_FEBIO_CLASS(FEPointConstraint, FENLCONSTRAINT_ID, "point"          );
REGISTER_FEBIO_CLASS(FEInSituStretch  , FENLCONSTRAINT_ID, "in-situ stretch");

//-----------------------------------------------------------------------------
// classes derived from FEContactInterface
REGISTER_FEBIO_CLASS(FEFacet2FacetSliding   , FESURFACEPAIRINTERACTION_ID, "facet-to-facet sliding"     );
REGISTER_FEBIO_CLASS(FEPeriodicBoundary     , FESURFACEPAIRINTERACTION_ID, "periodic boundary"          );
REGISTER_FEBIO_CLASS(FERigidWallInterface   , FESURFACEPAIRINTERACTION_ID, "rigid_wall"                 );
REGISTER_FEBIO_CLASS(FESlidingInterface     , FESURFACEPAIRINTERACTION_ID, "sliding_with_gaps"          );
REGISTER_FEBIO_CLASS(FESlidingInterfaceBW   , FESURFACEPAIRINTERACTION_ID, "sliding-tension-compression");
REGISTER_FEBIO_CLASS(FESurfaceConstraint    , FESURFACEPAIRINTERACTION_ID, "surface constraint"         );
REGISTER_FEBIO_CLASS(FETiedInterface        , FESURFACEPAIRINTERACTION_ID, "tied"                       );
REGISTER_FEBIO_CLASS(FEFacet2FacetTied      , FESURFACEPAIRINTERACTION_ID, "facet-to-facet tied"        );

//-----------------------------------------------------------------------------
// classes derived from FEPlotData
REGISTER_FEBIO_CLASS(FEPlotStrainEnergyDensity          , FEPLOTDATA_ID, "strain energy density"           );
REGISTER_FEBIO_CLASS(FEPlotSpecificStrainEnergy         , FEPLOTDATA_ID, "specific strain energy"          );
REGISTER_FEBIO_CLASS(FEPlotDensity                      , FEPLOTDATA_ID, "density"                         );
REGISTER_FEBIO_CLASS(FEPlotElementStress                , FEPLOTDATA_ID, "stress"                          );
REGISTER_FEBIO_CLASS(FEPlotRelativeVolume               , FEPLOTDATA_ID, "relative volume"                 );
REGISTER_FEBIO_CLASS(FEPlotFiberVector                  , FEPLOTDATA_ID, "fiber vector"                    );
REGISTER_FEBIO_CLASS(FEPlotShellThickness               , FEPLOTDATA_ID, "shell thickness"                 );
REGISTER_FEBIO_CLASS(FEPlotDamage                       , FEPLOTDATA_ID, "damage"                          );
REGISTER_FEBIO_CLASS(FEPlotMixtureVolumeFraction        , FEPLOTDATA_ID, "volume fraction"                 );
REGISTER_FEBIO_CLASS(FEPlotUT4NodalStresses             , FEPLOTDATA_ID, "ut4 nodal stress"                );
REGISTER_FEBIO_CLASS(FEPlotFiberPreStretch				, FEPLOTDATA_ID, "in-situ fiber stretch"           );

//-----------------------------------------------------------------------------
REGISTER_FEBIO_CLASS(FEPlotNodeDisplacement  , FEPLOTDATA_ID, "displacement"   );
REGISTER_FEBIO_CLASS(FEPlotNodeVelocity      , FEPLOTDATA_ID, "velocity"       );
REGISTER_FEBIO_CLASS(FEPlotNodeAcceleration  , FEPLOTDATA_ID, "acceleration"   );
REGISTER_FEBIO_CLASS(FEPlotNodeReactionForces, FEPLOTDATA_ID, "reaction forces");

//-----------------------------------------------------------------------------
// Derived from FENodeLogData
REGISTER_FEBIO_CLASS(FENodeXPos  , FENODELOGDATA_ID, "x");
REGISTER_FEBIO_CLASS(FENodeYPos  , FENODELOGDATA_ID, "y");
REGISTER_FEBIO_CLASS(FENodeZPos  , FENODELOGDATA_ID, "z");
REGISTER_FEBIO_CLASS(FENodeXDisp , FENODELOGDATA_ID, "ux");
REGISTER_FEBIO_CLASS(FENodeYDisp , FENODELOGDATA_ID, "uy");
REGISTER_FEBIO_CLASS(FENodeZDisp , FENODELOGDATA_ID, "uz");
REGISTER_FEBIO_CLASS(FENodeXVel  , FENODELOGDATA_ID, "vx");
REGISTER_FEBIO_CLASS(FENodeYVel  , FENODELOGDATA_ID, "vy");
REGISTER_FEBIO_CLASS(FENodeZVel  , FENODELOGDATA_ID, "vz");
REGISTER_FEBIO_CLASS(FENodeForceX, FENODELOGDATA_ID, "Rx");
REGISTER_FEBIO_CLASS(FENodeForceY, FENODELOGDATA_ID, "Ry");
REGISTER_FEBIO_CLASS(FENodeForceZ, FENODELOGDATA_ID, "Rz");

//-----------------------------------------------------------------------------
// Derived from FELogElemData
REGISTER_FEBIO_CLASS(FELogElemPosX    , FEELEMLOGDATA_ID, "x");
REGISTER_FEBIO_CLASS(FELogElemPosY    , FEELEMLOGDATA_ID, "y");
REGISTER_FEBIO_CLASS(FELogElemPosZ    , FEELEMLOGDATA_ID, "z");
REGISTER_FEBIO_CLASS(FELogElemJacobian, FEELEMLOGDATA_ID, "J");
REGISTER_FEBIO_CLASS(FELogElemStrainX , FEELEMLOGDATA_ID, "Ex");
REGISTER_FEBIO_CLASS(FELogElemStrainY , FEELEMLOGDATA_ID, "Ey");
REGISTER_FEBIO_CLASS(FELogElemStrainZ , FEELEMLOGDATA_ID, "Ez");
REGISTER_FEBIO_CLASS(FELogElemStrainXY, FEELEMLOGDATA_ID, "Exy");
REGISTER_FEBIO_CLASS(FELogElemStrainYZ, FEELEMLOGDATA_ID, "Eyz");
REGISTER_FEBIO_CLASS(FELogElemStrainXZ, FEELEMLOGDATA_ID, "Exz");
REGISTER_FEBIO_CLASS(FELogElemStrain1 , FEELEMLOGDATA_ID, "E1");
REGISTER_FEBIO_CLASS(FELogElemStrain2 , FEELEMLOGDATA_ID, "E2");
REGISTER_FEBIO_CLASS(FELogElemStrain3 , FEELEMLOGDATA_ID, "E3");
REGISTER_FEBIO_CLASS(FELogElemStressX , FEELEMLOGDATA_ID, "sx");
REGISTER_FEBIO_CLASS(FELogElemStressY , FEELEMLOGDATA_ID, "sy");
REGISTER_FEBIO_CLASS(FELogElemStressZ , FEELEMLOGDATA_ID, "sz");
REGISTER_FEBIO_CLASS(FELogElemStressXY, FEELEMLOGDATA_ID, "sxy");
REGISTER_FEBIO_CLASS(FELogElemStressYZ, FEELEMLOGDATA_ID, "syz");
REGISTER_FEBIO_CLASS(FELogElemStressXZ, FEELEMLOGDATA_ID, "sxz");
REGISTER_FEBIO_CLASS(FELogElemStress1 , FEELEMLOGDATA_ID, "s1");
REGISTER_FEBIO_CLASS(FELogElemStress2 , FEELEMLOGDATA_ID, "s2");
REGISTER_FEBIO_CLASS(FELogElemStress3 , FEELEMLOGDATA_ID, "s3");
REGISTER_FEBIO_CLASS(FELogElemDeformationGradientXX, FEELEMLOGDATA_ID, "Fxx");
REGISTER_FEBIO_CLASS(FELogElemDeformationGradientXY, FEELEMLOGDATA_ID, "Fxy");
REGISTER_FEBIO_CLASS(FELogElemDeformationGradientXZ, FEELEMLOGDATA_ID, "Fxz");
REGISTER_FEBIO_CLASS(FELogElemDeformationGradientYX, FEELEMLOGDATA_ID, "Fyx");
REGISTER_FEBIO_CLASS(FELogElemDeformationGradientYY, FEELEMLOGDATA_ID, "Fyy");
REGISTER_FEBIO_CLASS(FELogElemDeformationGradientYZ, FEELEMLOGDATA_ID, "Fyz");
REGISTER_FEBIO_CLASS(FELogElemDeformationGradientZX, FEELEMLOGDATA_ID, "Fzx");
REGISTER_FEBIO_CLASS(FELogElemDeformationGradientZY, FEELEMLOGDATA_ID, "Fzy");
REGISTER_FEBIO_CLASS(FELogElemDeformationGradientZZ, FEELEMLOGDATA_ID, "Fzz");

//-----------------------------------------------------------------------------
// Derived from FELogObjectData
REGISTER_FEBIO_CLASS(FELogRigidBodyPosX   , FEOBJLOGDATA_ID, "x");
REGISTER_FEBIO_CLASS(FELogRigidBodyPosY   , FEOBJLOGDATA_ID, "y");
REGISTER_FEBIO_CLASS(FELogRigidBodyPosZ   , FEOBJLOGDATA_ID, "z");
REGISTER_FEBIO_CLASS(FELogRigidBodyQuatX  , FEOBJLOGDATA_ID, "qx");
REGISTER_FEBIO_CLASS(FELogRigidBodyQuatY  , FEOBJLOGDATA_ID, "qy");
REGISTER_FEBIO_CLASS(FELogRigidBodyQuatZ  , FEOBJLOGDATA_ID, "qz");
REGISTER_FEBIO_CLASS(FELogRigidBodyQuatW  , FEOBJLOGDATA_ID, "qw");
REGISTER_FEBIO_CLASS(FELogRigidBodyForceX , FEOBJLOGDATA_ID, "Fx");
REGISTER_FEBIO_CLASS(FELogRigidBodyForceY , FEOBJLOGDATA_ID, "Fy");
REGISTER_FEBIO_CLASS(FELogRigidBodyForceZ , FEOBJLOGDATA_ID, "Fz");
REGISTER_FEBIO_CLASS(FELogRigidBodyTorqueX, FEOBJLOGDATA_ID, "Mx");
REGISTER_FEBIO_CLASS(FELogRigidBodyTorqueY, FEOBJLOGDATA_ID, "My");
REGISTER_FEBIO_CLASS(FELogRigidBodyTorqueZ, FEOBJLOGDATA_ID, "Mz");
}
