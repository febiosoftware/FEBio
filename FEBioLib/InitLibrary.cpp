#include "stdafx.h"
#include "FECore/febio.h"

#include "FEConstBodyForce.h"
#include "FEPointBodyForce.h"

#include "FE2DFiberNeoHookean.h"
#include "FE2DTransIsoMooneyRivlin.h"
#include "FE2DTransIsoVerondaWestmann.h"
#include "FEArrudaBoyce.h"
#include "FEBiphasic.h"
#include "FEBiphasicSolute.h"
#include "FECellGrowth.h"
#include "FEDamageMooneyRivlin.h"
#include "FEDamageNeoHookean.h"
#include "FEDamageTransIsoMooneyRivlin.h"
#include "FEDiffConstIso.h"
#include "FEDiffConstOrtho.h"
#include "FEDiffRefIso.h"
#include "FEDiscreteMaterial.h"
#include "FEDonnanEquilibrium.h"
#include "FEEFDDonnanEquilibrium.h"
#include "FEEFDMooneyRivlin.h"
#include "FEEFDNeoHookean.h"
#include "FEEFDUncoupled.h"
#include "FEEFDVerondaWestmann.h"
#include "FEElasticMixture.h"
#include "FEEllipsoidalFiberDistribution.h"
#include "FEFiberExpPow.h"
#include "FEFiberExpPowUncoupled.h"
#include "FEFiberNeoHookean.h"
#include "FEFungOrthoCompressible.h"
#include "FEFungOrthotropic.h"
#include "FEGasserOgdenHolzapfel.h"
#include "FEHolmesMow.h"
#include "FEIncompNeoHookean.h"
#include "FEIsotropicElastic.h"
#include "FEIsotropicFourier.h"
#include "FELinearElastic.h"
#include "FELinearOrthotropic.h"
#include "FELinearTransIso.h"
#include "FEMooneyRivlin.h"
#include "FEMultiphasic.h"
#include "FEMuscleMaterial.h"
#include "FENeoHookean.h"
#include "FENeoHookeanTransIso.h"
#include "FEOgdenMaterial.h"
#include "FEOgdenUnconstrained.h"
#include "FEOrthoElastic.h"
#include "FEOsmCoefConst.h"
#include "FEPerfectOsmometer.h"
#include "FEPermConstIso.h"
#include "FEPermHolmesMow.h"
#include "FEPermRefIso.h"
#include "FEPermRefOrtho.h"
#include "FEPermRefTransIso.h"
#include "FERigid.h"
#include "FEStVenantKirchhoff.h"
#include "FESolubConst.h"
#include "FESupplyBinding.h"
#include "FESupplyConst.h"
#include "FESupplySynthesisBinding.h"
#include "FESolute.h"
#include "FESphericalFiberDistribution.h"
#include "FETCNonlinearOrthotropic.h"
#include "FETendonMaterial.h"
#include "FETransIsoMooneyRivlin.h"
#include "FETransIsoVerondaWestmann.h"
#include "FETriphasic.h"
#include "FETrussMaterial.h"
#include "FEUncoupledElasticMixture.h"
#include "FEVerondaWestmann.h"
#include "FEViscoElasticMaterial.h"
#include "FEVonMisesPlasticity.h"

#include "FESurfaceConstraint.h"
#include "FEPeriodicBoundary.h"
#include "FETiedInterface.h"
#include "FETiedBiphasicInterface.h"

#include "FEBioPlot/FEPlotDomainData.h"
#include "FEDamageNeoHookean.h"
#include "FEDamageTransIsoMooneyRivlin.h"
#include "FEBiphasicSoluteDomain.h"
#include "FEBiphasicSolidDomain.h"
#include "FETriphasicDomain.h"
#include "FERigidSolidDomain.h"
#include "FERigidShellDomain.h"
#include "FEElasticMixture.h"
#include "FEBiphasicSolute.h"
#include "FETriphasic.h"
#include "FEBioPlot/FEPlotNodeData.h"
#include "FEBioPlot/FEPlotSurfaceData.h"

void InitFEBioLibrary()
{
//-----------------------------------------------------------------------------
// classes derived from FEBodyForce
REGISTER_FEBIO_CLASS(FEConstBodyForce      , FEBodyForce, "const"      );
REGISTER_FEBIO_CLASS(FENonConstBodyForce   , FEBodyForce, "non-const"  );
REGISTER_FEBIO_CLASS(FECentrifugalBodyForce, FEBodyForce, "centrifugal");
REGISTER_FEBIO_CLASS(FEPointBodyForce      , FEBodyForce, "point"      );

//-----------------------------------------------------------------------------
// material classes
REGISTER_MATERIAL(FE2DFiberNeoHookean           , "2D fiber neo-Hookean"          );
REGISTER_MATERIAL(FE2DTransIsoMooneyRivlin      , "2D trans iso Mooney-Rivlin"    );
REGISTER_MATERIAL(FE2DTransIsoVerondaWestmann   , "2D trans iso Veronda-Westmann" );
REGISTER_MATERIAL(FEArrudaBoyce                 , "Arruda-Boyce"                  );
REGISTER_MATERIAL(FEBiphasic                    , "biphasic"                      );
REGISTER_MATERIAL(FEBiphasicSolute              , "biphasic-solute"               );
REGISTER_MATERIAL(FECellGrowth                  , "cell growth"                   );
REGISTER_MATERIAL(FEDamageMooneyRivlin          , "damage Mooney-Rivlin"          );
REGISTER_MATERIAL(FEDamageNeoHookean            , "damage neo-Hookean"            );
REGISTER_MATERIAL(FEDamageTransIsoMooneyRivlin  , "damage trans iso Mooney-Rivlin");
REGISTER_MATERIAL(FEDiffConstIso                , "diff-const-iso"                );
REGISTER_MATERIAL(FEDiffConstOrtho              , "diff-const-ortho"              );
REGISTER_MATERIAL(FEDiffRefIso                  , "diff-ref-iso"                  );
REGISTER_MATERIAL(FEDonnanEquilibrium           , "Donnan equilibrium"            );
REGISTER_MATERIAL(FEEFDDonnanEquilibrium        , "EFD Donnan equilibrium"        );
REGISTER_MATERIAL(FEEFDMooneyRivlin             , "EFD Mooney-Rivlin"             );
REGISTER_MATERIAL(FEEFDNeoHookean               , "EFD neo-Hookean"               );
REGISTER_MATERIAL(FEEFDUncoupled                , "EFD uncoupled"                 );
REGISTER_MATERIAL(FEEFDVerondaWestmann          , "EFD Veronda-Westmann"          );
REGISTER_MATERIAL(FEElasticMixture              , "solid mixture"                 );
REGISTER_MATERIAL(FEEllipsoidalFiberDistribution, "ellipsoidal fiber distribution");
REGISTER_MATERIAL(FEFiberExpPow                 , "fiber-exp-pow"                 );
REGISTER_MATERIAL(FEFiberExpPowUncoupled        , "fiber-exp-pow-uncoupled"       );
REGISTER_MATERIAL(FEFiberNeoHookean             , "fiber neo-Hookean"             );
REGISTER_MATERIAL(FEFungOrthoCompressible       , "Fung-ortho-compressible"       );
REGISTER_MATERIAL(FEFungOrthotropic             , "Fung orthotropic"              );
REGISTER_MATERIAL(FEGasserOgdenHolzapfel        , "Gasser-Ogden-Holzapfel"        );
REGISTER_MATERIAL(FEHolmesMow                   , "Holmes-Mow"                    );
REGISTER_MATERIAL(FEIncompNeoHookean            , "incomp neo-Hookean"            );
REGISTER_MATERIAL(FEIsotropicElastic            , "isotropic elastic"             );
REGISTER_MATERIAL(FEIsotropicFourier            , "isotropic Fourier"             );
REGISTER_MATERIAL(FELinearElastic               , "linear elastic"                );
REGISTER_MATERIAL(FELinearOrthotropic           , "linear orthotropic"            );
REGISTER_MATERIAL(FELinearSpring                , "linear spring"                 );
REGISTER_MATERIAL(FELinearTransIso              , "linear trans iso"              );
REGISTER_MATERIAL(FEMooneyRivlin                , "Mooney-Rivlin"                 );
REGISTER_MATERIAL(FEMultiphasic                 , "multiphasic"                   );
REGISTER_MATERIAL(FEMuscleMaterial              , "muscle material"               );
REGISTER_MATERIAL(FENeoHookean                  , "neo-Hookean"                   );
REGISTER_MATERIAL(FENeoHookeanTransIso          , "neo-Hookean transiso"          );
REGISTER_MATERIAL(FENonLinearSpring             , "nonlinear spring"              );
REGISTER_MATERIAL(FEOgdenMaterial               , "Ogden"                         );
REGISTER_MATERIAL(FEOgdenUnconstrained          , "Ogden unconstrained"           );
REGISTER_MATERIAL(FEOrthoElastic                , "orthotropic elastic"           );
REGISTER_MATERIAL(FEOsmCoefConst                , "osm-coef-const"                );
REGISTER_MATERIAL(FEPerfectOsmometer            , "perfect osmometer"             );
REGISTER_MATERIAL(FEPermConstIso                , "perm-const-iso"                );
REGISTER_MATERIAL(FEPermHolmesMow               , "perm-Holmes-Mow"               );
REGISTER_MATERIAL(FEPermRefIso                  , "perm-ref-iso"                  );
REGISTER_MATERIAL(FEPermRefOrtho                , "perm-ref-ortho"                );
REGISTER_MATERIAL(FEPermRefTransIso             , "perm-ref-trans-iso"            );
REGISTER_MATERIAL(FERigidMaterial               , "rigid body"                    );
REGISTER_MATERIAL(FESphericalFiberDistribution  , "spherical fiber distribution"  );
REGISTER_MATERIAL(FEStVenantKirchhoff           , "St.Venant-Kirchhoff"           );
REGISTER_MATERIAL(FESolubConst                  , "solub-const"                   );
REGISTER_MATERIAL(FESolute                      , "solute"                        );
REGISTER_MATERIAL(FESupplyBinding               , "supply-binding"                );
REGISTER_MATERIAL(FESupplyConst                 , "supply-const"                  );
REGISTER_MATERIAL(FESupplySynthesisBinding      , "supply-synthesis-binding"      );
REGISTER_MATERIAL(FETCNonlinearOrthotropic      , "TC nonlinear orthotropic"      );
REGISTER_MATERIAL(FETendonMaterial              , "tendon material"               );
REGISTER_MATERIAL(FETensionOnlyLinearSpring     , "tension only linear spring"    );
REGISTER_MATERIAL(FETransIsoMooneyRivlin        , "trans iso Mooney-Rivlin"       );
REGISTER_MATERIAL(FETransIsoVerondaWestmann     , "trans iso Veronda-Westmann"    );
REGISTER_MATERIAL(FETriphasic                   , "triphasic"                     );
REGISTER_MATERIAL(FETrussMaterial               , "linear truss"                  );
REGISTER_MATERIAL(FEUncoupledElasticMixture     , "uncoupled solid mixture"       );
REGISTER_MATERIAL(FEVerondaWestmann             , "Veronda-Westmann"              );
REGISTER_MATERIAL(FEViscoElasticMaterial        , "viscoelastic"                  );
REGISTER_MATERIAL(FEVonMisesPlasticity          , "von-Mises plasticity"          );

//-----------------------------------------------------------------------------
// classes derived from FEContactInterface
REGISTER_FEBIO_CLASS(FEPeriodicBoundary     , FEContactInterface, "periodic boundary" );
REGISTER_FEBIO_CLASS(FESurfaceConstraint    , FEContactInterface, "surface constraint");
REGISTER_FEBIO_CLASS(FETiedInterface        , FEContactInterface, "tied"              );
REGISTER_FEBIO_CLASS(FETiedBiphasicInterface, FEContactInterface, "tied-biphasic"     );

//-----------------------------------------------------------------------------
// classes derived from FEPlotData
REGISTER_FEBIO_CLASS(FEPlotEffectiveFluidPressure		, FEPlotData, "effective fluid pressure"        );
REGISTER_FEBIO_CLASS(FEPlotActualFluidPressure          , FEPlotData, "fluid pressure"                  );
REGISTER_FEBIO_CLASS(FEPlotElementStress                , FEPlotData, "stress"                          );
REGISTER_FEBIO_CLASS(FEPlotRelativeVolume               , FEPlotData, "relative volume"                 );
REGISTER_FEBIO_CLASS(FEPlotFluidFlux                    , FEPlotData, "fluid flux"                      );
REGISTER_FEBIO_CLASS(FEPlotFiberVector                  , FEPlotData, "fiber vector"                    );
REGISTER_FEBIO_CLASS(FEPlotEffectiveSoluteConcentration , FEPlotData, "effective solute concentration"  );
REGISTER_FEBIO_CLASS(FEPlotShellThickness               , FEPlotData, "shell thickness"                 );
REGISTER_FEBIO_CLASS(FEPlotActualSoluteConcentration    , FEPlotData, "solute concentration"            );
REGISTER_FEBIO_CLASS(FEPlotSoluteFlux                   , FEPlotData, "solute flux"                     );
REGISTER_FEBIO_CLASS(FEPlotDamage                       , FEPlotData, "damage"                          );
REGISTER_FEBIO_CLASS(FEPlotMixtureVolumeFraction        , FEPlotData, "volume fraction"                 );
REGISTER_FEBIO_CLASS(FEPlotReceptorLigandConcentration  , FEPlotData, "receptor-ligand concentration"   );
REGISTER_FEBIO_CLASS_T(FEPlotEffectiveSolConcentrationT , FEPlotData, 0, "effective solute 1 concentration");
REGISTER_FEBIO_CLASS_T(FEPlotActualSolConcentrationT    , FEPlotData, 0, "solute 1 concentration"          );
REGISTER_FEBIO_CLASS_T(FEPlotSolFluxT                   , FEPlotData, 0, "solute 1 flux"                   );
REGISTER_FEBIO_CLASS_T(FEPlotEffectiveSolConcentrationT , FEPlotData, 1, "effective solute 2 concentration");
REGISTER_FEBIO_CLASS_T(FEPlotActualSolConcentrationT    , FEPlotData, 1, "solute 2 concentration"          );
REGISTER_FEBIO_CLASS_T(FEPlotSolFluxT                   , FEPlotData, 1, "solute 2 flux"                   );
REGISTER_FEBIO_CLASS_T(FEPlotEffectiveSolConcentrationT , FEPlotData, 2, "effective solute 3 concentration");
REGISTER_FEBIO_CLASS_T(FEPlotActualSolConcentrationT    , FEPlotData, 2, "solute 3 concentration"          );
REGISTER_FEBIO_CLASS_T(FEPlotSolFluxT                   , FEPlotData, 2, "solute 3 flux"                   );
REGISTER_FEBIO_CLASS_T(FEPlotEffectiveSolConcentrationT , FEPlotData, 3, "effective solute 4 concentration");
REGISTER_FEBIO_CLASS_T(FEPlotActualSolConcentrationT    , FEPlotData, 3, "solute 4 concentration"          );
REGISTER_FEBIO_CLASS_T(FEPlotSolFluxT                   , FEPlotData, 3, "solute 4 flux"                   );
REGISTER_FEBIO_CLASS_T(FEPlotEffectiveSolConcentrationT , FEPlotData, 4, "effective solute 5 concentration");
REGISTER_FEBIO_CLASS_T(FEPlotActualSolConcentrationT    , FEPlotData, 4, "solute 5 concentration"          );
REGISTER_FEBIO_CLASS_T(FEPlotSolFluxT                   , FEPlotData, 4, "solute 5 flux"                   );
REGISTER_FEBIO_CLASS_T(FEPlotEffectiveSolConcentrationT , FEPlotData, 5, "effective solute 6 concentration");
REGISTER_FEBIO_CLASS_T(FEPlotActualSolConcentrationT    , FEPlotData, 5, "solute 6 concentration"          );
REGISTER_FEBIO_CLASS_T(FEPlotSolFluxT                   , FEPlotData, 5, "solute 6 flux"                   );
REGISTER_FEBIO_CLASS(FEPlotElectricPotential            , FEPlotData, "electric potential"              );
REGISTER_FEBIO_CLASS(FEPlotCurrentDensity               , FEPlotData, "current density"                 );
REGISTER_FEBIO_CLASS(FEPlotFixedChargeDensity           , FEPlotData, "fixed charge density"            );
REGISTER_FEBIO_CLASS(FEPlotNodalFluidFlux               , FEPlotData, "nodal fluid flux"                );

//-----------------------------------------------------------------------------
REGISTER_FEBIO_CLASS(FEPlotNodeDisplacement  , FEPlotData, "displacement"   );
REGISTER_FEBIO_CLASS(FEPlotNodeVelocity      , FEPlotData, "velocity"       );
REGISTER_FEBIO_CLASS(FEPlotNodeAcceleration  , FEPlotData, "acceleration"   );
REGISTER_FEBIO_CLASS(FEPlotNodeTemperature   , FEPlotData, "temperature"    );
REGISTER_FEBIO_CLASS(FEPlotNodeReactionForces, FEPlotData, "reaction forces");

//-----------------------------------------------------------------------------
REGISTER_FEBIO_CLASS(FEPlotContactGap      , FEPlotData, "contact gap"     );
REGISTER_FEBIO_CLASS(FEPlotContactPressure , FEPlotData, "contact pressure");
REGISTER_FEBIO_CLASS(FEPlotContactTraction , FEPlotData, "contact traction");
}
