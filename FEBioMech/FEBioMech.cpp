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

#include "FEConstBodyForce.h"
#include "FEPointBodyForce.h"

#include "FEFacet2FacetSliding.h"
#include "FEPeriodicBoundary.h"
#include "FERigidWallInterface.h"
#include "FESlidingInterface.h"
#include "FESlidingInterfaceBW.h"
#include "FESurfaceConstraint.h"
#include "FETiedInterface.h"

//-----------------------------------------------------------------------------
//! Register all the classes of the FEBioMech module with the FEBio framework.
void FEBioMech::InitModule()
{
//-----------------------------------------------------------------------------
// material classes
REGISTER_MATERIAL(FE2DFiberNeoHookean            , "2D fiber neo-Hookean"          );
REGISTER_MATERIAL(FE2DTransIsoMooneyRivlin       , "2D trans iso Mooney-Rivlin"    );
REGISTER_MATERIAL(FE2DTransIsoVerondaWestmann    , "2D trans iso Veronda-Westmann" );
REGISTER_MATERIAL(FEArrudaBoyce                  , "Arruda-Boyce"                  );
REGISTER_MATERIAL(FECellGrowth                   , "cell growth"                   );
REGISTER_MATERIAL(FEDamageMooneyRivlin           , "damage Mooney-Rivlin"          );
REGISTER_MATERIAL(FEDamageNeoHookean             , "damage neo-Hookean"            );
REGISTER_MATERIAL(FEDamageTransIsoMooneyRivlin   , "damage trans iso Mooney-Rivlin");
REGISTER_MATERIAL(FEDonnanEquilibrium            , "Donnan equilibrium"            );
REGISTER_MATERIAL(FEEFD                          , "EFD"                           );
REGISTER_MATERIAL(FEEFDDonnanEquilibrium         , "EFD Donnan equilibrium"        );
REGISTER_MATERIAL(FEEFDMooneyRivlin              , "EFD Mooney-Rivlin"             );
//REGISTER_MATERIAL(FEEFDNeoHookean                , "EFD neo-Hookean"               );
REGISTER_MATERIAL(FEEFDNeoHookeanOld             , "EFD neo-Hookean"               );
REGISTER_MATERIAL(FEEFDUncoupled                 , "EFD uncoupled"                 );
REGISTER_MATERIAL(FEEFDVerondaWestmann           , "EFD Veronda-Westmann"          );
REGISTER_MATERIAL(FEElasticMixture               , "solid mixture"                 );
REGISTER_MATERIAL(FEEllipsoidalFiberDistribution , "ellipsoidal fiber distribution");
REGISTER_MATERIAL(FEFiberExpPow                  , "fiber-exp-pow"                 );
REGISTER_MATERIAL(FEFiberExpPowUncoupled         , "fiber-exp-pow-uncoupled"       );
REGISTER_MATERIAL(FEFiberNeoHookean              , "fiber neo-Hookean"             );
REGISTER_MATERIAL(FEFungOrthoCompressible        , "Fung-ortho-compressible"       );
REGISTER_MATERIAL(FEFungOrthotropic              , "Fung orthotropic"              );
REGISTER_MATERIAL(FEGasserOgdenHolzapfel         , "Gasser-Ogden-Holzapfel"        );
REGISTER_MATERIAL(FEHolmesMow                    , "Holmes-Mow"                    );
REGISTER_MATERIAL(FEIncompNeoHookean             , "incomp neo-Hookean"            );
REGISTER_MATERIAL(FEIsotropicElastic             , "isotropic elastic"             );
REGISTER_MATERIAL(FELinearElastic                , "linear elastic"                );
REGISTER_MATERIAL(FELinearOrthotropic            , "linear orthotropic"            );
REGISTER_MATERIAL(FELinearSpring                 , "linear spring"                 );
REGISTER_MATERIAL(FELinearTransIso               , "linear trans iso"              );
REGISTER_MATERIAL(FEMooneyRivlin                 , "Mooney-Rivlin"                 );
REGISTER_MATERIAL(FEMuscleMaterial               , "muscle material"               );
REGISTER_MATERIAL(FENeoHookean                   , "neo-Hookean"                   );
REGISTER_MATERIAL(FENeoHookeanTransIso           , "neo-Hookean transiso"          );
REGISTER_MATERIAL(FENonLinearSpring              , "nonlinear spring"              );
REGISTER_MATERIAL(FEOgdenMaterial                , "Ogden"                         );
REGISTER_MATERIAL(FEOgdenUnconstrained           , "Ogden unconstrained"           );
REGISTER_MATERIAL(FEOrthoElastic                 , "orthotropic elastic"           );
REGISTER_MATERIAL(FEPerfectOsmometer             , "perfect osmometer"             );
REGISTER_MATERIAL(FERigidMaterial                , "rigid body"                    );
REGISTER_MATERIAL(FESphericalFiberDistribution   , "spherical fiber distribution"  );
REGISTER_MATERIAL(FEStVenantKirchhoff            , "St.Venant-Kirchhoff"           );
REGISTER_MATERIAL(FETCNonlinearOrthotropic       , "TC nonlinear orthotropic"      );
REGISTER_MATERIAL(FETendonMaterial               , "tendon material"               );
REGISTER_MATERIAL(FETensionOnlyLinearSpring      , "tension only linear spring"    );
REGISTER_MATERIAL(FETransIsoMooneyRivlin         , "trans iso Mooney-Rivlin"       );
REGISTER_MATERIAL(FETransIsoVerondaWestmann      , "trans iso Veronda-Westmann"    );
REGISTER_MATERIAL(FETrussMaterial                , "linear truss"                  );
REGISTER_MATERIAL(FEUncoupledElasticMixture      , "uncoupled solid mixture"       );
REGISTER_MATERIAL(FEVerondaWestmann              , "Veronda-Westmann"              );
REGISTER_MATERIAL(FEViscoElasticMaterial         , "viscoelastic"                  );
REGISTER_MATERIAL(FEUncoupledViscoElasticMaterial, "uncoupled viscoelastic"        );
REGISTER_MATERIAL(FEVonMisesPlasticity           , "von-Mises plasticity"          );
REGISTER_MATERIAL(FEElasticMultigeneration       , "multigeneration"               );
REGISTER_MATERIAL(FEMRVonMisesFibers             , "Mooney-Rivlin von Mises Fibers");
REGISTER_MATERIAL(FEUncoupledActiveContraction   , "uncoupled active contraction"  );
REGISTER_MATERIAL(FEHuiskesSupply                , "Huiskes-supply"                );
REGISTER_MATERIAL(FERemodelingElasticMaterial    , "remodeling solid"              );
REGISTER_MATERIAL(FECarterHayesOld               , "Carter-Hayes (old)"            );

//-----------------------------------------------------------------------------
// classes derived from FEBodyForce
REGISTER_FEBIO_CLASS(FEConstBodyForce      , FEBodyForce, "const"      );
REGISTER_FEBIO_CLASS(FENonConstBodyForce   , FEBodyForce, "non-const"  );
REGISTER_FEBIO_CLASS(FECentrifugalBodyForce, FEBodyForce, "centrifugal");
REGISTER_FEBIO_CLASS(FEPointBodyForce      , FEBodyForce, "point"      );

//-----------------------------------------------------------------------------
// classes derived from FEContactInterface
REGISTER_FEBIO_CLASS(FEFacet2FacetSliding   , FEContactInterface, "facet-to-facet sliding");
REGISTER_FEBIO_CLASS(FEPeriodicBoundary     , FEContactInterface, "periodic boundary"     );
REGISTER_FEBIO_CLASS(FERigidWallInterface   , FEContactInterface, "rigid_wall"            );
REGISTER_FEBIO_CLASS(FESlidingInterface     , FEContactInterface, "sliding_with_gaps"     );
REGISTER_FEBIO_CLASS(FESlidingInterfaceBW   , FEContactInterface, "sliding-Bonet-Wood"    );
REGISTER_FEBIO_CLASS(FESurfaceConstraint    , FEContactInterface, "surface constraint"    );
REGISTER_FEBIO_CLASS(FETiedInterface        , FEContactInterface, "tied"                  );

}
