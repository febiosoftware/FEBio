#include "stdafx.h"
#include "FECore/febio.h"

#include "FEBioLib/FEConstBodyForce.h"
#include "FEBioLib/FEPointBodyForce.h"

#include "FEBioLib/FE2DFiberNeoHookean.h"
#include "FEBioLib/FE2DTransIsoMooneyRivlin.h"
#include "FEBioLib/FE2DTransIsoVerondaWestmann.h"
#include "FEBioLib/FEArrudaBoyce.h"
#include "FEBioLib/FECellGrowth.h"
#include "FEBioLib/FEDamageMooneyRivlin.h"
#include "FEBioLib/FEDamageNeoHookean.h"
#include "FEBioLib/FEDamageTransIsoMooneyRivlin.h"
#include "FEBioLib/FEDiscreteMaterial.h"
#include "FEBioLib/FEDonnanEquilibrium.h"
#include "FEBioLib/FEEFDDonnanEquilibrium.h"
#include "FEBioLib/FEEFDMooneyRivlin.h"
#include "FEBioLib/FEEFDNeoHookean.h"
#include "FEBioLib/FEEFDUncoupled.h"
#include "FEBioLib/FEEFDVerondaWestmann.h"
#include "FEBioLib/FEElasticMixture.h"
#include "FEBioLib/FEEllipsoidalFiberDistribution.h"
#include "FEBioLib/FEFiberExpPow.h"
#include "FEBioLib/FEFiberExpPowUncoupled.h"
#include "FEBioLib/FEFiberNeoHookean.h"
#include "FEBioLib/FEFungOrthoCompressible.h"
#include "FEBioLib/FEFungOrthotropic.h"
#include "FEBioLib/FEGasserOgdenHolzapfel.h"
#include "FEBioLib/FEHolmesMow.h"
#include "FEBioLib/FEIncompNeoHookean.h"
#include "FEBioLib/FEIsotropicElastic.h"
#include "FEBioLib/FEIsotropicFourier.h"
#include "FEBioLib/FELinearElastic.h"
#include "FEBioLib/FELinearOrthotropic.h"
#include "FEBioLib/FELinearTransIso.h"
#include "FEBioLib/FEMooneyRivlin.h"
#include "FEBioLib/FEMuscleMaterial.h"
#include "FEBioLib/FENeoHookean.h"
#include "FEBioLib/FEOgdenMaterial.h"
#include "FEBioLib/FEStVenantKirchhoff.h"
#include "FEBioLib/FETCNonlinearOrthotropic.h"
#include "FEBioLib/FETendonMaterial.h"
#include "FEBioLib/FETransIsoMooneyRivlin.h"
#include "FEBioLib/FETransIsoVerondaWestmann.h"
#include "FEBioLib/FETrussMaterial.h"
#include "FEBioLib/FEUncoupledElasticMixture.h"
#include "FEBioLib/FEVerondaWestmann.h"
#include "FEBioLib/FEViscoElasticMaterial.h"
#include "FEBioLib/FEVonMisesPlasticity.h"

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
REGISTER_MATERIAL(FECellGrowth                  , "cell growth"                   );
REGISTER_MATERIAL(FEDamageMooneyRivlin          , "damage Mooney-Rivlin"          );
REGISTER_MATERIAL(FEDamageNeoHookean            , "damage neo-Hookean"            );
REGISTER_MATERIAL(FEDamageTransIsoMooneyRivlin  , "damage trans iso Mooney-Rivlin");
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
REGISTER_MATERIAL(FEMuscleMaterial              , "muscle material"               );
REGISTER_MATERIAL(FENeoHookean                  , "neo-Hookean"                   );
REGISTER_MATERIAL(FEOgdenMaterial               , "Ogden"                         );
REGISTER_MATERIAL(FEStVenantKirchhoff           , "St.Venant-Kirchhoff"           );
REGISTER_MATERIAL(FETCNonlinearOrthotropic      , "TC nonlinear orthotropic"      );
REGISTER_MATERIAL(FETendonMaterial              , "tendon material"               );
REGISTER_MATERIAL(FETransIsoMooneyRivlin        , "trans iso Mooney-Rivlin"       );
REGISTER_MATERIAL(FETransIsoVerondaWestmann     , "trans iso Veronda-Westmann"    );
REGISTER_MATERIAL(FETrussMaterial               , "linear truss"                  );
REGISTER_MATERIAL(FEUncoupledElasticMixture     , "uncoupled solid mixture"       );
REGISTER_MATERIAL(FEVerondaWestmann             , "Veronda-Westmann"              );
REGISTER_MATERIAL(FEViscoElasticMaterial        , "viscoelastic"                  );
REGISTER_MATERIAL(FEVonMisesPlasticity          , "von-Mises plasticity"          );
