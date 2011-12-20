#include "stdafx.h"
#include "FECore/febio.h"

#include "FEBioLib/FEConstBodyForce.h"
#include "FEBioLib/FEPointBodyForce.h"

#include "FEBioLib/FE2DFiberNeoHookean.h"
#include "FEBioLib/FE2DTransIsoMooneyRivlin.h"
#include "FEBioLib/FE2DTransIsoVerondaWestmann.h"
#include "FEBioLib/FEArrudaBoyce.h"
#include "FEBioLib/FEHolmesMow.h"
#include "FEBioLib/FEIncompNeoHookean.h"
#include "FEBioLib/FEIsotropicElastic.h"
#include "FEBioLib/FELinearElastic.h"
#include "FEBioLib/FELinearOrthotropic.h"
#include "FEBioLib/FEMooneyRivlin.h"
#include "FEBioLib/FENeoHookean.h"
#include "FEBioLib/FEOgdenMaterial.h"
#include "FEBioLib/FEVerondaWestmann.h"

//-----------------------------------------------------------------------------
// classes derived from FEBodyForce
REGISTER_FEBIO_CLASS(FEConstBodyForce      , FEBodyForce, "const"      );
REGISTER_FEBIO_CLASS(FENonConstBodyForce   , FEBodyForce, "non-const"  );
REGISTER_FEBIO_CLASS(FECentrifugalBodyForce, FEBodyForce, "centrifugal");
REGISTER_FEBIO_CLASS(FEPointBodyForce      , FEBodyForce, "point"      );

//-----------------------------------------------------------------------------
// material classes
REGISTER_MATERIAL(FE2DFiberNeoHookean        , "2D fiber neo-Hookean"         );
REGISTER_MATERIAL(FE2DTransIsoMooneyRivlin   , "2D trans iso Mooney-Rivlin"   );
REGISTER_MATERIAL(FE2DTransIsoVerondaWestmann, "2D trans iso Veronda-Westmann");
REGISTER_MATERIAL(FEArrudaBoyce              , "Arruda-Boyce"                 );
REGISTER_MATERIAL(FEHolmesMow                , "Holmes-Mow"                   );
REGISTER_MATERIAL(FEIncompNeoHookean         , "incomp neo-Hookean"           );
REGISTER_MATERIAL(FEIsotropicElastic         , "isotropic elastic"            );
REGISTER_MATERIAL(FELinearElastic            , "linear elastic"               );
REGISTER_MATERIAL(FELinearOrthotropic        , "linear orthotropic"           );
REGISTER_MATERIAL(FEMooneyRivlin             , "Mooney-Rivlin"                );
REGISTER_MATERIAL(FENeoHookean               , "neo-Hookean"                  );
REGISTER_MATERIAL(FEOgdenMaterial            , "Ogden"                        );
REGISTER_MATERIAL(FEVerondaWestmann          , "Veronda-Westmann"             );
