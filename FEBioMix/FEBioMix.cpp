#include "FEBioMix.h"
#include "FEDiffConstIso.h"
#include "FEDiffConstOrtho.h"
#include "FEDiffRefIso.h"
#include "FEPermConstIso.h"
#include "FEPermHolmesMow.h"
#include "FEPermRefIso.h"
#include "FEPermRefOrtho.h"
#include "FEPermRefTransIso.h"
#include "FEOsmCoefConst.h"
#include "FESolventSupplyStarling.h"
#include "FESolubConst.h"
#include "FESupplyBinding.h"
#include "FESupplyConst.h"
#include "FESupplySynthesisBinding.h"
#include "FESupplyMichaelisMenten.h"

//-----------------------------------------------------------------------------
//! Initialization of the FEBioMix module. This function registers all the classes
//! in this module with the FEBio framework.
void FEBioMix::InitModule()
{
	// Materials
	REGISTER_MATERIAL(FEDiffConstIso                 , "diff-const-iso"    );
	REGISTER_MATERIAL(FEDiffConstOrtho               , "diff-const-ortho"  );
	REGISTER_MATERIAL(FEDiffRefIso                   , "diff-ref-iso"      );
	REGISTER_MATERIAL(FEPermConstIso                 , "perm-const-iso"    );
	REGISTER_MATERIAL(FEPermHolmesMow                , "perm-Holmes-Mow"   );
	REGISTER_MATERIAL(FEPermRefIso                   , "perm-ref-iso"      );
	REGISTER_MATERIAL(FEPermRefOrtho                 , "perm-ref-ortho"    );
	REGISTER_MATERIAL(FEPermRefTransIso              , "perm-ref-trans-iso");
	REGISTER_MATERIAL(FEOsmCoefConst                 , "osm-coef-const"    );
	REGISTER_MATERIAL(FESolventSupplyStarling        , "Starling"          );
	REGISTER_MATERIAL(FESolubConst                   , "solub-const"       );
	REGISTER_MATERIAL(FESupplyBinding                , "supply-binding"          );
	REGISTER_MATERIAL(FESupplyConst                  , "supply-const"            );
	REGISTER_MATERIAL(FESupplySynthesisBinding       , "supply-synthesis-binding");
	REGISTER_MATERIAL(FESupplyMichaelisMenten        , "supply-Michaelis-Menten" );
}
