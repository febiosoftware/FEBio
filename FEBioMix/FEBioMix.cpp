#include "FEBioMix.h"
#include "FEBiphasic.h"
#include "FEBiphasicSolute.h"
#include "FEMultiphasic.h"
#include "FESolute.h"
#include "FETriphasic.h"
#include "FEDiffConstIso.h"
#include "FEDiffConstOrtho.h"
#include "FEDiffRefIso.h"
#include "FEPermConstIso.h"
#include "FEPermHolmesMow.h"
#include "FEPermRefIso.h"
#include "FEPermRefOrtho.h"
#include "FEPermRefTransIso.h"
#include "FEOsmCoefConst.h"
#include "FESFDSBM.h"
#include "FESolventSupplyStarling.h"
#include "FESolubConst.h"
#include "FESupplyBinding.h"
#include "FESupplyConst.h"
#include "FESupplySynthesisBinding.h"
#include "FESupplyMichaelisMenten.h"
#include "FECarterHayes.h"

#include "FESlidingInterface2.h"
#include "FESlidingInterface3.h"
#include "FETiedBiphasicInterface.h"

#include "FEBioMixPlot.h"
#include "FEBioMixData.h"

//-----------------------------------------------------------------------------
//! Initialization of the FEBioMix module. This function registers all the classes
//! in this module with the FEBio framework.
void FEBioMix::InitModule()
{

//-----------------------------------------------------------------------------
// Materials
REGISTER_MATERIAL(FEBiphasic                     , "biphasic"          );
REGISTER_MATERIAL(FEBiphasicSolute               , "biphasic-solute"   );
REGISTER_MATERIAL(FEMultiphasic                  , "multiphasic"       );
REGISTER_MATERIAL(FESolute                       , "solute"            );
REGISTER_MATERIAL(FETriphasic                    , "triphasic"         );

REGISTER_MATERIAL(FEDiffConstIso                 , "diff-const-iso"    );
REGISTER_MATERIAL(FEDiffConstOrtho               , "diff-const-ortho"  );
REGISTER_MATERIAL(FEDiffRefIso                   , "diff-ref-iso"      );
REGISTER_MATERIAL(FEPermConstIso                 , "perm-const-iso"    );
REGISTER_MATERIAL(FEPermHolmesMow                , "perm-Holmes-Mow"   );
REGISTER_MATERIAL(FEPermRefIso                   , "perm-ref-iso"      );
REGISTER_MATERIAL(FEPermRefOrtho                 , "perm-ref-ortho"    );
REGISTER_MATERIAL(FEPermRefTransIso              , "perm-ref-trans-iso");
REGISTER_MATERIAL(FEOsmCoefConst                 , "osm-coef-const"    );
REGISTER_MATERIAL(FESFDSBM                       , "spherical fiber distribution sbm");
REGISTER_MATERIAL(FESolventSupplyStarling        , "Starling"          );
REGISTER_MATERIAL(FESolubConst                   , "solub-const"       );
REGISTER_MATERIAL(FESupplyBinding                , "supply-binding"          );
REGISTER_MATERIAL(FESupplyConst                  , "supply-const"            );
REGISTER_MATERIAL(FESupplySynthesisBinding       , "supply-synthesis-binding");
REGISTER_MATERIAL(FESupplyMichaelisMenten        , "supply-Michaelis-Menten" );
REGISTER_MATERIAL(FECarterHayes					 , "Carter-Hayes"             );

//-----------------------------------------------------------------------------
// Contact interfaces
REGISTER_FEBIO_CLASS(FESlidingInterface2    , FEContactInterface, "sliding2"              );
REGISTER_FEBIO_CLASS(FESlidingInterface3    , FEContactInterface, "sliding3"              );
REGISTER_FEBIO_CLASS(FETiedBiphasicInterface, FEContactInterface, "tied-biphasic"         );

//-----------------------------------------------------------------------------
// classes derived from FEPlotData
REGISTER_FEBIO_CLASS(FEPlotEffectiveFluidPressure		, FEPlotData, "effective fluid pressure"        );
REGISTER_FEBIO_CLASS(FEPlotActualFluidPressure          , FEPlotData, "fluid pressure"                  );
REGISTER_FEBIO_CLASS(FEPlotFluidFlux                    , FEPlotData, "fluid flux"                      );
REGISTER_FEBIO_CLASS(FEPlotEffectiveSoluteConcentration , FEPlotData, "effective solute concentration"  );
REGISTER_FEBIO_CLASS(FEPlotActualSoluteConcentration    , FEPlotData, "solute concentration"            );
REGISTER_FEBIO_CLASS(FEPlotSoluteFlux                   , FEPlotData, "solute flux"                     );
REGISTER_FEBIO_CLASS(FEPlotReceptorLigandConcentration  , FEPlotData, "receptor-ligand concentration"   );
REGISTER_FEBIO_CLASS_T(FEPlotEffectiveSolConcentrationT , FEPlotData, 0, "effective solute 1 concentration");
REGISTER_FEBIO_CLASS_T(FEPlotActualSolConcentrationT    , FEPlotData, 0, "solute 1 concentration"          );
REGISTER_FEBIO_CLASS_T(FEPlotSolFluxT                   , FEPlotData, 0, "solute 1 flux"                   );
REGISTER_FEBIO_CLASS_T(FEPlotSBMConcentrationT          , FEPlotData, 0, "sbm 1 concentration"			   );
REGISTER_FEBIO_CLASS_T(FEPlotEffectiveSolConcentrationT , FEPlotData, 1, "effective solute 2 concentration");
REGISTER_FEBIO_CLASS_T(FEPlotActualSolConcentrationT    , FEPlotData, 1, "solute 2 concentration"          );
REGISTER_FEBIO_CLASS_T(FEPlotSolFluxT                   , FEPlotData, 1, "solute 2 flux"                   );
REGISTER_FEBIO_CLASS_T(FEPlotSBMConcentrationT          , FEPlotData, 1, "sbm 2 concentration"			   );
REGISTER_FEBIO_CLASS_T(FEPlotEffectiveSolConcentrationT , FEPlotData, 2, "effective solute 3 concentration");
REGISTER_FEBIO_CLASS_T(FEPlotActualSolConcentrationT    , FEPlotData, 2, "solute 3 concentration"          );
REGISTER_FEBIO_CLASS_T(FEPlotSolFluxT                   , FEPlotData, 2, "solute 3 flux"                   );
REGISTER_FEBIO_CLASS_T(FEPlotSBMConcentrationT          , FEPlotData, 2, "sbm 3 concentration"			   );
REGISTER_FEBIO_CLASS_T(FEPlotEffectiveSolConcentrationT , FEPlotData, 3, "effective solute 4 concentration");
REGISTER_FEBIO_CLASS_T(FEPlotActualSolConcentrationT    , FEPlotData, 3, "solute 4 concentration"          );
REGISTER_FEBIO_CLASS_T(FEPlotSolFluxT                   , FEPlotData, 3, "solute 4 flux"                   );
REGISTER_FEBIO_CLASS_T(FEPlotSBMConcentrationT          , FEPlotData, 3, "sbm 4 concentration"			   );
REGISTER_FEBIO_CLASS_T(FEPlotEffectiveSolConcentrationT , FEPlotData, 4, "effective solute 5 concentration");
REGISTER_FEBIO_CLASS_T(FEPlotActualSolConcentrationT    , FEPlotData, 4, "solute 5 concentration"          );
REGISTER_FEBIO_CLASS_T(FEPlotSolFluxT                   , FEPlotData, 4, "solute 5 flux"                   );
REGISTER_FEBIO_CLASS_T(FEPlotSBMConcentrationT          , FEPlotData, 4, "sbm 5 concentration"			   );
REGISTER_FEBIO_CLASS_T(FEPlotEffectiveSolConcentrationT , FEPlotData, 5, "effective solute 6 concentration");
REGISTER_FEBIO_CLASS_T(FEPlotActualSolConcentrationT    , FEPlotData, 5, "solute 6 concentration"          );
REGISTER_FEBIO_CLASS_T(FEPlotSolFluxT                   , FEPlotData, 5, "solute 6 flux"                   );
REGISTER_FEBIO_CLASS_T(FEPlotSBMConcentrationT          , FEPlotData, 5, "sbm 6 concentration"			   );
REGISTER_FEBIO_CLASS(FEPlotReferentialSolidVolumeFraction, FEPlotData, "referential solid volume fraction");
REGISTER_FEBIO_CLASS(FEPlotElectricPotential            , FEPlotData, "electric potential"  );
REGISTER_FEBIO_CLASS(FEPlotCurrentDensity               , FEPlotData, "current density"     );
REGISTER_FEBIO_CLASS(FEPlotFixedChargeDensity           , FEPlotData, "fixed charge density");
REGISTER_FEBIO_CLASS(FEPlotReferentialFixedChargeDensity, FEPlotData, "referential fixed charge density");
REGISTER_FEBIO_CLASS(FEPlotNodalFluidFlux               , FEPlotData, "nodal fluid flux"    );
REGISTER_FEBIO_CLASS(FEPlotOsmolarity                   , FEPlotData,  "osmolarity"         );
REGISTER_FEBIO_CLASS_T(FEPlotSBMRefAppDensityT          , FEPlotData, 0, "sbm 1 referential apparent density");
REGISTER_FEBIO_CLASS_T(FEPlotSBMRefAppDensityT          , FEPlotData, 1, "sbm 2 referential apparent density");
REGISTER_FEBIO_CLASS_T(FEPlotSBMRefAppDensityT          , FEPlotData, 2, "sbm 3 referential apparent density");
REGISTER_FEBIO_CLASS_T(FEPlotSBMRefAppDensityT          , FEPlotData, 3, "sbm 4 referential apparent density");
REGISTER_FEBIO_CLASS_T(FEPlotSBMRefAppDensityT          , FEPlotData, 4, "sbm 5 referential apparent density");
REGISTER_FEBIO_CLASS_T(FEPlotSBMRefAppDensityT          , FEPlotData, 5, "sbm 6 referential apparent density");
REGISTER_FEBIO_CLASS_T(FEPlotSBMRefAppDensityT          , FEPlotData, 6, "sbm 7 referential apparent density");
REGISTER_FEBIO_CLASS_T(FEPlotSBMRefAppDensityT          , FEPlotData, 7, "sbm 8 referential apparent density");

//-----------------------------------------------------------------------------
REGISTER_FEBIO_CLASS(FENodeTemp, FENodeLogData, "T");
REGISTER_FEBIO_CLASS(FENodePressure, FENodeLogData, "p");
REGISTER_FEBIO_CLASS(FENodeConcentration, FENodeLogData, "c");
REGISTER_FEBIO_CLASS_T(FENodeConcentration_T, FENodeLogData, 0, "c1");
REGISTER_FEBIO_CLASS_T(FENodeConcentration_T, FENodeLogData, 1, "c2");
REGISTER_FEBIO_CLASS_T(FENodeConcentration_T, FENodeLogData, 2, "c3");
REGISTER_FEBIO_CLASS_T(FENodeConcentration_T, FENodeLogData, 3, "c4");
REGISTER_FEBIO_CLASS_T(FENodeConcentration_T, FENodeLogData, 4, "c5");
REGISTER_FEBIO_CLASS_T(FENodeConcentration_T, FENodeLogData, 5, "c6");
REGISTER_FEBIO_CLASS_T(FENodeConcentration_T, FENodeLogData, 6, "c7");
REGISTER_FEBIO_CLASS_T(FENodeConcentration_T, FENodeLogData, 7, "c8");

//-----------------------------------------------------------------------------
REGISTER_FEBIO_CLASS(FELogElemFluidPressure, FELogElemData, "p");
REGISTER_FEBIO_CLASS(FELogElemFluidFluxX, FELogElemData, "wx");
REGISTER_FEBIO_CLASS(FELogElemFluidFluxY, FELogElemData, "wy");
REGISTER_FEBIO_CLASS(FELogElemFluidFluxZ, FELogElemData, "wz");
REGISTER_FEBIO_CLASS(FELogElemSoluteConcentration, FELogElemData, "c");
REGISTER_FEBIO_CLASS(FELogElemSoluteFluxX, FELogElemData, "jx");
REGISTER_FEBIO_CLASS(FELogElemSoluteFluxY, FELogElemData, "jy");
REGISTER_FEBIO_CLASS(FELogElemSoluteFluxZ, FELogElemData, "jz");
REGISTER_FEBIO_CLASS(FELogElemSoluteRefConcentration, FELogElemData, "crc");
REGISTER_FEBIO_CLASS(FELogElemElectricPotential, FELogElemData, "psi");
REGISTER_FEBIO_CLASS(FELogElemCurrentDensityX, FELogElemData, "Iex");
REGISTER_FEBIO_CLASS(FELogElemCurrentDensityY, FELogElemData, "Iey");
REGISTER_FEBIO_CLASS(FELogElemCurrentDensityZ, FELogElemData, "Iez");

REGISTER_FEBIO_CLASS_T(FELogElemSoluteConcentration_T, FELogElemData, 0, "c1");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteConcentration_T, FELogElemData, 1, "c2");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteConcentration_T, FELogElemData, 2, "c3");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteConcentration_T, FELogElemData, 3, "c4");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteConcentration_T, FELogElemData, 4, "c5");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteConcentration_T, FELogElemData, 5, "c6");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteConcentration_T, FELogElemData, 6, "c7");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteConcentration_T, FELogElemData, 7, "c8");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxX_T, FELogElemData, 0, "j1x");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxY_T, FELogElemData, 0, "j1y");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxZ_T, FELogElemData, 0, "j1z");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxX_T, FELogElemData, 1, "j2x");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxY_T, FELogElemData, 1, "j2y");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxZ_T, FELogElemData, 1, "j2z");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxX_T, FELogElemData, 2, "j3x");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxY_T, FELogElemData, 2, "j3y");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxZ_T, FELogElemData, 2, "j3z");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxX_T, FELogElemData, 3, "j4x");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxY_T, FELogElemData, 3, "j4y");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxZ_T, FELogElemData, 3, "j4z");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxX_T, FELogElemData, 4, "j5x");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxY_T, FELogElemData, 4, "j5y");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxZ_T, FELogElemData, 4, "j5z");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxX_T, FELogElemData, 5, "j6x");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxY_T, FELogElemData, 5, "j6y");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxZ_T, FELogElemData, 5, "j6z");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxX_T, FELogElemData, 6, "j7x");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxY_T, FELogElemData, 6, "j7y");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxZ_T, FELogElemData, 6, "j7z");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxX_T, FELogElemData, 7, "j8x");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxY_T, FELogElemData, 7, "j8y");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxZ_T, FELogElemData, 7, "j8z");
REGISTER_FEBIO_CLASS_T(FELogElemSBMConcentration_T, FELogElemData, 0, "sbm1");
REGISTER_FEBIO_CLASS_T(FELogElemSBMConcentration_T, FELogElemData, 1, "sbm2");
REGISTER_FEBIO_CLASS_T(FELogElemSBMConcentration_T, FELogElemData, 2, "sbm3");
REGISTER_FEBIO_CLASS_T(FELogElemSBMConcentration_T, FELogElemData, 3, "sbm4");
REGISTER_FEBIO_CLASS_T(FELogElemSBMConcentration_T, FELogElemData, 4, "sbm5");
REGISTER_FEBIO_CLASS_T(FELogElemSBMConcentration_T, FELogElemData, 5, "sbm6");
REGISTER_FEBIO_CLASS_T(FELogElemSBMConcentration_T, FELogElemData, 6, "sbm7");
REGISTER_FEBIO_CLASS_T(FELogElemSBMConcentration_T, FELogElemData, 7, "sbm8");
}
