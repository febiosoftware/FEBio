#include "FEBioMix.h"
#include "FEBiphasic.h"
#include "FEBiphasicSolute.h"
#include "FEMultiphasic.h"
#include "FESolute.h"
#include "FETriphasic.h"
#include "FEDiffConstIso.h"
#include "FEDiffConstOrtho.h"
#include "FEDiffRefIso.h"
#include "FEDiffAlbroIso.h"
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
#include "FEReactionRateConst.h"
#include "FEReactionRateHuiskes.h"
#include "FEMassActionForward.h"
#include "FEMichaelisMenten.h"
#include "FEMassActionReversible.h"
#include "FEConcentrationIndependentReaction.h"

#include "FEPoroTraction.h"
#include "FEFluidFlux.h"
#include "FESoluteFlux.h"

#include "FESlidingInterface2.h"
#include "FESlidingInterface3.h"
#include "FESlidingInterfaceMP.h"
#include "FETiedBiphasicInterface.h"

#include "FEBiphasicAnalysis.h"
#include "FEBiphasicSoluteAnalysis.h"
#include "FEMultiphasicAnalysis.h"
#include "FEBiphasicSolver.h"
#include "FEBiphasicSoluteSolver.h"
#include "FEMultiphasicSolver.h"

#include "FEBioMixPlot.h"
#include "FEBioMixData.h"

#include "FEMixDomainFactory.h"

//-----------------------------------------------------------------------------
//! Initialization of the FEBioMix module. This function registers all the classes
//! in this module with the FEBio framework.
void FEBioMix::InitModule()
{
//-----------------------------------------------------------------------------
// Domain factory
	FEBioKernel& febio = FEBioKernel::GetInstance();
	febio.RegisterDomain(new FEMixDomainFactory);

//-----------------------------------------------------------------------------
// Analysis classes
REGISTER_FEBIO_CLASS(FEBiphasicAnalysis      , FEANALYSIS_ID, "biphasic"       );
REGISTER_FEBIO_CLASS(FEBiphasicSoluteAnalysis, FEANALYSIS_ID, "biphasic-solute");
REGISTER_FEBIO_CLASS(FEMultiphasicAnalysis   , FEANALYSIS_ID, "multiphasic"    );

//-----------------------------------------------------------------------------
// solver classes
REGISTER_FEBIO_CLASS(FEBiphasicSolver      , FESOLVER_ID, "biphasic"       );
REGISTER_FEBIO_CLASS(FEBiphasicSoluteSolver, FESOLVER_ID, "biphasic-solute");
REGISTER_FEBIO_CLASS(FEMultiphasicSolver   , FESOLVER_ID, "multiphasic"    );

//-----------------------------------------------------------------------------
// Materials
REGISTER_FEBIO_CLASS(FEBiphasic                          ,FEMATERIAL_ID, "biphasic"          );
REGISTER_FEBIO_CLASS(FEBiphasicSolute                    ,FEMATERIAL_ID, "biphasic-solute"   );
REGISTER_FEBIO_CLASS(FEMultiphasic                       ,FEMATERIAL_ID, "multiphasic"       );
REGISTER_FEBIO_CLASS(FESolute                            ,FEMATERIAL_ID, "solute"            );
REGISTER_FEBIO_CLASS(FETriphasic                         ,FEMATERIAL_ID, "triphasic"         );
REGISTER_FEBIO_CLASS(FEDiffConstIso                      ,FEMATERIAL_ID, "diff-const-iso"    );
REGISTER_FEBIO_CLASS(FEDiffConstOrtho                    ,FEMATERIAL_ID, "diff-const-ortho"  );
REGISTER_FEBIO_CLASS(FEDiffRefIso                        ,FEMATERIAL_ID, "diff-ref-iso"      );
REGISTER_FEBIO_CLASS(FEDiffAlbroIso                      ,FEMATERIAL_ID, "diff-Albro-iso"    );
REGISTER_FEBIO_CLASS(FEPermConstIso                      ,FEMATERIAL_ID, "perm-const-iso"    );
REGISTER_FEBIO_CLASS(FEPermHolmesMow                     ,FEMATERIAL_ID, "perm-Holmes-Mow"   );
REGISTER_FEBIO_CLASS(FEPermRefIso                        ,FEMATERIAL_ID, "perm-ref-iso"      );
REGISTER_FEBIO_CLASS(FEPermRefOrtho                      ,FEMATERIAL_ID, "perm-ref-ortho"    );
REGISTER_FEBIO_CLASS(FEPermRefTransIso                   ,FEMATERIAL_ID, "perm-ref-trans-iso");
REGISTER_FEBIO_CLASS(FEOsmCoefConst                      ,FEMATERIAL_ID, "osm-coef-const"    );
REGISTER_FEBIO_CLASS(FESFDSBM                            ,FEMATERIAL_ID, "spherical fiber distribution sbm");
REGISTER_FEBIO_CLASS(FESolventSupplyStarling             ,FEMATERIAL_ID, "Starling"          );
REGISTER_FEBIO_CLASS(FESolubConst                        ,FEMATERIAL_ID, "solub-const"       );
REGISTER_FEBIO_CLASS(FESupplyBinding                     ,FEMATERIAL_ID, "supply-binding"          );
REGISTER_FEBIO_CLASS(FESupplyConst                       ,FEMATERIAL_ID, "supply-const"            );
REGISTER_FEBIO_CLASS(FESupplySynthesisBinding            ,FEMATERIAL_ID, "supply-synthesis-binding");
REGISTER_FEBIO_CLASS(FESupplyMichaelisMenten             ,FEMATERIAL_ID, "supply-Michaelis-Menten" );
REGISTER_FEBIO_CLASS(FECarterHayes				    	 ,FEMATERIAL_ID, "Carter-Hayes"            );
REGISTER_FEBIO_CLASS(FEReactionRateConst		    	 ,FEMATERIAL_ID, "constant reaction rate"  );
REGISTER_FEBIO_CLASS(FEReactionRateHuiskes		    	 ,FEMATERIAL_ID, "Huiskes reaction rate"   );
REGISTER_FEBIO_CLASS(FEMassActionForward		    	 ,FEMATERIAL_ID, "mass-action-forward"     );
REGISTER_FEBIO_CLASS(FEConcentrationIndependentReaction  ,FEMATERIAL_ID, "concentration-independent");
REGISTER_FEBIO_CLASS(FEMassActionReversible              ,FEMATERIAL_ID, "mass-action-reversible");
REGISTER_FEBIO_CLASS(FEMichaelisMenten                   ,FEMATERIAL_ID, "Michaelis-Menten"      );
REGISTER_FEBIO_CLASS(FESolidBoundMolecule                ,FEMATERIAL_ID, "solid_bound"           );

//-----------------------------------------------------------------------------
// Surface loads
REGISTER_FEBIO_CLASS(FEPoroNormalTraction, FESURFACELOAD_ID, "normal_traction");
REGISTER_FEBIO_CLASS(FEFluidFlux         , FESURFACELOAD_ID, "fluidflux"      );
REGISTER_FEBIO_CLASS(FESoluteFlux        , FESURFACELOAD_ID, "soluteflux"     );

//-----------------------------------------------------------------------------
// Contact interfaces
REGISTER_FEBIO_CLASS(FESlidingInterface2    , FESURFACEPAIRINTERACTION_ID, "sliding2"              );
REGISTER_FEBIO_CLASS(FESlidingInterface3    , FESURFACEPAIRINTERACTION_ID, "sliding3"              );
REGISTER_FEBIO_CLASS(FESlidingInterfaceMP   , FESURFACEPAIRINTERACTION_ID, "sliding-multiphasic"   );
REGISTER_FEBIO_CLASS(FETiedBiphasicInterface, FESURFACEPAIRINTERACTION_ID, "tied-biphasic"         );

//-----------------------------------------------------------------------------
// classes derived from FEPlotData
REGISTER_FEBIO_CLASS(FEPlotEffectiveFluidPressure		, FEPLOTDATA_ID, "effective fluid pressure"        );
REGISTER_FEBIO_CLASS(FEPlotActualFluidPressure          , FEPLOTDATA_ID, "fluid pressure"                  );
REGISTER_FEBIO_CLASS(FEPlotFluidFlux                    , FEPLOTDATA_ID, "fluid flux"                      );
REGISTER_FEBIO_CLASS(FEPlotEffectiveSoluteConcentration , FEPLOTDATA_ID, "effective solute concentration"  );
REGISTER_FEBIO_CLASS(FEPlotActualSoluteConcentration    , FEPLOTDATA_ID, "solute concentration"            );
REGISTER_FEBIO_CLASS(FEPlotSoluteFlux                   , FEPLOTDATA_ID, "solute flux"                     );
REGISTER_FEBIO_CLASS(FEPlotReceptorLigandConcentration  , FEPLOTDATA_ID, "receptor-ligand concentration"   );
REGISTER_FEBIO_CLASS_T(FEPlotEffectiveSolConcentrationT , FEPLOTDATA_ID, 0, "effective solute 1 concentration");
REGISTER_FEBIO_CLASS_T(FEPlotActualSolConcentrationT    , FEPLOTDATA_ID, 0, "solute 1 concentration"          );
REGISTER_FEBIO_CLASS_T(FEPlotSolFluxT                   , FEPLOTDATA_ID, 0, "solute 1 flux"                   );
REGISTER_FEBIO_CLASS_T(FEPlotSBMConcentrationT          , FEPLOTDATA_ID, 0, "sbm 1 concentration"			   );
REGISTER_FEBIO_CLASS_T(FEPlotEffectiveSolConcentrationT , FEPLOTDATA_ID, 1, "effective solute 2 concentration");
REGISTER_FEBIO_CLASS_T(FEPlotActualSolConcentrationT    , FEPLOTDATA_ID, 1, "solute 2 concentration"          );
REGISTER_FEBIO_CLASS_T(FEPlotSolFluxT                   , FEPLOTDATA_ID, 1, "solute 2 flux"                   );
REGISTER_FEBIO_CLASS_T(FEPlotSBMConcentrationT          , FEPLOTDATA_ID, 1, "sbm 2 concentration"			   );
REGISTER_FEBIO_CLASS_T(FEPlotEffectiveSolConcentrationT , FEPLOTDATA_ID, 2, "effective solute 3 concentration");
REGISTER_FEBIO_CLASS_T(FEPlotActualSolConcentrationT    , FEPLOTDATA_ID, 2, "solute 3 concentration"          );
REGISTER_FEBIO_CLASS_T(FEPlotSolFluxT                   , FEPLOTDATA_ID, 2, "solute 3 flux"                   );
REGISTER_FEBIO_CLASS_T(FEPlotSBMConcentrationT          , FEPLOTDATA_ID, 2, "sbm 3 concentration"			   );
REGISTER_FEBIO_CLASS_T(FEPlotEffectiveSolConcentrationT , FEPLOTDATA_ID, 3, "effective solute 4 concentration");
REGISTER_FEBIO_CLASS_T(FEPlotActualSolConcentrationT    , FEPLOTDATA_ID, 3, "solute 4 concentration"          );
REGISTER_FEBIO_CLASS_T(FEPlotSolFluxT                   , FEPLOTDATA_ID, 3, "solute 4 flux"                   );
REGISTER_FEBIO_CLASS_T(FEPlotSBMConcentrationT          , FEPLOTDATA_ID, 3, "sbm 4 concentration"			   );
REGISTER_FEBIO_CLASS_T(FEPlotEffectiveSolConcentrationT , FEPLOTDATA_ID, 4, "effective solute 5 concentration");
REGISTER_FEBIO_CLASS_T(FEPlotActualSolConcentrationT    , FEPLOTDATA_ID, 4, "solute 5 concentration"          );
REGISTER_FEBIO_CLASS_T(FEPlotSolFluxT                   , FEPLOTDATA_ID, 4, "solute 5 flux"                   );
REGISTER_FEBIO_CLASS_T(FEPlotSBMConcentrationT          , FEPLOTDATA_ID, 4, "sbm 5 concentration"			   );
REGISTER_FEBIO_CLASS_T(FEPlotEffectiveSolConcentrationT , FEPLOTDATA_ID, 5, "effective solute 6 concentration");
REGISTER_FEBIO_CLASS_T(FEPlotActualSolConcentrationT    , FEPLOTDATA_ID, 5, "solute 6 concentration"          );
REGISTER_FEBIO_CLASS_T(FEPlotSolFluxT                   , FEPLOTDATA_ID, 5, "solute 6 flux"                   );
REGISTER_FEBIO_CLASS_T(FEPlotSBMConcentrationT          , FEPLOTDATA_ID, 5, "sbm 6 concentration"			   );
REGISTER_FEBIO_CLASS(FEPlotReferentialSolidVolumeFraction, FEPLOTDATA_ID, "referential solid volume fraction");
REGISTER_FEBIO_CLASS(FEPlotElectricPotential            , FEPLOTDATA_ID, "electric potential"  );
REGISTER_FEBIO_CLASS(FEPlotCurrentDensity               , FEPLOTDATA_ID, "current density"     );
REGISTER_FEBIO_CLASS(FEPlotFixedChargeDensity           , FEPLOTDATA_ID, "fixed charge density");
REGISTER_FEBIO_CLASS(FEPlotReferentialFixedChargeDensity, FEPLOTDATA_ID, "referential fixed charge density");
REGISTER_FEBIO_CLASS(FEPlotNodalFluidFlux               , FEPLOTDATA_ID, "nodal fluid flux"    );
REGISTER_FEBIO_CLASS(FEPlotOsmolarity                   , FEPLOTDATA_ID,  "osmolarity"         );
REGISTER_FEBIO_CLASS(FEPlotFluidForce                   , FEPLOTDATA_ID, "fluid force"         );
REGISTER_FEBIO_CLASS_T(FEPlotSBMRefAppDensityT          , FEPLOTDATA_ID, 0, "sbm 1 referential apparent density");
REGISTER_FEBIO_CLASS_T(FEPlotSBMRefAppDensityT          , FEPLOTDATA_ID, 1, "sbm 2 referential apparent density");
REGISTER_FEBIO_CLASS_T(FEPlotSBMRefAppDensityT          , FEPLOTDATA_ID, 2, "sbm 3 referential apparent density");
REGISTER_FEBIO_CLASS_T(FEPlotSBMRefAppDensityT          , FEPLOTDATA_ID, 3, "sbm 4 referential apparent density");
REGISTER_FEBIO_CLASS_T(FEPlotSBMRefAppDensityT          , FEPLOTDATA_ID, 4, "sbm 5 referential apparent density");
REGISTER_FEBIO_CLASS_T(FEPlotSBMRefAppDensityT          , FEPLOTDATA_ID, 5, "sbm 6 referential apparent density");
REGISTER_FEBIO_CLASS_T(FEPlotSBMRefAppDensityT          , FEPLOTDATA_ID, 6, "sbm 7 referential apparent density");
REGISTER_FEBIO_CLASS_T(FEPlotSBMRefAppDensityT          , FEPLOTDATA_ID, 7, "sbm 8 referential apparent density");

//-----------------------------------------------------------------------------
REGISTER_FEBIO_CLASS(FENodeTemp, FENODELOGDATA_ID, "T");
REGISTER_FEBIO_CLASS(FENodePressure, FENODELOGDATA_ID, "p");
REGISTER_FEBIO_CLASS(FENodeConcentration, FENODELOGDATA_ID, "c");
REGISTER_FEBIO_CLASS_T(FENodeConcentration_T, FENODELOGDATA_ID, 0, "c1");
REGISTER_FEBIO_CLASS_T(FENodeConcentration_T, FENODELOGDATA_ID, 1, "c2");
REGISTER_FEBIO_CLASS_T(FENodeConcentration_T, FENODELOGDATA_ID, 2, "c3");
REGISTER_FEBIO_CLASS_T(FENodeConcentration_T, FENODELOGDATA_ID, 3, "c4");
REGISTER_FEBIO_CLASS_T(FENodeConcentration_T, FENODELOGDATA_ID, 4, "c5");
REGISTER_FEBIO_CLASS_T(FENodeConcentration_T, FENODELOGDATA_ID, 5, "c6");
REGISTER_FEBIO_CLASS_T(FENodeConcentration_T, FENODELOGDATA_ID, 6, "c7");
REGISTER_FEBIO_CLASS_T(FENodeConcentration_T, FENODELOGDATA_ID, 7, "c8");

//-----------------------------------------------------------------------------
REGISTER_FEBIO_CLASS(FELogElemFluidPressure, FEELEMLOGDATA_ID, "p");
REGISTER_FEBIO_CLASS(FELogElemFluidFluxX, FEELEMLOGDATA_ID, "wx");
REGISTER_FEBIO_CLASS(FELogElemFluidFluxY, FEELEMLOGDATA_ID, "wy");
REGISTER_FEBIO_CLASS(FELogElemFluidFluxZ, FEELEMLOGDATA_ID, "wz");
REGISTER_FEBIO_CLASS(FELogElemSoluteConcentration, FEELEMLOGDATA_ID, "c");
REGISTER_FEBIO_CLASS(FELogElemSoluteFluxX, FEELEMLOGDATA_ID, "jx");
REGISTER_FEBIO_CLASS(FELogElemSoluteFluxY, FEELEMLOGDATA_ID, "jy");
REGISTER_FEBIO_CLASS(FELogElemSoluteFluxZ, FEELEMLOGDATA_ID, "jz");
REGISTER_FEBIO_CLASS(FELogElemSoluteRefConcentration, FEELEMLOGDATA_ID, "crc");
REGISTER_FEBIO_CLASS(FELogElemElectricPotential, FEELEMLOGDATA_ID, "psi");
REGISTER_FEBIO_CLASS(FELogElemCurrentDensityX, FEELEMLOGDATA_ID, "Iex");
REGISTER_FEBIO_CLASS(FELogElemCurrentDensityY, FEELEMLOGDATA_ID, "Iey");
REGISTER_FEBIO_CLASS(FELogElemCurrentDensityZ, FEELEMLOGDATA_ID, "Iez");

REGISTER_FEBIO_CLASS_T(FELogElemSoluteConcentration_T, FEELEMLOGDATA_ID, 0, "c1");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteConcentration_T, FEELEMLOGDATA_ID, 1, "c2");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteConcentration_T, FEELEMLOGDATA_ID, 2, "c3");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteConcentration_T, FEELEMLOGDATA_ID, 3, "c4");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteConcentration_T, FEELEMLOGDATA_ID, 4, "c5");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteConcentration_T, FEELEMLOGDATA_ID, 5, "c6");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteConcentration_T, FEELEMLOGDATA_ID, 6, "c7");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteConcentration_T, FEELEMLOGDATA_ID, 7, "c8");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxX_T, FEELEMLOGDATA_ID, 0, "j1x");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxY_T, FEELEMLOGDATA_ID, 0, "j1y");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxZ_T, FEELEMLOGDATA_ID, 0, "j1z");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxX_T, FEELEMLOGDATA_ID, 1, "j2x");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxY_T, FEELEMLOGDATA_ID, 1, "j2y");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxZ_T, FEELEMLOGDATA_ID, 1, "j2z");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxX_T, FEELEMLOGDATA_ID, 2, "j3x");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxY_T, FEELEMLOGDATA_ID, 2, "j3y");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxZ_T, FEELEMLOGDATA_ID, 2, "j3z");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxX_T, FEELEMLOGDATA_ID, 3, "j4x");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxY_T, FEELEMLOGDATA_ID, 3, "j4y");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxZ_T, FEELEMLOGDATA_ID, 3, "j4z");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxX_T, FEELEMLOGDATA_ID, 4, "j5x");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxY_T, FEELEMLOGDATA_ID, 4, "j5y");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxZ_T, FEELEMLOGDATA_ID, 4, "j5z");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxX_T, FEELEMLOGDATA_ID, 5, "j6x");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxY_T, FEELEMLOGDATA_ID, 5, "j6y");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxZ_T, FEELEMLOGDATA_ID, 5, "j6z");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxX_T, FEELEMLOGDATA_ID, 6, "j7x");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxY_T, FEELEMLOGDATA_ID, 6, "j7y");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxZ_T, FEELEMLOGDATA_ID, 6, "j7z");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxX_T, FEELEMLOGDATA_ID, 7, "j8x");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxY_T, FEELEMLOGDATA_ID, 7, "j8y");
REGISTER_FEBIO_CLASS_T(FELogElemSoluteFluxZ_T, FEELEMLOGDATA_ID, 7, "j8z");
REGISTER_FEBIO_CLASS_T(FELogElemSBMConcentration_T, FEELEMLOGDATA_ID, 0, "sbm1");
REGISTER_FEBIO_CLASS_T(FELogElemSBMConcentration_T, FEELEMLOGDATA_ID, 1, "sbm2");
REGISTER_FEBIO_CLASS_T(FELogElemSBMConcentration_T, FEELEMLOGDATA_ID, 2, "sbm3");
REGISTER_FEBIO_CLASS_T(FELogElemSBMConcentration_T, FEELEMLOGDATA_ID, 3, "sbm4");
REGISTER_FEBIO_CLASS_T(FELogElemSBMConcentration_T, FEELEMLOGDATA_ID, 4, "sbm5");
REGISTER_FEBIO_CLASS_T(FELogElemSBMConcentration_T, FEELEMLOGDATA_ID, 5, "sbm6");
REGISTER_FEBIO_CLASS_T(FELogElemSBMConcentration_T, FEELEMLOGDATA_ID, 6, "sbm7");
REGISTER_FEBIO_CLASS_T(FELogElemSBMConcentration_T, FEELEMLOGDATA_ID, 7, "sbm8");
}
