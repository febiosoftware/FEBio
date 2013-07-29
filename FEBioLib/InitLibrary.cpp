#include "stdafx.h"
#include "FECore/febio.h"

#include "FEBioPlot/FEPlotDomainData.h"
#include "FEBioPlot/FEPlotNodeData.h"
#include "FEBioPlot/FEPlotSurfaceData.h"
#include "FECore/FECoordSysMap.h"

#include "FEBioMech/FEBioMech.h"
#include "FEBioMix/FEBioMix.h"
#include "FEBioHeat/FEBioHeat.h"

void InitFEBioLibrary()
{
//-----------------------------------------------------------------------------
// import all modules
FEBioMech::InitModule();
FEBioMix::InitModule();
FEBioHeat::InitModule();

//-----------------------------------------------------------------------------
// classes derived from FEPlotData
REGISTER_FEBIO_CLASS(FEPlotEffectiveFluidPressure		, FEPlotData, "effective fluid pressure"        );
REGISTER_FEBIO_CLASS(FEPlotActualFluidPressure          , FEPlotData, "fluid pressure"                  );
REGISTER_FEBIO_CLASS(FEPlotStrainEnergyDensity          , FEPlotData, "strain energy density"           );
REGISTER_FEBIO_CLASS(FEPlotSpecificStrainEnergy         , FEPlotData, "specific strain energy"          );
REGISTER_FEBIO_CLASS(FEPlotDensity                      , FEPlotData, "density"                         );
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
REGISTER_FEBIO_CLASS(FEPlotUT4NodalStresses             , FEPlotData, "ut4 nodal stress"    );

//-----------------------------------------------------------------------------
REGISTER_FEBIO_CLASS(FEPlotNodeDisplacement  , FEPlotData, "displacement"   );
REGISTER_FEBIO_CLASS(FEPlotNodeVelocity      , FEPlotData, "velocity"       );
REGISTER_FEBIO_CLASS(FEPlotNodeAcceleration  , FEPlotData, "acceleration"   );
REGISTER_FEBIO_CLASS(FEPlotNodeReactionForces, FEPlotData, "reaction forces");

//-----------------------------------------------------------------------------
REGISTER_FEBIO_CLASS(FEPlotContactGap      , FEPlotData, "contact gap"     );
REGISTER_FEBIO_CLASS(FEPlotContactPressure , FEPlotData, "contact pressure");
REGISTER_FEBIO_CLASS(FEPlotContactTraction , FEPlotData, "contact traction");

//-----------------------------------------------------------------------------
// Classes derived from FECoordSysMap
REGISTER_FEBIO_CLASS(FELocalMap      , FECoordSysMap, "local"      );
REGISTER_FEBIO_CLASS(FESphericalMap  , FECoordSysMap, "spherical"  );
REGISTER_FEBIO_CLASS(FECylindricalMap, FECoordSysMap, "cylindrical");
REGISTER_FEBIO_CLASS(FEVectorMap     , FECoordSysMap, "vector"     );
}
