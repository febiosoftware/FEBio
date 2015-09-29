#include "FEBioFluid.h"
#include "FEFluid.h"
#include "FEIdealGas.h"
#include "FEIdealFluid.h"
#include "FENeoHookeanFluid.h"
#include "FENewtonianFluid.h"
#include "FECarreauFluid.h"

#include "FEFluidAnalysis.h"
#include "FEFluidSolver.h"
#include "FEFluidDomain.h"

#include "FEFluidTractionLoad.h"

#include "FEBioFluidPlot.h"
#include "FEBioFluidData.h"

#include "FEFluidDomainFactory.h"

//-----------------------------------------------------------------------------
//! Initialization of the FEBioFluid module. This function registers all the classes
//! in this module with the FEBio framework.
void FEBioFluid::InitModule()
{
//-----------------------------------------------------------------------------
// Domain factory
	FECoreKernel& febio = FECoreKernel::GetInstance();
	febio.RegisterDomain(new FEFluidDomainFactory);

//-----------------------------------------------------------------------------
// Analysis classes
REGISTER_FECORE_CLASS(FEFluidAnalysis      , FEANALYSIS_ID, "fluid"       );

//-----------------------------------------------------------------------------
// solver classes
REGISTER_FECORE_CLASS(FEFluidSolver        , FESOLVER_ID, "fluid"         );

//-----------------------------------------------------------------------------
// Materials
REGISTER_FECORE_CLASS(FEFluid                       ,FEMATERIAL_ID, "fluid"             );
REGISTER_FECORE_CLASS(FEIdealGas                    ,FEMATERIAL_ID, "ideal gas"         );
REGISTER_FECORE_CLASS(FEIdealFluid                  ,FEMATERIAL_ID, "ideal fluid"       );
REGISTER_FECORE_CLASS(FENeoHookeanFluid             ,FEMATERIAL_ID, "neo-Hookean fluid" );
REGISTER_FECORE_CLASS(FENewtonianFluid              ,FEMATERIAL_ID, "Newtonian fluid"   );
REGISTER_FECORE_CLASS(FECarreauFluid                ,FEMATERIAL_ID, "Carreau fluid"     );

//-----------------------------------------------------------------------------
// Domain classes
REGISTER_FECORE_CLASS(FEFluidDomain                 , FEDOMAIN_ID, "fluid"              );

//-----------------------------------------------------------------------------
// Surface loads
REGISTER_FECORE_CLASS(FEFluidTractionLoad           , FESURFACELOAD_ID, "fluid traction");
    

//-----------------------------------------------------------------------------
// classes derived from FEPlotData
REGISTER_FECORE_CLASS(FEPlotFluidDilatation         , FEPLOTDATA_ID, "fluid dilatation"   );
REGISTER_FECORE_CLASS(FEPlotElasticFluidPressure	, FEPLOTDATA_ID, "elastic fluid pressure"   );
REGISTER_FECORE_CLASS(FEPlotFluidVolumeRatio		, FEPLOTDATA_ID, "fluid volume ratio"       );
REGISTER_FECORE_CLASS(FEPlotFluidDensity            , FEPLOTDATA_ID, "fluid density"            );
REGISTER_FECORE_CLASS(FEPlotFluidVelocity           , FEPLOTDATA_ID, "fluid velocity"           );
REGISTER_FECORE_CLASS(FEPlotFluidAcceleration       , FEPLOTDATA_ID, "fluid acceleration"       );
REGISTER_FECORE_CLASS(FEPlotFluidVorticity          , FEPLOTDATA_ID, "fluid vorticity"          );
REGISTER_FECORE_CLASS(FEPlotElementFluidStress      , FEPLOTDATA_ID, "fluid stress"             );
REGISTER_FECORE_CLASS(FEPlotElementFluidRateOfDef   , FEPLOTDATA_ID, "fluid rate of deformation");

//-----------------------------------------------------------------------------
REGISTER_FECORE_CLASS(FELogElasticFluidPressure, FEELEMLOGDATA_ID, "fp");
REGISTER_FECORE_CLASS(FELogFluidVolumeRatio,     FEELEMLOGDATA_ID, "fJ");
REGISTER_FECORE_CLASS(FELogFluidAccelerationX,   FEELEMLOGDATA_ID, "fax");
REGISTER_FECORE_CLASS(FELogFluidAccelerationY,   FEELEMLOGDATA_ID, "fay");
REGISTER_FECORE_CLASS(FELogFluidAccelerationZ,   FEELEMLOGDATA_ID, "faz");
REGISTER_FECORE_CLASS(FELogFluidVorticityX,      FEELEMLOGDATA_ID, "fwx");
REGISTER_FECORE_CLASS(FELogFluidVorticityY,      FEELEMLOGDATA_ID, "fwy");
REGISTER_FECORE_CLASS(FELogFluidStressXX,        FEELEMLOGDATA_ID, "fsxx");
REGISTER_FECORE_CLASS(FELogFluidStressYY,        FEELEMLOGDATA_ID, "fsyy");
REGISTER_FECORE_CLASS(FELogFluidStressZZ,        FEELEMLOGDATA_ID, "fszz");
REGISTER_FECORE_CLASS(FELogFluidStressXY,        FEELEMLOGDATA_ID, "fsxy");
REGISTER_FECORE_CLASS(FELogFluidStressYZ,        FEELEMLOGDATA_ID, "fsyz");
REGISTER_FECORE_CLASS(FELogFluidStressXZ,        FEELEMLOGDATA_ID, "fsxz");
REGISTER_FECORE_CLASS(FELogFluidRateOfDefXX,     FEELEMLOGDATA_ID, "fdxx");
REGISTER_FECORE_CLASS(FELogFluidRateOfDefYY,     FEELEMLOGDATA_ID, "fdyy");
REGISTER_FECORE_CLASS(FELogFluidRateOfDefZZ,     FEELEMLOGDATA_ID, "fdzz");
REGISTER_FECORE_CLASS(FELogFluidRateOfDefXY,     FEELEMLOGDATA_ID, "fdxy");
REGISTER_FECORE_CLASS(FELogFluidRateOfDefYZ,     FEELEMLOGDATA_ID, "fdyz");
REGISTER_FECORE_CLASS(FELogFluidRateOfDefXZ,     FEELEMLOGDATA_ID, "fdxz");

}
