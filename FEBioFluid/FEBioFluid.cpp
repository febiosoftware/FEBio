#include "FEBioFluid.h"
#include "FEFluid.h"
#include "FENewtonianFluid.h"
#include "FECarreauFluid.h"
#include "FECarreauYasudaFluid.h"
#include "FEPowellEyringFluid.h"
#include "FECrossFluid.h"
#include "FEFluidFSI.h"

#include "FEFluidSolver.h"
#include "FEFluidDomain3D.h"
#include "FEFluidDomain2D.h"

#include "FEFluidTractionLoad.h"
#include "FEFluidNormalTraction.h"
#include "FEFluidNormalVelocity.h"
#include "FEFluidVelocity.h"
#include "FEFluidRotationalVelocity.h"
#include "FEFluidResistanceBC.h"
#include "FETangentialDamping.h"
#include "FETangentialFlowStabilization.h"
#include "FEBackFlowStabilization.h"

#include "FEConstraintFrictionlessWall.h"
#include "FEConstraintNormalFlow.h"

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

	// define the fluid module
	febio.CreateModule("fluid");

//-----------------------------------------------------------------------------
// solver classes
REGISTER_FECORE_CLASS(FEFluidSolver        , FESOLVER_ID, "fluid"         );

//-----------------------------------------------------------------------------
// Materials
REGISTER_FECORE_CLASS(FEFluid                       ,FEMATERIAL_ID, "fluid"         );
REGISTER_FECORE_CLASS(FENewtonianFluid              ,FEMATERIAL_ID, "Newtonian fluid");
REGISTER_FECORE_CLASS(FECarreauFluid                ,FEMATERIAL_ID, "Carreau"       );
REGISTER_FECORE_CLASS(FECarreauYasudaFluid          ,FEMATERIAL_ID, "Carreau-Yasuda");
REGISTER_FECORE_CLASS(FEPowellEyringFluid           ,FEMATERIAL_ID, "Powell-Eyring" );
REGISTER_FECORE_CLASS(FECrossFluid                  ,FEMATERIAL_ID, "Cross"         );
    
//-----------------------------------------------------------------------------
// Domain classes
REGISTER_FECORE_CLASS(FEFluidDomain3D               , FEDOMAIN_ID, "fluid-3D"              );
REGISTER_FECORE_CLASS(FEFluidDomain2D               , FEDOMAIN_ID, "fluid-2D"              );

//-----------------------------------------------------------------------------
// Surface loads
REGISTER_FECORE_CLASS(FEFluidTractionLoad           , FESURFACELOAD_ID, "fluid viscous traction");
REGISTER_FECORE_CLASS(FEFluidNormalTraction         , FESURFACELOAD_ID, "fluid normal traction");
REGISTER_FECORE_CLASS(FEFluidNormalVelocity         , FESURFACELOAD_ID, "fluid normal velocity");
REGISTER_FECORE_CLASS(FEFluidVelocity               , FESURFACELOAD_ID, "fluid velocity");
REGISTER_FECORE_CLASS(FEFluidRotationalVelocity     , FESURFACELOAD_ID, "fluid rotational velocity");
REGISTER_FECORE_CLASS(FEFluidResistanceBC           , FESURFACELOAD_ID, "fluid resistance");
REGISTER_FECORE_CLASS(FETangentialDamping           , FESURFACELOAD_ID, "fluid tangential damping");
REGISTER_FECORE_CLASS(FETangentialFlowStabilization , FESURFACELOAD_ID, "fluid tangential stabilization");
REGISTER_FECORE_CLASS(FEBackFlowStabilization       , FESURFACELOAD_ID, "fluid backflow stabilization");
    
//-----------------------------------------------------------------------------
// constraint classes
REGISTER_FECORE_CLASS(FEConstraintFrictionlessWall  , FENLCONSTRAINT_ID, "frictionless fluid wall");
REGISTER_FECORE_CLASS(FEConstraintNormalFlow        , FENLCONSTRAINT_ID, "normal fluid flow"      );

//-----------------------------------------------------------------------------
// classes derived from FEPlotData
REGISTER_FECORE_CLASS(FEPlotDisplacement            , FEPLOTDATA_ID, "displacement"             );
REGISTER_FECORE_CLASS(FEPlotNodalFluidVelocity      , FEPLOTDATA_ID, "nodal fluid velocity"     );
REGISTER_FECORE_CLASS(FEPlotNodalRelativeFluidVelocity, FEPLOTDATA_ID, "nodal relative fluid velocity"     );
REGISTER_FECORE_CLASS(FEPlotFluidDilatation         , FEPLOTDATA_ID, "fluid dilatation"         );
REGISTER_FECORE_CLASS(FEPlotElasticFluidPressure	, FEPLOTDATA_ID, "elastic fluid pressure"   );
REGISTER_FECORE_CLASS(FEPlotFluidVolumeRatio		, FEPLOTDATA_ID, "fluid volume ratio"       );
REGISTER_FECORE_CLASS(FEPlotFluidDensity            , FEPLOTDATA_ID, "fluid density"            );
REGISTER_FECORE_CLASS(FEPlotFluidDensityRate        , FEPLOTDATA_ID, "fluid density rate"       );
REGISTER_FECORE_CLASS(FEPlotFluidVelocity           , FEPLOTDATA_ID, "fluid velocity"           );
REGISTER_FECORE_CLASS(FEPlotRelativeFluidVelocity   , FEPLOTDATA_ID, "relative fluid velocity"  );
REGISTER_FECORE_CLASS(FEPlotFluidAcceleration       , FEPLOTDATA_ID, "fluid acceleration"       );
REGISTER_FECORE_CLASS(FEPlotFluidVorticity          , FEPLOTDATA_ID, "fluid vorticity"          );
REGISTER_FECORE_CLASS(FEPlotElementFluidStress      , FEPLOTDATA_ID, "fluid stress"             );
REGISTER_FECORE_CLASS(FEPlotElementFluidRateOfDef   , FEPLOTDATA_ID, "fluid rate of deformation");
REGISTER_FECORE_CLASS(FEPlotFluidStressPowerDensity , FEPLOTDATA_ID, "fluid stress power density");
REGISTER_FECORE_CLASS(FEPlotFluidHeatSupplyDensity  , FEPLOTDATA_ID, "fluid heat supply density");
REGISTER_FECORE_CLASS(FEPlotFluidSurfaceForce       , FEPLOTDATA_ID, "fluid surface force"      );
REGISTER_FECORE_CLASS(FEPlotFluidSurfaceTractionPower, FEPLOTDATA_ID, "fluid surface traction power");
REGISTER_FECORE_CLASS(FEPlotFluidSurfaceEnergyFlux  , FEPLOTDATA_ID, "fluid surface energy flux");
REGISTER_FECORE_CLASS(FEPlotFluidShearViscosity     , FEPLOTDATA_ID, "fluid shear viscosity"    );
REGISTER_FECORE_CLASS(FEPlotFluidMassFlowRate       , FEPLOTDATA_ID, "fluid mass flow rate"     );
REGISTER_FECORE_CLASS(FEPlotFluidStrainEnergyDensity, FEPLOTDATA_ID, "fluid strain energy density");
REGISTER_FECORE_CLASS(FEPlotFluidKineticEnergyDensity, FEPLOTDATA_ID,"fluid kinetic energy density");
REGISTER_FECORE_CLASS(FEPlotFluidEnergyDensity      , FEPLOTDATA_ID, "fluid energy density"     );
REGISTER_FECORE_CLASS(FEPlotFluidElementStrainEnergy, FEPLOTDATA_ID, "fluid element strain energy");
REGISTER_FECORE_CLASS(FEPlotFluidElementKineticEnergy, FEPLOTDATA_ID, "fluid element kinetic energy");
REGISTER_FECORE_CLASS(FEPlotFluidElementLinearMomentum, FEPLOTDATA_ID, "fluid element linear momentum");
REGISTER_FECORE_CLASS(FEPlotFluidElementAngularMomentum, FEPLOTDATA_ID, "fluid element angular momentum");
REGISTER_FECORE_CLASS(FEPlotFluidElementCenterOfMass, FEPLOTDATA_ID, "fluid element center of mass");
REGISTER_FECORE_CLASS(FEPlotFluidFlowRate           , FEPLOTDATA_ID, "fluid flow rate");
REGISTER_FECORE_CLASS(FEPlotFluidPressure           , FEPLOTDATA_ID, "fluid pressure");

//-----------------------------------------------------------------------------
REGISTER_FECORE_CLASS(FENodeFluidXVel,           FEELEMLOGDATA_ID, "nfvx");
REGISTER_FECORE_CLASS(FENodeFluidYVel,           FEELEMLOGDATA_ID, "nfvy");
REGISTER_FECORE_CLASS(FENodeFluidZVel,           FEELEMLOGDATA_ID, "nfvz");
REGISTER_FECORE_CLASS(FELogElemFluidPosX,        FEELEMLOGDATA_ID, "fx");
REGISTER_FECORE_CLASS(FELogElemFluidPosY,        FEELEMLOGDATA_ID, "fy");
REGISTER_FECORE_CLASS(FELogElemFluidPosZ,        FEELEMLOGDATA_ID, "fz");
REGISTER_FECORE_CLASS(FELogElasticFluidPressure, FEELEMLOGDATA_ID, "fp");
REGISTER_FECORE_CLASS(FELogFluidVolumeRatio,     FEELEMLOGDATA_ID, "fJ");
REGISTER_FECORE_CLASS(FELogFluidDensity,         FEELEMLOGDATA_ID, "fd");
REGISTER_FECORE_CLASS(FELogFluidStressPower,     FEELEMLOGDATA_ID, "fsp");
REGISTER_FECORE_CLASS(FELogFluidVelocityX,       FEELEMLOGDATA_ID, "fvx");
REGISTER_FECORE_CLASS(FELogFluidVelocityY,       FEELEMLOGDATA_ID, "fvy");
REGISTER_FECORE_CLASS(FELogFluidVelocityZ,       FEELEMLOGDATA_ID, "fvz");
REGISTER_FECORE_CLASS(FELogFluidAccelerationX,   FEELEMLOGDATA_ID, "fax");
REGISTER_FECORE_CLASS(FELogFluidAccelerationY,   FEELEMLOGDATA_ID, "fay");
REGISTER_FECORE_CLASS(FELogFluidAccelerationZ,   FEELEMLOGDATA_ID, "faz");
REGISTER_FECORE_CLASS(FELogFluidVorticityX,      FEELEMLOGDATA_ID, "fwx");
REGISTER_FECORE_CLASS(FELogFluidVorticityY,      FEELEMLOGDATA_ID, "fwy");
REGISTER_FECORE_CLASS(FELogFluidVorticityZ,      FEELEMLOGDATA_ID, "fwz");
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

	febio.SetActiveModule(0);
}
