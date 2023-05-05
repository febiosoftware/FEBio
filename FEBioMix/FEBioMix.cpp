/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FEBioMix.h"
#include "FEBiphasic.h"
#include "FEBiphasicSolute.h"
#include "FEMultiphasicStandard.h"
#include "FEMultiphasicMultigeneration.h"
#include "FESolute.h"
#include "FETriphasic.h"
#include "FEDiffConstIso.h"
#include "FEDiffConstOrtho.h"
#include "FEDiffRefIso.h"
#include "FEDiffAlbroIso.h"
#include "FEPermConstIso.h"
#include "FEPermHolmesMow.h"
#include "FEPermExpIso.h"
#include "FEPermRefIso.h"
#include "FEPermRefOrtho.h"
#include "FEPermRefTransIso.h"
#include "FEOsmCoefConst.h"
#include "FEOsmCoefManning.h"
#include "FESFDSBM.h"
#include "FEFiberExpPowSBM.h"
#include "FEFiberPowLinearSBM.h"
#include "FESolventSupplyStarling.h"
#include "FESolubConst.h"
#include "FESolubManning.h"
#include "FESupplyBinding.h"
#include "FESupplyConst.h"
#include "FESupplySynthesisBinding.h"
#include "FESupplyMichaelisMenten.h"
#include "FECarterHayes.h"
#include "FEReactionRateConst.h"
#include "FEReactionRateHuiskes.h"
#include "FEReactionRateNims.h"
#include "FEReactionRateExpSED.h"
#include "FEMembraneReactionRateConst.h"
#include "FEMembraneReactionRateIonChannel.h"
#include "FEMembraneReactionRateVoltageGated.h"
#include "FEMassActionForward.h"
#include "FEMassActionForwardEffective.h"
#include "FEMichaelisMenten.h"
#include "FEMassActionReversible.h"
#include "FEMassActionReversibleEffective.h"
#include "FEConcentrationIndependentReaction.h"
#include "FEMembraneMassActionForward.h"
#include "FEMembraneMassActionReversible.h"
#include "FEActiveConstantSupply.h"
#include "FEPorousNeoHookean.h"

#include "FEPoroTraction.h"
#include "FEFluidFlux.h"
#include "FESoluteFlux.h"
#include "FESoluteNaturalFlux.h"
#include "FEPressureStabilization.h"
#include "FEMatchingOsmoticCoefficientLoad.h"
#include "FEMatchingOsmoticCoefficientBC.h"
#include "FEMultiphasicFluidPressureLoad.h"

#include "FESlidingInterface2.h"
#include "FESlidingInterfaceBiphasic.h"
#include "FESlidingInterfaceBiphasicMixed.h"
#include "FESlidingInterface3.h"
#include "FESlidingInterfaceMP.h"
#include "FETiedBiphasicInterface.h"
#include "FETiedMultiphasicInterface.h"

#include "FEBiphasicSolver.h"
#include "FEBiphasicSoluteSolver.h"
#include "FEMultiphasicSolver.h"

#include "FEBioMixPlot.h"
#include "FEBioMixData.h"

#include "FEMixDomainFactory.h"
#include "FEBiphasicSolidDomain.h"
#include "FEBiphasicShellDomain.h"
#include "FEBiphasicSoluteSolidDomain.h"
#include "FEBiphasicSoluteShellDomain.h"
#include "FETriphasicDomain.h"
#include "FEMultiphasicSolidDomain.h"
#include "FEMultiphasicShellDomain.h"

#include "FESBMPointSource.h"
#include "FESolutePointSource.h"

#include "FEFixedFluidPressure.h"
#include "FEPrescribedNodalFluidPressure.h"
#include "FEFixedConcentration.h"
#include "FEPrescribedConcentration.h"

#include "FEInitialMixtureFluidPressure.h"
#include "FEInitialConcentration.h"
#include "FENodalFluidFlux.h"

#include "FEBiphasicModule.h"
#include "FEBiphasicAnalysis.h"
#include "FEBiphasicSoluteAnalysis.h"
#include "FEMultiphasicAnalysis.h"
#include <FECore/FEModelUpdate.h>

//-----------------------------------------------------------------------------
const char* FEBioMix::GetVariableName(FEBioMix::FEBIOMIX_VARIABLE var)
{
	switch (var)
	{
	case FLUID_PRESSURE: return "fluid pressure"; break;
	}

	assert(false);
	return nullptr;
}

//-----------------------------------------------------------------------------
//! Initialization of the FEBioMix module. This function registers all the classes
//! in this module with the FEBio framework.
void FEBioMix::InitModule()
{
//-----------------------------------------------------------------------------
// Domain factory
	FECoreKernel& febio = FECoreKernel::GetInstance();
	febio.RegisterDomain(new FEMixDomainFactory);

	// extensions to "solid" module
	febio.SetActiveModule("solid");
	REGISTER_FECORE_CLASS(FECarterHayes     , "Carter-Hayes");
	REGISTER_FECORE_CLASS(FEPorousNeoHookean, "porous neo-Hookean");
	REGISTER_FECORE_CLASS(FEPoroNormalTraction, "normal_traction");

//======================================================================
// setup the "biphasic" module
	febio.CreateModule(new FEBiphasicModule, "biphasic", 
		"{"
		"   \"title\" : \"Biphasic Analysis\","
		"   \"info\"  : \"Transient or quasi-static biphasic analysis.\""
		"}");

	febio.AddModuleDependency("solid");

	//-----------------------------------------------------------------------------
	// analyis classes (default type must match module name!)
	REGISTER_FECORE_CLASS(FEBiphasicAnalysis, "biphasic");

	//-----------------------------------------------------------------------------
	// solver classes
	REGISTER_FECORE_CLASS(FEBiphasicSolver      , "biphasic");

	//-----------------------------------------------------------------------------
	// Domain classes
	REGISTER_FECORE_CLASS(FEBiphasicSolidDomain, "biphasic-solid");
	REGISTER_FECORE_CLASS(FEBiphasicShellDomain, "biphasic-shell");

	//-----------------------------------------------------------------------------
	// Materials
	REGISTER_FECORE_CLASS(FEBiphasic             , "biphasic");
	REGISTER_FECORE_CLASS(FEPermConstIso         , "perm-const-iso");
	REGISTER_FECORE_CLASS(FEPermHolmesMow        , "perm-Holmes-Mow");
	REGISTER_FECORE_CLASS(FEPermExpIso           , "perm-exp-iso");
	REGISTER_FECORE_CLASS(FEPermRefIso           , "perm-ref-iso");
	REGISTER_FECORE_CLASS(FEPermRefOrtho         , "perm-ref-ortho");
	REGISTER_FECORE_CLASS(FEPermRefTransIso      , "perm-ref-trans-iso");
	REGISTER_FECORE_CLASS(FESolventSupplyStarling, "Starling");
	REGISTER_FECORE_CLASS(FEActiveConstantSupply , "active-const-supply");

	//-----------------------------------------------------------------------------
	// Boundary conditions
	REGISTER_FECORE_CLASS(FEFixedFluidPressure          , "zero fluid pressure");
	REGISTER_FECORE_CLASS(FEPrescribedNodalFluidPressure, "prescribed fluid pressure");
	REGISTER_FECORE_CLASS(FEPrescribedShellFluidPressure, "prescribed shell fluid pressure");

	//-----------------------------------------------------------------------------
	// Initial conditions
	REGISTER_FECORE_CLASS(FEInitialMixtureFluidPressure, "initial fluid pressure");
	REGISTER_FECORE_CLASS(FEInitialShellMixtureFluidPressure, "initial shell fluid pressure");
	REGISTER_FECORE_CLASS(FEInitialConcentration     , "initial concentration");
	REGISTER_FECORE_CLASS(FEInitialShellConcentration, "initial shell concentration");

	//-----------------------------------------------------------------------------
	// Nodal loads
	REGISTER_FECORE_CLASS(FENodalFluidFlux, "nodal fluidflux");

	//-----------------------------------------------------------------------------
	// Surface loads
	REGISTER_FECORE_CLASS(FEFluidFlux         , "fluidflux");
	REGISTER_FECORE_CLASS(FEPressureStabilization, "pressure_stabilization");

	//-----------------------------------------------------------------------------
	// Contact interfaces
	REGISTER_FECORE_CLASS(FESlidingInterface2            , "sliding2");
	REGISTER_FECORE_CLASS(FESlidingInterfaceBiphasic     , "sliding-biphasic");
	REGISTER_FECORE_CLASS(FESlidingInterfaceBiphasicMixed, "sliding-biphasic-mixed");
	REGISTER_FECORE_CLASS(FETiedBiphasicInterface        , "tied-biphasic");

	//-----------------------------------------------------------------------------
	// classes derived from FEPlotData
    REGISTER_FECORE_CLASS(FEPlotMPSpecificStrainEnergy           , "specific strain energy");
	REGISTER_FECORE_CLASS(FEPlotEffectiveElasticity		         , "effective elasticity"            );
	REGISTER_FECORE_CLASS(FEPlotEffectiveFluidPressure		     , "effective fluid pressure"        );
	REGISTER_FECORE_CLASS(FEPlotEffectiveShellFluidPressure      , "effective shell fluid pressure"  );
	REGISTER_FECORE_CLASS(FEPlotActualFluidPressure              , "fluid pressure"                  );
	REGISTER_FECORE_CLASS(FEPlotFluidFlux                        , "fluid flux"                      );
	REGISTER_FECORE_CLASS(FEPlotNodalFluidFlux                   , "nodal fluid flux");
	REGISTER_FECORE_CLASS(FEPlotPressureGap					     , "pressure gap"        );
	REGISTER_FECORE_CLASS(FEPlotFluidForce                       , "fluid force"         );
	REGISTER_FECORE_CLASS(FEPlotFluidForce2                      , "fluid force2"        );
	REGISTER_FECORE_CLASS(FEPlotFluidLoadSupport                 , "fluid load support"  );
	REGISTER_FECORE_CLASS(FEPlotMixtureFluidFlowRate             , "fluid flow rate"     );
	REGISTER_FECORE_CLASS(FEPlotReferentialSolidVolumeFraction   , "referential solid volume fraction");
	REGISTER_FECORE_CLASS(FEPlotPorosity                         , "porosity"            );
	REGISTER_FECORE_CLASS(FEPlotPerm                             , "permeability"        );
	REGISTER_FECORE_CLASS(FEPlotSolidStress                      , "solid stress"        );
	REGISTER_FECORE_CLASS(FEPlotLocalFluidLoadSupport            , "local fluid load support"        );
	REGISTER_FECORE_CLASS(FEPlotEffectiveFrictionCoeff           , "effective friction coefficient"  );

	//-----------------------------------------------------------------------------
	// Element log data
	REGISTER_FECORE_CLASS(FELogElemFluidPressure         , "p");
	REGISTER_FECORE_CLASS(FELogElemFluidFluxX            , "wx");
	REGISTER_FECORE_CLASS(FELogElemFluidFluxY            , "wy");
	REGISTER_FECORE_CLASS(FELogElemFluidFluxZ            , "wz");
	REGISTER_FECORE_CLASS(FELogElemPorosity              , "porosity");
	REGISTER_FECORE_CLASS_T(FELogElemPermeability_T, 0     , "Kpxx");
	REGISTER_FECORE_CLASS_T(FELogElemPermeability_T, 1     , "Kpyy");
	REGISTER_FECORE_CLASS_T(FELogElemPermeability_T, 2     , "Kpzz");
	REGISTER_FECORE_CLASS_T(FELogElemPermeability_T, 3     , "Kpxy");
	REGISTER_FECORE_CLASS_T(FELogElemPermeability_T, 4     , "Kpyz");
	REGISTER_FECORE_CLASS_T(FELogElemPermeability_T, 5     , "Kpxz");
	REGISTER_FECORE_CLASS_T(FELogElemSolidStress_T, 0, "esxx");
	REGISTER_FECORE_CLASS_T(FELogElemSolidStress_T, 1, "esyy");
	REGISTER_FECORE_CLASS_T(FELogElemSolidStress_T, 2, "eszz");
	REGISTER_FECORE_CLASS_T(FELogElemSolidStress_T, 3, "esxy");
	REGISTER_FECORE_CLASS_T(FELogElemSolidStress_T, 4, "esyz");
	REGISTER_FECORE_CLASS_T(FELogElemSolidStress_T, 5, "esxz");

	//-----------------------------------------------------------------------------
	// Model update requests (for biphasic module)
	febio.OnCreateEvent(UpdateModelWhenCreating<FEBiphasicAnalysis>([](FEModelUpdate& fem) {
		fem.AddPlotVariable("effective fluid pressure");
		fem.AddPlotVariable("fluid flux");
		})
	);


//======================================================================
// setup the "solute" module (i.e. biphasic-solute)
	febio.CreateModule(new FEBiphasicSoluteModule, "solute",
		"{"
		"   \"title\" : \"Biphasic Solute Analysis\","
		"   \"info\"  : \"Transient or quasi-static biphasic analysis with a single solute.\""
		"}");

	febio.AddModuleDependency("biphasic");

	//-----------------------------------------------------------------------------
	// Global data classes
	REGISTER_FECORE_CLASS(FESoluteData, "solute");

	//-----------------------------------------------------------------------------
	// analyis classes (default type must match module name!)
	REGISTER_FECORE_CLASS(FEBiphasicSoluteAnalysis, "solute");

	//-----------------------------------------------------------------------------
	// solver classes (default type must match module name!)
	REGISTER_FECORE_CLASS(FEBiphasicSoluteSolver, "solute");

	//-----------------------------------------------------------------------------
	// Domain classes
	REGISTER_FECORE_CLASS(FEBiphasicSoluteSolidDomain, "biphasic-solute-solid");
	REGISTER_FECORE_CLASS(FEBiphasicSoluteShellDomain, "biphasic-solute-shell");
	REGISTER_FECORE_CLASS(FETriphasicDomain, "triphasic-solid");

	//-----------------------------------------------------------------------------
	// Materials
	REGISTER_FECORE_CLASS(FEBiphasicSolute        , "biphasic-solute");
	REGISTER_FECORE_CLASS(FESoluteMaterial        , "solute");
	REGISTER_FECORE_CLASS(FETriphasic             , "triphasic");
	REGISTER_FECORE_CLASS(FEDiffConstIso          , "diff-const-iso");
	REGISTER_FECORE_CLASS(FEDiffConstOrtho        , "diff-const-ortho");
	REGISTER_FECORE_CLASS(FEDiffRefIso            , "diff-ref-iso");
	REGISTER_FECORE_CLASS(FEDiffAlbroIso          , "diff-Albro-iso");
	REGISTER_FECORE_CLASS(FEOsmCoefConst          , "osm-coef-const");
	REGISTER_FECORE_CLASS(FEOsmCoefManning        , "osm-coef-Manning");
	REGISTER_FECORE_CLASS(FESolubConst            , "solub-const");
	REGISTER_FECORE_CLASS(FESolubManning          , "solub-Manning");
	REGISTER_FECORE_CLASS(FESupplyBinding         , "supply-binding");
	REGISTER_FECORE_CLASS(FESupplyConst           , "supply-const");
	REGISTER_FECORE_CLASS(FESupplySynthesisBinding, "supply-synthesis-binding");
	REGISTER_FECORE_CLASS(FESupplyMichaelisMenten , "supply-Michaelis-Menten");

	//-----------------------------------------------------------------------------
	// Surface loads
	REGISTER_FECORE_CLASS(FESoluteFlux, "soluteflux");
    REGISTER_FECORE_CLASS(FESoluteNaturalFlux, "solute natural flux");
    REGISTER_FECORE_CLASS(FEMultiphasicFluidPressureLoad, "fluid pressure");

	//-----------------------------------------------------------------------------
	// boundary conditions
	REGISTER_FECORE_CLASS(FEFixedConcentration, "zero concentration");
	REGISTER_FECORE_CLASS(FEPrescribedConcentration, "prescribed concentration");

	//-----------------------------------------------------------------------------
	// Contact interfaces
	REGISTER_FECORE_CLASS(FESlidingInterface3, "sliding-biphasic-solute");

	//-----------------------------------------------------------------------------
	// classes derived from FEPlotData
	REGISTER_FECORE_CLASS(FEPlotEffectiveSoluteConcentration     , "effective solute concentration");
	REGISTER_FECORE_CLASS(FEPlotEffectiveShellSoluteConcentration, "effective shell solute concentration");
	REGISTER_FECORE_CLASS(FEPlotActualSoluteConcentration        , "solute concentration");
    REGISTER_FECORE_CLASS(FEPlotPartitionCoefficient             , "partition coefficient");
	REGISTER_FECORE_CLASS(FEPlotSoluteFlux		                 , "solute flux"                     );
    REGISTER_FECORE_CLASS(FEPlotSoluteVolumetricFlux             , "solute volumetric flux"          );
	REGISTER_FECORE_CLASS(FEPlotOsmolarity                       , "osmolarity");
	REGISTER_FECORE_CLASS(FEPlotCurrentDensity                   , "current density"     );
	REGISTER_FECORE_CLASS(FEPlotFixedChargeDensity               , "fixed charge density");
	REGISTER_FECORE_CLASS(FEPlotReferentialFixedChargeDensity    , "referential fixed charge density");
	REGISTER_FECORE_CLASS(FEPlotElectricPotential                , "electric potential"  );

	//-----------------------------------------------------------------------------
	// classes derived from FELogNodeData
	REGISTER_FECORE_CLASS(FENodeConcentration, "c");

	//-----------------------------------------------------------------------------
	// Element log data
	REGISTER_FECORE_CLASS(FELogElemSoluteConcentration   , "c");
	REGISTER_FECORE_CLASS(FELogElemSoluteFluxX           , "jx");
	REGISTER_FECORE_CLASS(FELogElemSoluteFluxY           , "jy");
	REGISTER_FECORE_CLASS(FELogElemSoluteFluxZ           , "jz");
	REGISTER_FECORE_CLASS(FELogElemSoluteRefConcentration, "crc");
	REGISTER_FECORE_CLASS(FELogFixedChargeDensity        , "cF");

	REGISTER_FECORE_CLASS_T(FELogElemSoluteConcentration_T, 0, "c1");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteConcentration_T, 1, "c2");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteConcentration_T, 2, "c3");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteConcentration_T, 3, "c4");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteConcentration_T, 4, "c5");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteConcentration_T, 5, "c6");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteConcentration_T, 6, "c7");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteConcentration_T, 7, "c8");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxX_T, 0, "j1x");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxY_T, 0, "j1y");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxZ_T, 0, "j1z");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxX_T, 1, "j2x");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxY_T, 1, "j2y");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxZ_T, 1, "j2z");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxX_T, 2, "j3x");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxY_T, 2, "j3y");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxZ_T, 2, "j3z");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxX_T, 3, "j4x");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxY_T, 3, "j4y");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxZ_T, 3, "j4z");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxX_T, 4, "j5x");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxY_T, 4, "j5y");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxZ_T, 4, "j5z");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxX_T, 5, "j6x");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxY_T, 5, "j6y");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxZ_T, 5, "j6z");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxX_T, 6, "j7x");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxY_T, 6, "j7y");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxZ_T, 6, "j7z");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxX_T, 7, "j8x");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxY_T, 7, "j8y");
	REGISTER_FECORE_CLASS_T(FELogElemSoluteFluxZ_T, 7, "j8z");

//======================================================================
// setup the "multiphasic" module
	febio.CreateModule(new FEMultiphasicModule, "multiphasic",
		"{"
		"   \"title\" : \"Multiphasic Analysis\","
		"   \"info\"  : \"Transient or quasi-static analysis with solutes.\""
		"}");

	febio.AddModuleDependency("solute");

	//-----------------------------------------------------------------------------
	// Global data classes
	REGISTER_FECORE_CLASS(FESBMData, "solid_bound");

	//-----------------------------------------------------------------------------
	// analyis classes (default type must match module name!)
	REGISTER_FECORE_CLASS(FEMultiphasicAnalysis, "multiphasic");

	//-----------------------------------------------------------------------------
	// solver classes
	REGISTER_FECORE_CLASS(FEMultiphasicSolver   , "multiphasic");

	//-----------------------------------------------------------------------------
	// Domain classes
	REGISTER_FECORE_CLASS(FEMultiphasicSolidDomain, "multiphasic-solid");
	REGISTER_FECORE_CLASS(FEMultiphasicShellDomain, "multiphasic-shell");

	//-----------------------------------------------------------------------------
	// Materials
	REGISTER_FECORE_CLASS(FEMultiphasicStandard               , "multiphasic"       );
	REGISTER_FECORE_CLASS(FEMultiphasicMultigeneration        , "multiphasic-multigeneration");
	REGISTER_FECORE_CLASS(FESFDSBM                            , "spherical fiber distribution sbm");
	REGISTER_FECORE_CLASS(FEFiberExpPowSBM                    , "fiber-exp-pow sbm" );
	REGISTER_FECORE_CLASS(FEFiberPowLinearSBM                 , "fiber-pow-linear sbm");
	REGISTER_FECORE_CLASS(FEReactionRateConst		    	  , "constant reaction rate"    );
	REGISTER_FECORE_CLASS(FEReactionRateHuiskes		    	  , "Huiskes reaction rate"     );
	REGISTER_FECORE_CLASS(FEReactionRateNims		    	  , "Nims reaction rate"        );
	REGISTER_FECORE_CLASS(FEReactionRateExpSED                , "exp-sed reaction rate"     );
	REGISTER_FECORE_CLASS(FEMembraneReactionRateConst         , "membrane constant reaction rate");
	REGISTER_FECORE_CLASS(FEMembraneReactionRateIonChannel    , "membrane ion channel reaction rate");
	REGISTER_FECORE_CLASS(FEMembraneReactionRateVoltageGated  , "membrane voltage-gated reaction rate");
	REGISTER_FECORE_CLASS(FEMassActionForward		    	  , "mass-action-forward"       );
	REGISTER_FECORE_CLASS(FEMassActionForwardEffective		  , "mass-action-forward-effective");
	REGISTER_FECORE_CLASS(FEMembraneMassActionForward         , "membrane-mass-action-forward");
	REGISTER_FECORE_CLASS(FEConcentrationIndependentReaction  , "concentration-independent");
	REGISTER_FECORE_CLASS(FEMassActionReversible              , "mass-action-reversible"   );
	REGISTER_FECORE_CLASS(FEMassActionReversibleEffective     , "mass-action-reversible-effective");
	REGISTER_FECORE_CLASS(FEMembraneMassActionReversible      , "membrane-mass-action-reversible");
	REGISTER_FECORE_CLASS(FEMichaelisMenten                   , "Michaelis-Menten"         );
	REGISTER_FECORE_CLASS(FESolidBoundMolecule                , "solid_bound"              );

	REGISTER_FECORE_CLASS(FEReactantSpeciesRef, "vR");
	REGISTER_FECORE_CLASS(FEProductSpeciesRef , "vP");
	REGISTER_FECORE_CLASS(FEInternalReactantSpeciesRef, "vRi");
	REGISTER_FECORE_CLASS(FEInternalProductSpeciesRef , "vPi");
	REGISTER_FECORE_CLASS(FEExternalReactantSpeciesRef, "vRe");
	REGISTER_FECORE_CLASS(FEExternalProductSpeciesRef , "vPe");
    
	//-----------------------------------------------------------------------------
	// Surface loads
	REGISTER_FECORE_CLASS(FEMatchingOsmoticCoefficientLoad, "matching_osm_coef", 0x0300); // deprecated, use BC version

	//-----------------------------------------------------------------------------
	// Boundary conditions
	REGISTER_FECORE_CLASS(FEMatchingOsmoticCoefficientBC, "matching_osm_coef");

	//-----------------------------------------------------------------------------
	// Body loads
	REGISTER_FECORE_CLASS(FESBMPointSource   , "sbm point source");
	REGISTER_FECORE_CLASS(FESolutePointSource, "solute point source");

	//-----------------------------------------------------------------------------
	// Contact interfaces
	REGISTER_FECORE_CLASS(FESlidingInterfaceMP      , "sliding-multiphasic"    );
	REGISTER_FECORE_CLASS(FEAmbientConcentration    , "ambient_concentration"  );
	REGISTER_FECORE_CLASS(FETiedMultiphasicInterface, "tied-multiphasic"       );

	//-----------------------------------------------------------------------------
	// classes derived from FEPlotData
	REGISTER_FECORE_CLASS(FEPlotReceptorLigandConcentration      , "receptor-ligand concentration"   );
	REGISTER_FECORE_CLASS(FEPlotSBMConcentration                 , "sbm concentration"			     );
    REGISTER_FECORE_CLASS(FEPlotSBMArealConcentration            , "sbm areal concentration"         );
	REGISTER_FECORE_CLASS(FEPlotSBMRefAppDensity			     , "sbm referential apparent density");
	REGISTER_FECORE_CLASS(FEPlotOsmoticCoefficient               , "osmotic coefficient"             );

	//-----------------------------------------------------------------------------
	// Element log data
	REGISTER_FECORE_CLASS(FELogElemElectricPotential     , "psi");
	REGISTER_FECORE_CLASS(FELogElemCurrentDensityX       , "Iex");
	REGISTER_FECORE_CLASS(FELogElemCurrentDensityY       , "Iey");
	REGISTER_FECORE_CLASS(FELogElemCurrentDensityZ       , "Iez");

	REGISTER_FECORE_CLASS_T(FELogElemSBMConcentration_T, 0, "sbm1");
	REGISTER_FECORE_CLASS_T(FELogElemSBMConcentration_T, 1, "sbm2");
	REGISTER_FECORE_CLASS_T(FELogElemSBMConcentration_T, 2, "sbm3");
	REGISTER_FECORE_CLASS_T(FELogElemSBMConcentration_T, 3, "sbm4");
	REGISTER_FECORE_CLASS_T(FELogElemSBMConcentration_T, 4, "sbm5");
	REGISTER_FECORE_CLASS_T(FELogElemSBMConcentration_T, 5, "sbm6");
	REGISTER_FECORE_CLASS_T(FELogElemSBMConcentration_T, 6, "sbm7");
	REGISTER_FECORE_CLASS_T(FELogElemSBMConcentration_T, 7, "sbm8");

	//-----------------------------------------------------------------------------
	// domain log data
	REGISTER_FECORE_CLASS_T(FELogDomainIntegralSBMConcentration_T, 0, "sbm1_integral");
	REGISTER_FECORE_CLASS_T(FELogDomainIntegralSBMConcentration_T, 1, "sbm2_integral");
	REGISTER_FECORE_CLASS_T(FELogDomainIntegralSBMConcentration_T, 2, "sbm3_integral");
	REGISTER_FECORE_CLASS_T(FELogDomainIntegralSBMConcentration_T, 3, "sbm4_integral");
	REGISTER_FECORE_CLASS_T(FELogDomainIntegralSBMConcentration_T, 4, "sbm5_integral");
	REGISTER_FECORE_CLASS_T(FELogDomainIntegralSBMConcentration_T, 5, "sbm6_integral");
	REGISTER_FECORE_CLASS_T(FELogDomainIntegralSBMConcentration_T, 6, "sbm7_integral");
	REGISTER_FECORE_CLASS_T(FELogDomainIntegralSBMConcentration_T, 7, "sbm8_integral");

	REGISTER_FECORE_CLASS_T(FELogDomainIntegralSoluteConcentration_T, 0, "c1_integral");
	REGISTER_FECORE_CLASS_T(FELogDomainIntegralSoluteConcentration_T, 1, "c2_integral");
	REGISTER_FECORE_CLASS_T(FELogDomainIntegralSoluteConcentration_T, 2, "c3_integral");
	REGISTER_FECORE_CLASS_T(FELogDomainIntegralSoluteConcentration_T, 3, "c4_integral");
	REGISTER_FECORE_CLASS_T(FELogDomainIntegralSoluteConcentration_T, 4, "c5_integral");
	REGISTER_FECORE_CLASS_T(FELogDomainIntegralSoluteConcentration_T, 5, "c6_integral");
	REGISTER_FECORE_CLASS_T(FELogDomainIntegralSoluteConcentration_T, 6, "c7_integral");
	REGISTER_FECORE_CLASS_T(FELogDomainIntegralSoluteConcentration_T, 7, "c8_integral");

	febio.SetActiveModule(0);
}
