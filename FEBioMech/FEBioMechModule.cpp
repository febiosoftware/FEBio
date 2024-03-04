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

#include "FEBioMech.h"
#include "FEBioMechModule.h"
#include "FE2DFiberNeoHookean.h"
#include "FE2DTransIsoMooneyRivlin.h"
#include "FE2DTransIsoVerondaWestmann.h"
#include "FEABUnconstrained.h"
#include "FEActiveFiberContraction.h"
#include "FEArrudaBoyce.h"
#include "FECarreauYasudaViscousSolid.h"
#include "FECarterHayesOld.h"
#include "FECellGrowth.h"
#include "FECubicCLE.h"
#include "FEDamageMooneyRivlin.h"
#include "FEDamageNeoHookean.h"
#include "FEDamageTransIsoMooneyRivlin.h"
#include "FESpringMaterial.h"
#include "FETorsionalSpring.h"
#include "FEDonnanEquilibrium.h"
#include "FEEFDDonnanEquilibrium.h"
#include "FEEFDMooneyRivlin.h"
#include "FEEFDNeoHookean.h"
#include "FEEFDUncoupled.h"
#include "FEEFDVerondaWestmann.h"
#include "FEElasticMixture.h"
#include "FEElasticMultigeneration.h"
#include "FEEllipsoidalFiberDistribution.h"
#include "FEFiberCDF.h"
#include "FEFiberCDFUncoupled.h"
#include "FEFiberEntropyChain.h"
#include "FEFiberEntropyChainUC.h"
#include "FEFiberExpPow.h"
#include "FEFiberExpPowUncoupled.h"
#include "FEFiberNaturalNeoHookean.h"
#include "FEFiberNeoHookean.h"
#include "FEFiberPowLinear.h"
#include "FEFiberPowLinearUncoupled.h"
#include "FEFiberEFDNeoHookean.h"
#include "FEFiberExponentialPowerUC.h"
#include "FEFiberNHUC.h"
#include "FEFiberKiousisUncoupled.h"
#include "FEForceVelocityContraction.h"
#include "FEFungOrthoCompressible.h"
#include "FEFungOrthotropic.h"
#include "FEHolmesMow.h"
#include "FEHolmesMowUC.h"
#include "FEHolzapfelGasserOgden.h"
#include "FEHolzapfelUnconstrained.h"
#include "FEHuiskesSupply.h"
#include "FEIncompNeoHookean.h"
#include "FEIsotropicElastic.h"
#include "FEMooneyRivlin.h"
#include "FEMRVonMisesFibers.h"
#include "FEMuscleMaterial.h"
#include "FENaturalNeoHookean.h"
#include "FENeoHookean.h"
#include "FENeoHookeanTransIso.h"
#include "FENewtonianViscousSolid.h"
#include "FENewtonianViscousSolidUC.h"
#include "FEOgdenMaterial.h"
#include "FEOgdenUnconstrained.h"
#include "FEOrthoElastic.h"
#include "FEOrthotropicCLE.h"
#include "FEOsmoticVirialExpansion.h"
#include "FEPerfectOsmometer.h"
#include "FEReactiveFatigue.h"
#include "FEUncoupledReactiveFatigue.h"
#include "FERemodelingElasticMaterial.h"
#include "FERigidMaterial.h"
#include "FESphericalFiberDistribution.h"
#include "FEStVenantKirchhoff.h"
#include "FETCNonlinearOrthotropic.h"
#include "FETendonMaterial.h"
#include "FETraceFreeNeoHookean.h"
#include "FETransIsoMooneyRivlin.h"
#include "FETransIsoMREstrada.h"
#include "FETransIsoVerondaWestmann.h"
#include "FETrussMaterial.h"
#include "FEUncoupledActiveContraction.h"
#include "FEUncoupledElasticMixture.h"
#include "FEUncoupledViscoElasticMaterial.h"
#include "FEUncoupledViscoElasticDamage.h"
#include "FEVerondaWestmann.h"
#include "FEViscoElasticDamage.h"
#include "FEViscoElasticMaterial.h"
#include "FEVonMisesPlasticity.h"
#include "FEElasticFiberMaterial.h"
#include "FEElasticFiberMaterialUC.h"
#include "FEFiberDensityDistribution.h"
#include "FEContinuousFiberDistribution.h"
#include "FEContinuousFiberDistributionUC.h"
#include "FEFiberIntegrationGauss.h"
#include "FEFiberIntegrationTrapezoidal.h"
#include "FEFiberIntegrationGeodesic.h"
#include "FEFiberIntegrationGaussKronrod.h"
#include "FEFiberIntegrationTriangle.h"
#include "FECoupledTransIsoMooneyRivlin.h"
#include "FECoupledTransIsoVerondaWestmann.h"
#include "FEHGOCoronary.h"
#include "FENonlinearSpring.h"
#include "FEDiscreteElementMaterial.h"
#include "FEPRLig.h"
#include "FECoupledMooneyRivlin.h"
#include "FECoupledVerondaWestmann.h"
#include "FEReactivePlasticity.h"
#include "FEReactivePlasticDamage.h"
#include "FEReactiveViscoelastic.h"
#include "FEUncoupledReactiveViscoelastic.h"
#include "FEBondRelaxation.h"
#include "FEBondRecruitment.h"
#include "FEDamageMaterial.h"
#include "FEDamageMaterialUC.h"
#include "FERVEDamageMaterial.h"
#include "FERVEFatigueMaterial.h"
#include "FEDamageCDF.h"
#include "FEDamageCriterion.h"
#include "FEPlasticFlowCurve.h"
#include "FEFiberExpLinear.h"
#include "FEUncoupledFiberExpLinear.h"
#include "FEPrescribedActiveContractionUniaxial.h"
#include "FEPrescribedActiveContractionUniaxialUC.h"
#include "FEPrescribedActiveContractionTransIso.h"
#include "FEPrescribedActiveContractionTransIsoUC.h"
#include "FEPrescribedActiveContractionIsotropic.h"
#include "FEPrescribedActiveContractionIsotropicUC.h"
#include "FEGentMaterial.h"
#include "FEWrinkleOgdenMaterial.h"
#include "FEGenericHyperelastic.h"
#include "FEGenericHyperelasticUC.h"
#include "FEGenericTransIsoHyperelastic.h"
#include "FEGenericTransIsoHyperelasticUC.h"
#include "FEActiveFiberStress.h"
#include "FEActiveFiberStressUC.h"
#include "FEContinuousElasticDamage.h"
#include "FEIsotropicLeeSacks.h"
#include "FEIsotropicLeeSacksUncoupled.h"
#include "FEPolynomialHyperElastic.h"
#include "FEShenoyMaterial.h"
#include "FELungMaterial.h"
#include "FEGrowthTensor.h"
#include "FEKinematicGrowth.h"
#include "FEYeoh.h"

#include "FEPressureLoad.h"
#include "FETractionLoad.h"
#include "FESurfaceForceUniform.h"
#include "FEBearingLoad.h"
#include "FEGenericBodyForce.h"
#include "FECentrifugalBodyForce.h"
#include "FEPointBodyForce.h"
#include "FESurfaceAttractionBodyForce.h"
#include "FEMassDamping.h"

#include "FEFacet2FacetSliding.h"
#include "FEPeriodicBoundary.h"
#include "FERigidWallInterface.h"
#include "FESlidingInterface.h"
#include "FESlidingElasticInterface.h"
#include "FEPeriodicSurfaceConstraint.h"
#include "FETiedInterface.h"
#include "FETiedElasticInterface.h"
#include "FEStickyInterface.h"
#include "FEPointConstraint.h"
#include "FEAzimuthConstraint.h"
#include "FEFacet2FacetTied.h"
#include "FEVolumeConstraint.h"
#include "FEDistanceConstraint.h"
#include "FEMortarSlidingContact.h"
#include "FEMortarTiedContact.h"
#include "FEContactPotential.h"

#include "FESymmetryPlane.h"
#include "FERigidJoint.h"
#include "FEGenericRigidJoint.h"
#include "FERigidSphericalJoint.h"
#include "FERigidRevoluteJoint.h"
#include "FERigidPrismaticJoint.h"
#include "FERigidCylindricalJoint.h"
#include "FERigidPlanarJoint.h"
#include "FERigidLock.h"
#include "FERigidSpring.h"
#include "FERigidDamper.h"
#include "FERigidAngularDamper.h"
#include "FERigidContractileForce.h"
#include "FERigidForce.h"
#include "FERigidCable.h"
#include "FEDiscreteContact.h"
#include "FERigidFollowerForce.h"
#include "FERigidFollowerMoment.h"
#include "FEFixedNormalDisplacement.h"

#include "FESolidSolver.h"
#include "FESolidSolver2.h"
#include "FEExplicitSolidSolver.h"
#include "FECGSolidSolver.h"

#include "FEBioMechPlot.h"
#include "FEBioMechData.h"

#include "FESolidDomainFactory.h"
#include "FEElasticSolidDomain.h"
#include "FEElasticShellDomain.h"
#include "FEElasticShellDomainOld.h"
#include "FEElasticEASShellDomain.h"
#include "FEElasticANSShellDomain.h"
#include "FEElasticTrussDomain.h"
#include "FELinearTrussDomain.h"
#include "FERigidSolidDomain.h"
#include "FERigidShellDomain.h"
#include "FEElasticBeamDomain.h"
#include "FERemodelingElasticDomain.h"
#include "FEUDGHexDomain.h"
#include "FEUT4Domain.h"
#include "FESRIElasticSolidDomain.h"
#include "FE3FieldElasticSolidDomain.h"
#include "FE3FieldElasticShellDomain.h"
#include "FEDiscreteElasticDomain.h"
#include "FEDeformableSpringDomain.h"
#include "RigidBC.h"
#include "FERigidNodeSet.h"
#include "FEFixedDisplacement.h"
#include "FEFixedShellDisplacement.h"
#include "FEFixedRotation.h"
#include "FEPrescribedDisplacement.h"
#include "FEPrescribedShellDisplacement.h"
#include "FEPrescribedRotation.h"
#include "FEBCPrescribedDeformation.h"
#include "FEBCRigidDeformation.h"
#include "FEPrescribedNormalDisplacement.h"
#include "FEMaxStressCriterion.h"
#include "FEMaxDamageCriterion.h"
#include "FESpringRuptureCriterion.h"

#include "FEInitialVelocity.h"
#include "FENodalForce.h"

#include "FEPreStrainElastic.h"
#include "FEPreStrainUncoupledElastic.h"
#include "FEConstPrestrain.h"
#include "FEInSituStretchGradient.h"
#include "FEPreStrainConstraint.h"
#include "FEInitialPreStrain.h"

#include "FENodeToNodeConstraint.h"

#include "FEDeformationMapGenerator.h"

#include "FESolidModule.h"

#include "FESolidAnalysis.h"
#include <FECore/FEModelUpdate.h>

#include "FEElasticBeamMaterial.h"

//-----------------------------------------------------------------------------
//! Register all the classes of the FEBioMech module with the FEBio framework.
void FEBioMech::InitModule()
{
	//-----------------------------------------------------------------------------
	// Domain factory
	FECoreKernel& febio = FECoreKernel::GetInstance();
	febio.RegisterDomain(new FESolidDomainFactory);

	//-----------------------------------------------------------------------------
	// create module
	febio.CreateModule(new FESolidModule, "solid",
		"{"
		"   \"title\" : \"Structural Mechanics\","
		"   \"info\"  : \"Quasi-static or dynamical structural mechanics analysis.\""
		"}");

	//-----------------------------------------------------------------------------
	// analyis classes (default type must match module name!)
	REGISTER_FECORE_CLASS(FESolidAnalysis, "solid");

	//-----------------------------------------------------------------------------
	// Solver classes (default type must match module name!)
	REGISTER_FECORE_CLASS(FESolidSolver2, "solid");
	REGISTER_FECORE_CLASS(FEExplicitSolidSolver, "explicit-solid");
	REGISTER_FECORE_CLASS(FESolidSolver, "solid_old");
	REGISTER_FECORE_CLASS(FECGSolidSolver, "CG-solid");

	//-----------------------------------------------------------------------------
	// material classes

	// elastic materials (derived from FEElasticMaterial)
	REGISTER_FECORE_CLASS(FE2DFiberNeoHookean, "2D fiber neo-Hookean");
    REGISTER_FECORE_CLASS(FEABUnconstrained, "Arruda-Boyce unconstrained");
	REGISTER_FECORE_CLASS(FECarreauYasudaViscousSolid, "Carreau-Yasuda viscous solid");
	REGISTER_FECORE_CLASS(FECellGrowth, "cell growth");
	REGISTER_FECORE_CLASS(FECubicCLE, "cubic CLE");
	REGISTER_FECORE_CLASS(FEDamageNeoHookean, "damage neo-Hookean");
	REGISTER_FECORE_CLASS(FEDonnanEquilibrium, "Donnan equilibrium");
	REGISTER_FECORE_CLASS(FEEFDDonnanEquilibrium, "EFD Donnan equilibrium");
	REGISTER_FECORE_CLASS(FEEFDNeoHookean, "EFD neo-Hookean (new)");
	REGISTER_FECORE_CLASS(FEEFDNeoHookeanOld, "EFD neo-Hookean");
	REGISTER_FECORE_CLASS(FEEllipsoidalFiberDistribution, "ellipsoidal fiber distribution");
	REGISTER_FECORE_CLASS(FEEllipsoidalFiberDistributionOld, "ellipsoidal fiber distribution (old)");
	REGISTER_FECORE_CLASS(FEFungOrthoCompressible, "Fung-ortho-compressible");
	REGISTER_FECORE_CLASS(FECompressibleGentMaterial, "compressible Gent");
	REGISTER_FECORE_CLASS(FEHolmesMow, "Holmes-Mow");
    REGISTER_FECORE_CLASS(FEHolmesMowUC, "uncoupled Holmes-Mow");
    REGISTER_FECORE_CLASS(FEHolzapfelUnconstrained, "HGO unconstrained");
	REGISTER_FECORE_CLASS(FEIsotropicElastic, "isotropic elastic");
	REGISTER_FECORE_CLASS(FECoupledMooneyRivlin, "coupled Mooney-Rivlin");
	REGISTER_FECORE_CLASS(FECoupledVerondaWestmann, "coupled Veronda-Westmann");
	REGISTER_FECORE_CLASS(FENaturalNeoHookean, "natural neo-Hookean");
	REGISTER_FECORE_CLASS(FENeoHookean, "neo-Hookean");
	REGISTER_FECORE_CLASS(FENeoHookeanTransIso, "neo-Hookean transiso");
    REGISTER_FECORE_CLASS(FETraceFreeNeoHookean, "trace-free neo-Hookean");
	REGISTER_FECORE_CLASS(FENewtonianViscousSolid, "Newtonian viscous solid");
	REGISTER_FECORE_CLASS(FEOgdenUnconstrained, "Ogden unconstrained");
	REGISTER_FECORE_CLASS(FEOrthoElastic, "orthotropic elastic");
	REGISTER_FECORE_CLASS(FEOrthotropicCLE, "orthotropic CLE");
	REGISTER_FECORE_CLASS(FEOsmoticVirialExpansion, "osmotic virial expansion");
	REGISTER_FECORE_CLASS(FEPerfectOsmometer, "perfect osmometer");
	REGISTER_FECORE_CLASS(FESphericalFiberDistribution, "spherical fiber distribution");
	REGISTER_FECORE_CLASS(FEStVenantKirchhoff, "St.Venant-Kirchhoff");
    REGISTER_FECORE_CLASS(FEViscoElasticDamage, "viscoelastic damage");
	REGISTER_FECORE_CLASS(FEViscoElasticMaterial, "viscoelastic");
	REGISTER_FECORE_CLASS(FEElasticMultigeneration, "multigeneration");
	REGISTER_FECORE_CLASS(FERemodelingElasticMaterial, "remodeling solid");
	REGISTER_FECORE_CLASS(FECarterHayesOld, "Carter-Hayes (old)");
	REGISTER_FECORE_CLASS(FEContinuousFiberDistribution, "continuous fiber distribution");
	REGISTER_FECORE_CLASS(FECoupledTransIsoVerondaWestmann, "coupled trans-iso Veronda-Westmann");
	REGISTER_FECORE_CLASS(FECoupledTransIsoMooneyRivlin, "coupled trans-iso Mooney-Rivlin");
	REGISTER_FECORE_CLASS(FEGenericHyperelastic, "hyperelastic");
	REGISTER_FECORE_CLASS(FEGenericTransIsoHyperelastic, "trans-iso hyperelastic");
	REGISTER_FECORE_CLASS(FEDamageFiberPower, "damage fiber power");
	REGISTER_FECORE_CLASS(FEDamageFiberExponential, "damage fiber exponential");
	REGISTER_FECORE_CLASS(FEDamageFiberExpLinear, "damage fiber exp-linear");
	REGISTER_FECORE_CLASS(FEGenerationMaterial, "generation");
	REGISTER_FECORE_CLASS(FEHGOCoronary, "HGO-coronary");
    REGISTER_FECORE_CLASS(FELungMaterial, "lung");
    REGISTER_FECORE_CLASS(FEKinematicGrowth, "kinematic growth");

	// These materials are derived from FEElasticMaterial and use FEElasticMaterials
	REGISTER_FECORE_CLASS(FEElasticMixture, "solid mixture");
	REGISTER_FECORE_CLASS(FEReactiveViscoelasticMaterial, "reactive viscoelastic");
	REGISTER_FECORE_CLASS(FEDamageMaterial, "elastic damage");
	REGISTER_FECORE_CLASS(FERVEDamageMaterial, "reactive viscoelastic damage");
	REGISTER_FECORE_CLASS(FEReactiveFatigue, "reactive fatigue", FECORE_EXPERIMENTAL);
    REGISTER_FECORE_CLASS(FERVEFatigueMaterial, "reactive viscoelastic fatigue", FECORE_EXPERIMENTAL);
	REGISTER_FECORE_CLASS(FEReactivePlasticity, "reactive plasticity");
	REGISTER_FECORE_CLASS(FEReactivePlasticDamage, "reactive plastic damage");

	// Uncoupled elastic materials (derived from FEUncoupledMaterial)
	REGISTER_FECORE_CLASS(FEArrudaBoyce, "Arruda-Boyce");
	REGISTER_FECORE_CLASS(FE2DTransIsoMooneyRivlin, "2D trans iso Mooney-Rivlin");
	REGISTER_FECORE_CLASS(FE2DTransIsoVerondaWestmann, "2D trans iso Veronda-Westmann");
	REGISTER_FECORE_CLASS(FEDamageMooneyRivlin, "damage Mooney-Rivlin");
	REGISTER_FECORE_CLASS(FEDamageTransIsoMooneyRivlin, "damage trans iso Mooney-Rivlin");
	REGISTER_FECORE_CLASS(FEEFDMooneyRivlin, "EFD Mooney-Rivlin");
	REGISTER_FECORE_CLASS(FEEFDUncoupled, "EFD uncoupled");
	REGISTER_FECORE_CLASS(FEEFDVerondaWestmann, "EFD Veronda-Westmann");
	REGISTER_FECORE_CLASS(FEFungOrthotropic, "Fung orthotropic");
	REGISTER_FECORE_CLASS(FEHolzapfelGasserOgden, "Holzapfel-Gasser-Ogden");
	REGISTER_FECORE_CLASS(FEGentMaterial, "Gent");
	REGISTER_FECORE_CLASS(FEIncompNeoHookean, "incomp neo-Hookean");
	REGISTER_FECORE_CLASS(FEMooneyRivlin, "Mooney-Rivlin");
	REGISTER_FECORE_CLASS(FEMuscleMaterial, "muscle material");
	REGISTER_FECORE_CLASS(FENewtonianViscousSolidUC, "Newtonian viscous solid uncoupled");
	REGISTER_FECORE_CLASS(FEOgdenMaterial, "Ogden");
	REGISTER_FECORE_CLASS(FETCNonlinearOrthotropic, "TC nonlinear orthotropic");
	REGISTER_FECORE_CLASS(FETendonMaterial, "tendon material");
	REGISTER_FECORE_CLASS(FETransIsoMooneyRivlin, "trans iso Mooney-Rivlin");
    REGISTER_FECORE_CLASS(FETransIsoMREstrada, "trans iso MR-Estrada");
	REGISTER_FECORE_CLASS(FETransIsoVerondaWestmann, "trans iso Veronda-Westmann");
	REGISTER_FECORE_CLASS(FEUncoupledElasticMixture, "uncoupled solid mixture");
	REGISTER_FECORE_CLASS(FEVerondaWestmann, "Veronda-Westmann");
    REGISTER_FECORE_CLASS(FEUncoupledViscoElasticDamage, "uncoupled viscoelastic damage");
	REGISTER_FECORE_CLASS(FEUncoupledViscoElasticMaterial, "uncoupled viscoelastic");
	REGISTER_FECORE_CLASS(FEMRVonMisesFibers, "Mooney-Rivlin von Mises Fibers");
	REGISTER_FECORE_CLASS(FEUncoupledActiveContraction, "uncoupled active contraction");
	REGISTER_FECORE_CLASS(FEContinuousFiberDistributionUC, "continuous fiber distribution uncoupled");
	REGISTER_FECORE_CLASS(FEPRLig, "PRLig");
	REGISTER_FECORE_CLASS(FEUncoupledReactiveViscoelasticMaterial, "uncoupled reactive viscoelastic");
	REGISTER_FECORE_CLASS(FEDamageMaterialUC, "uncoupled elastic damage");
    REGISTER_FECORE_CLASS(FEUncoupledReactiveFatigue, "uncoupled reactive fatigue", FECORE_EXPERIMENTAL);
	REGISTER_FECORE_CLASS(FEGenericHyperelasticUC, "uncoupled hyperelastic");
	REGISTER_FECORE_CLASS(FEGenericTransIsoHyperelasticUC, "uncoupled trans-iso hyperelastic");
	REGISTER_FECORE_CLASS(FEIsotropicLeeSacks, "isotropic Lee-Sacks");
	REGISTER_FECORE_CLASS(FEIsotropicLeeSacksUncoupled, "uncoupled isotropic Lee-Sacks");
	REGISTER_FECORE_CLASS(FEPolynomialHyperElastic, "polynomial");
	REGISTER_FECORE_CLASS(FEShenoyMaterial, "Shenoy");
	REGISTER_FECORE_CLASS(FEFiberEFDNeoHookean, "fiber neo-Hookean");
    REGISTER_FECORE_CLASS(FEYeoh, "Yeoh");

	// fiber materials (derived from FEFiberMaterial)
    REGISTER_FECORE_CLASS(FEFiberCDF         , "fiber-CDF"           );
	REGISTER_FECORE_CLASS(FEFiberNH          , "fiber-NH"            );
	REGISTER_FECORE_CLASS(FEFiberExpPow      , "fiber-exp-pow"       );
	REGISTER_FECORE_CLASS(FEFiberExpLinear   , "fiber-exp-linear"    );
	REGISTER_FECORE_CLASS(FEFiberPowLinear   , "fiber-pow-linear"    );
	REGISTER_FECORE_CLASS(FEFiberExpPowLinear, "fiber-exp-pow-linear");
	REGISTER_FECORE_CLASS(FEFiberNaturalNH   , "fiber-natural-NH"    );
    REGISTER_FECORE_CLASS(FEFiberEntropyChain, "fiber-entropy-chain" );
    REGISTER_FECORE_CLASS(FEVolumeGrowth     , "volume growth"       );
    REGISTER_FECORE_CLASS(FEAreaGrowth       , "area growth"         );
    REGISTER_FECORE_CLASS(FEFiberGrowth      , "fiber growth"        );

	// Elastic Fiber materials (derived from FEElasticFiberMaterial)
    REGISTER_FECORE_CLASS(FEElasticFiberCDF         , "fiber-CDF"           );
	REGISTER_FECORE_CLASS(FEElasticFiberNH          , "fiber-NH"            );
	REGISTER_FECORE_CLASS(FEElasticFiberExpPow      , "fiber-exp-pow"       );
	REGISTER_FECORE_CLASS(FEElasticFiberExpLinear   , "fiber-exp-linear"    );
	REGISTER_FECORE_CLASS(FEElasticFiberPowLinear   , "fiber-pow-linear"    );
	REGISTER_FECORE_CLASS(FEElasticFiberExpPowLinear, "fiber-exp-pow-linear");
	REGISTER_FECORE_CLASS(FEElasticFiberNaturalNH   , "fiber-natural-NH"    );
    REGISTER_FECORE_CLASS(FEElasticFiberEntropyChain, "fiber-entropy-chain" );

	// fiber materials for uncoupled formulation (derived from FEFiberMaterialUC)
    REGISTER_FECORE_CLASS(FEFiberCDFUncoupled  , "fiber-CDF-uncoupled"       );
	REGISTER_FECORE_CLASS(FEFiberExpLinearUC   , "uncoupled fiber-exp-linear");
	REGISTER_FECORE_CLASS(FEFiberNHUC          , "fiber-NH-uncoupled");
	REGISTER_FECORE_CLASS(FEFiberExpPowUC      , "fiber-exp-pow-uncoupled");
	REGISTER_FECORE_CLASS(FEFiberPowLinearUC   , "fiber-pow-linear-uncoupled");
    REGISTER_FECORE_CLASS(FEFiberEntropyChainUC, "uncoupled fiber-entropy-chain");

	// Uncoupled elastic fiber materials (derived from FEUncoupledFiberMaterial)
    REGISTER_FECORE_CLASS(FEElasticFiberCDFUncoupled    , "fiber-CDF-uncoupled"       );
	REGISTER_FECORE_CLASS(FEUncoupledFiberExpLinear     , "uncoupled fiber-exp-linear");
	REGISTER_FECORE_CLASS(FEUncoupledFiberNH            , "fiber-NH-uncoupled");
	REGISTER_FECORE_CLASS(FEUncoupledFiberExpPow        , "fiber-exp-pow-uncoupled");
	REGISTER_FECORE_CLASS(FEUncoupledFiberPowLinear     , "fiber-pow-linear-uncoupled");
    REGISTER_FECORE_CLASS(FEUncoupledFiberKiousis       , "fiber-Kiousis-uncoupled");
    REGISTER_FECORE_CLASS(FEUncoupledFiberEntropyChainUC, "uncoupled fiber-entropy-chain");

	// obsolete fiber materials
	REGISTER_FECORE_CLASS(FEFiberExponentialPower, "fiber-exponential-power-law");
	REGISTER_FECORE_CLASS(FEFiberExponentialPowerUC, "fiber-exponential-power-law-uncoupled");

	// solid materials (derived from FESolidMaterial)
	REGISTER_FECORE_CLASS(FERigidMaterial, "rigid body");
	REGISTER_FECORE_CLASS(FEVonMisesPlasticity, "von-Mises plasticity");

	// Fiber density distributions for CFD materials
	REGISTER_FECORE_CLASS(FESphericalFiberDensityDistribution, "spherical");
	REGISTER_FECORE_CLASS(FEEllipsoidalFiberDensityDistribution, "ellipsoidal");
	REGISTER_FECORE_CLASS(FEVonMises3DFiberDensityDistribution, "von-Mises-3d");
	REGISTER_FECORE_CLASS(FEVonMises3DTwoFDDAxisymmetric, "von-Mises-3d-two-axisym");
	REGISTER_FECORE_CLASS(FECircularFiberDensityDistribution, "circular");
	REGISTER_FECORE_CLASS(FEEllipticalFiberDensityDistribution, "elliptical");
	REGISTER_FECORE_CLASS(FEVonMises2DFiberDensityDistribution, "von-Mises-2d");
	REGISTER_FECORE_CLASS(FEStructureTensorDistribution, "structure-tensor");

	// Fiber distribution integration schemes for CFD materials
	REGISTER_FECORE_CLASS(FEFiberIntegrationGauss, "fibers-3d-gauss");
	REGISTER_FECORE_CLASS(FEFiberIntegrationGeodesic, "fibers-3d-geodesic");
	REGISTER_FECORE_CLASS(FEFiberIntegrationGaussKronrod, "fibers-3d-gkt");
	REGISTER_FECORE_CLASS(FEFiberIntegrationTriangle, "fibers-3d-fei");
	REGISTER_FECORE_CLASS(FEFiberIntegrationTrapezoidal, "fibers-2d-trapezoidal");

	// Other materials 
	REGISTER_FECORE_CLASS(FELinearTrussMaterial, "linear truss");
	REGISTER_FECORE_CLASS(FEHuiskesSupply, "Huiskes-supply");
	REGISTER_FECORE_CLASS(FEActiveFiberContraction, "active_contraction");
    REGISTER_FECORE_CLASS(FEForceVelocityContraction, "force-velocity-Estrada");
	REGISTER_FECORE_CLASS(FEWrinkleOgdenMaterial, "wrinkle Ogden");
	REGISTER_FECORE_CLASS(FEElasticMembrane, "elastic membrane");

	// active contraction materials
	REGISTER_FECORE_CLASS(FEPrescribedActiveContractionUniaxial, "prescribed uniaxial active contraction");
	REGISTER_FECORE_CLASS(FEPrescribedActiveContractionUniaxialUC, "uncoupled prescribed uniaxial active contraction");
	REGISTER_FECORE_CLASS(FEPrescribedActiveContractionTransIso, "prescribed trans iso active contraction");
	REGISTER_FECORE_CLASS(FEPrescribedActiveContractionTransIsoUC, "uncoupled prescribed trans iso active contraction");
	REGISTER_FECORE_CLASS(FEPrescribedActiveContractionIsotropic, "prescribed isotropic active contraction");
	REGISTER_FECORE_CLASS(FEPrescribedActiveContractionIsotropicUC, "uncoupled prescribed isotropic active contraction");
	REGISTER_FECORE_CLASS(FEPrescribedActiveContractionFiber, "prescribed fiber active contraction");
	REGISTER_FECORE_CLASS(FEPrescribedActiveContractionFiberUC, "uncoupled prescribed fiber active contraction");

	REGISTER_FECORE_CLASS(FEActiveFiberStress, "active fiber stress");
	REGISTER_FECORE_CLASS(FEActiveFiberStressUC, "uncoupled active fiber stress");

	// discrete materials
	REGISTER_FECORE_CLASS(FECompositeDiscreteMaterial, "discrete composite");
	REGISTER_FECORE_CLASS(FELinearSpring, "linear spring");
	REGISTER_FECORE_CLASS(FETensionOnlyLinearSpring, "tension-only linear spring");
	REGISTER_FECORE_CLASS(FENonlinearSpringMaterial, "nonlinear spring");
	REGISTER_FECORE_CLASS(FEExperimentalSpring, "experimental spring");
	REGISTER_FECORE_CLASS(FEDiscreteContractileMaterial, "Hill");
	REGISTER_FECORE_CLASS(FETorsionalSpring, "torsion spring");

	// bond relaxation materials (used by reactive visco-elastic materials)
	REGISTER_FECORE_CLASS(FEBondRelaxationExponential, "relaxation-exponential");
	REGISTER_FECORE_CLASS(FEBondRelaxationExpDistortion, "relaxation-exp-distortion");
    REGISTER_FECORE_CLASS(FEBondRelaxationExpDistUser, "relaxation-exp-dist-user");
	REGISTER_FECORE_CLASS(FEBondRelaxationFung, "relaxation-Fung");
	REGISTER_FECORE_CLASS(FEBondRelaxationPark, "relaxation-Park");
	REGISTER_FECORE_CLASS(FEBondRelaxationParkDistortion, "relaxation-Park-distortion");
    REGISTER_FECORE_CLASS(FEBondRelaxationParkDistUser, "relaxation-Park-dist-user");
	REGISTER_FECORE_CLASS(FEBondRelaxationPower, "relaxation-power");
	REGISTER_FECORE_CLASS(FEBondRelaxationPowerDistortion, "relaxation-power-distortion");
    REGISTER_FECORE_CLASS(FEBondRelaxationPowerDistUser, "relaxation-power-dist-user");
	REGISTER_FECORE_CLASS(FEBondRelaxationCarreau, "relaxation-Carreau");
    REGISTER_FECORE_CLASS(FEBondRelaxationProny, "relaxation-Prony");
    REGISTER_FECORE_CLASS(FEBondRelaxationMalkin, "relaxation-Malkin");
    REGISTER_FECORE_CLASS(FEBondRelaxationMalkinDist, "relaxation-Malkin-distortion");
    REGISTER_FECORE_CLASS(FEBondRelaxationMalkinDistUser, "relaxation-Malkin-dist-user");
    REGISTER_FECORE_CLASS(FEBondRelaxationCSexp, "relaxation-CSexp");
    REGISTER_FECORE_CLASS(FEBondRelaxationCSexpDistUser, "relaxation-CSexp-dist-user");

    // bond recruitment materials (used by reactive visco-elastic materials)
    REGISTER_FECORE_CLASS(FEBondRecruitmentUser, "recruitment user");
    REGISTER_FECORE_CLASS(FEBondRecruitmentPower, "recruitment power");
    REGISTER_FECORE_CLASS(FEBondRecruitmentExp, "recruitment exponential");
    REGISTER_FECORE_CLASS(FEBondRecruitmentPoly, "recruitment polynomial");
    REGISTER_FECORE_CLASS(FEBondRecruitmentLogNormal, "recruitment log-normal");
    REGISTER_FECORE_CLASS(FEBondRecruitmentWeibull, "recruitment Weibull");
    REGISTER_FECORE_CLASS(FEBondRecruitmentPQP, "recruitment quintic");
    REGISTER_FECORE_CLASS(FEBondRecruitmentGamma, "recruitment gamma");

	// damage cumulative distribution functions (used by damage materials)
	REGISTER_FECORE_CLASS(FEDamageCDFSimo, "CDF Simo");
	REGISTER_FECORE_CLASS(FEDamageCDFLogNormal, "CDF log-normal");
	REGISTER_FECORE_CLASS(FEDamageCDFWeibull, "CDF Weibull");
	REGISTER_FECORE_CLASS(FEDamageCDFStep, "CDF step");
	REGISTER_FECORE_CLASS(FEDamageCDFPQP, "CDF quintic");
	REGISTER_FECORE_CLASS(FEDamageCDFGamma, "CDF gamma");
	REGISTER_FECORE_CLASS(FEDamageCDFUser, "CDF user");

	// damage criterion (used by damage and plastic materials)
	REGISTER_FECORE_CLASS(FEDamageCriterionSimo, "DC Simo");
	REGISTER_FECORE_CLASS(FEDamageCriterionSED, "DC strain energy density");
	REGISTER_FECORE_CLASS(FEDamageCriterionSSE, "DC specific strain energy");
	REGISTER_FECORE_CLASS(FEDamageCriterionVMS, "DC von Mises stress");
    REGISTER_FECORE_CLASS(FEDamageCriterionDrucker, "DC Drucker shear stress");
	REGISTER_FECORE_CLASS(FEDamageCriterionMSS, "DC max shear stress");
	REGISTER_FECORE_CLASS(FEDamageCriterionMNS, "DC max normal stress");
	REGISTER_FECORE_CLASS(FEDamageCriterionMNLS, "DC max normal Lagrange strain");
	REGISTER_FECORE_CLASS(FEDamageCriterionOSS, "DC octahedral shear strain");
    REGISTER_FECORE_CLASS(FEDamageCriterionONS, "DC octahedral natural strain");

    // plastic flow curve (used by plastic materials)
    REGISTER_FECORE_CLASS(FEPlasticFlowCurvePaper, "PFC paper");
    REGISTER_FECORE_CLASS(FEPlasticFlowCurveUser , "PFC user");
    REGISTER_FECORE_CLASS(FEPlasticFlowCurveMath , "PFC math");

	// prestrain materials
	REGISTER_FECORE_CLASS(FEPrestrainElastic, "prestrain elastic");
	REGISTER_FECORE_CLASS(FEPreStrainUncoupledElastic, "uncoupled prestrain elastic");
	REGISTER_FECORE_CLASS(FEConstPrestrainGradient, "prestrain gradient");
	REGISTER_FECORE_CLASS(FEInSituStretchGradient, "in-situ stretch");

	// beam materials
	REGISTER_FECORE_CLASS(FEElasticBeamMaterial, "linear-beam");

	//-----------------------------------------------------------------------------
	// domain classes
	REGISTER_FECORE_CLASS(FERigidSolidDomain, "rigid-solid");
	REGISTER_FECORE_CLASS(FERigidShellDomain, "rigid-shell");
	REGISTER_FECORE_CLASS(FERigidShellDomainOld, "rigid-shell-old");
	REGISTER_FECORE_CLASS(FERemodelingElasticDomain, "remodeling-solid");
	REGISTER_FECORE_CLASS(FE3FieldElasticSolidDomain, "three-field-solid");
	REGISTER_FECORE_CLASS(FE3FieldElasticShellDomain, "three-field-shell");
	REGISTER_FECORE_CLASS(FEUDGHexDomain, "udg-hex");
	REGISTER_FECORE_CLASS(FESRIElasticSolidDomain, "sri-solid");
	REGISTER_FECORE_CLASS(FEUT4Domain, "ut4-solid");
	REGISTER_FECORE_CLASS(FEStandardElasticSolidDomain, "elastic-solid");
	REGISTER_FECORE_CLASS(FEElasticShellDomain, "elastic-shell");
	REGISTER_FECORE_CLASS(FEElasticShellDomainOld, "elastic-shell-old");
	REGISTER_FECORE_CLASS(FEElasticEASShellDomain, "elastic-shell-eas");
	REGISTER_FECORE_CLASS(FEElasticANSShellDomain, "elastic-shell-ans");
	REGISTER_FECORE_CLASS(FELinearTrussDomain, "linear-truss");
	REGISTER_FECORE_CLASS(FEElasticTrussDomain, "elastic-truss");
	REGISTER_FECORE_CLASS(FEElasticBeamDomain, "linear-beam");
	REGISTER_FECORE_CLASS(FEDiscreteElasticDomain, "discrete");
	REGISTER_FECORE_CLASS(FEDeformableSpringDomain, "deformable-spring");
	REGISTER_FECORE_CLASS(FEDeformableSpringDomain2, "deformable-spring2");

	//-----------------------------------------------------------------------------
	// classes derived from FEBoundaryCondition
	REGISTER_FECORE_CLASS(FEFixedDisplacement           , "zero displacement");
	REGISTER_FECORE_CLASS(FEFixedRotation               , "zero rotation");
	REGISTER_FECORE_CLASS(FEFixedShellDisplacement      , "zero shell displacement");
	REGISTER_FECORE_CLASS(FEPrescribedDisplacement      , "prescribed displacement");
	REGISTER_FECORE_CLASS(FEPrescribedRotation          , "prescribed rotation");
	REGISTER_FECORE_CLASS(FEPrescribedShellDisplacement , "prescribed shell displacement");
	REGISTER_FECORE_CLASS(FEBCPrescribedDeformation     , "prescribed deformation");
	REGISTER_FECORE_CLASS(FEPrescribedNormalDisplacement, "normal displacement");
	REGISTER_FECORE_CLASS(FEBCRigidDeformation          , "rigid deformation");
	REGISTER_FECORE_CLASS(FEBCPrescribedDeformation2O   , "prescribed deformation 2O");
	REGISTER_FECORE_CLASS(FERigidNodeSet                , "rigid");

	//-----------------------------------------------------------------------------
	// classes derived from FEInitialCondition
	REGISTER_FECORE_CLASS(FEInitialVelocity, "velocity");
	REGISTER_FECORE_CLASS(FEInitialShellVelocity, "shell velocity");
	REGISTER_FECORE_CLASS(FEInitialPreStrain, "prestrain");

	//-----------------------------------------------------------------------------
	// classes derived from FENodalLoad
	REGISTER_FECORE_CLASS(FENodalForce, "nodal_force");

	//-----------------------------------------------------------------------------
	// classes derived from FESurfaceLoad
	REGISTER_FECORE_CLASS(FEPressureLoad, "pressure");
	REGISTER_FECORE_CLASS(FETractionLoad, "traction");
    REGISTER_FECORE_CLASS(FESurfaceForceUniform, "force");
    REGISTER_FECORE_CLASS(FEBearingLoad, "bearing load");

	//-----------------------------------------------------------------------------
	// classes derived from FEBodyForce
	REGISTER_FECORE_CLASS(FEConstBodyForceOld, "const");	// obsolete in 3.0
	REGISTER_FECORE_CLASS(FENonConstBodyForceOld, "non-const");	// obsolete in 3.0

	REGISTER_FECORE_CLASS(FEGenericBodyForce, "body force");
	REGISTER_FECORE_CLASS(FECentrifugalBodyForce, "centrifugal");
	REGISTER_FECORE_CLASS(FEPointBodyForce, "point", FECORE_EXPERIMENTAL);
	REGISTER_FECORE_CLASS(FESurfaceAttractionBodyForce, "surface attraction");
	REGISTER_FECORE_CLASS(FEMassDamping, "mass damping");

	//-----------------------------------------------------------------------------
	// constraint classes
	REGISTER_FECORE_CLASS(FEPointConstraint, "point");
	REGISTER_FECORE_CLASS(FESymmetryPlane, "symmetry plane");
	REGISTER_FECORE_CLASS(FERigidJoint, "rigid joint");
	REGISTER_FECORE_CLASS(FEGenericRigidJoint, "generic rigid joint");
	REGISTER_FECORE_CLASS(FERigidSphericalJoint, "rigid spherical joint");
	REGISTER_FECORE_CLASS(FERigidRevoluteJoint, "rigid revolute joint");
	REGISTER_FECORE_CLASS(FERigidPrismaticJoint, "rigid prismatic joint");
	REGISTER_FECORE_CLASS(FERigidCylindricalJoint, "rigid cylindrical joint");
	REGISTER_FECORE_CLASS(FERigidPlanarJoint, "rigid planar joint");
	REGISTER_FECORE_CLASS(FERigidLock, "rigid lock");
	REGISTER_FECORE_CLASS(FERigidSpring, "rigid spring");
	REGISTER_FECORE_CLASS(FERigidDamper, "rigid damper");
	REGISTER_FECORE_CLASS(FERigidAngularDamper, "rigid angular damper");
	REGISTER_FECORE_CLASS(FERigidContractileForce, "rigid contractile force");
	REGISTER_FECORE_CLASS(FEVolumeConstraint, "volume");
	REGISTER_FECORE_CLASS(FEDiscreteContact, "discrete contact");
	REGISTER_FECORE_CLASS(FEDiscreteContact2, "discrete contact2");
	REGISTER_FECORE_CLASS(FEDistanceConstraint, "node distance");
	REGISTER_FECORE_CLASS(FEGPAConstraint, "prestrain");
	REGISTER_FECORE_CLASS(FEInSituStretchConstraint, "in-situ stretch");
	REGISTER_FECORE_CLASS(FEAzimuthConstraint, "azimuth constraint");
    REGISTER_FECORE_CLASS(FEFixedNormalDisplacement, "fixed normal displacement");

	// Lagrange multiplier constraints
	REGISTER_FECORE_CLASS(FENodeToNodeConstraint, "node-on-node");

	//-----------------------------------------------------------------------------
	// classes derived from FEContactInterface

	REGISTER_FECORE_CLASS(FESlidingInterface, "sliding-node-on-facet");
	REGISTER_FECORE_CLASS(FEFacet2FacetSliding, "sliding-facet-on-facet");
	REGISTER_FECORE_CLASS(FESlidingElasticInterface, "sliding-elastic");
	REGISTER_FECORE_CLASS(FETiedInterface, "tied-node-on-facet");
	REGISTER_FECORE_CLASS(FEFacet2FacetTied, "tied-facet-on-facet");
	REGISTER_FECORE_CLASS(FETiedElasticInterface, "tied-elastic");
	REGISTER_FECORE_CLASS(FEPeriodicBoundary, "periodic boundary");
	REGISTER_FECORE_CLASS(FERigidWallInterface, "rigid_wall");
	REGISTER_FECORE_CLASS(FEPeriodicSurfaceConstraint, "surface constraint");
	REGISTER_FECORE_CLASS(FEStickyInterface, "sticky");
	REGISTER_FECORE_CLASS(FEMortarSlidingContact, "mortar-sliding", FECORE_EXPERIMENTAL);
	REGISTER_FECORE_CLASS(FEMortarTiedContact, "mortar-tied", FECORE_EXPERIMENTAL);
	REGISTER_FECORE_CLASS(FEContactPotential, "contact potential");

	//-----------------------------------------------------------------------------
	// classes derived directly from FERigidBC
	REGISTER_FECORE_CLASS(FERigidFixedBCNew     , "rigid_fixed"           );
	REGISTER_FECORE_CLASS(FERigidDisplacement   , "rigid_displacement"    );
	REGISTER_FECORE_CLASS(FERigidRotation       , "rigid_rotation"        );

	REGISTER_FECORE_CLASS(FERigidFixedBCOld     , "rigid_fixed_old"     , 0x300);	// obsolete in 4.0
	REGISTER_FECORE_CLASS(FERigidPrescribedOld  , "rigid_prescribed_old", 0x300);	// obsolete in 4.0
	
	// classes derived directly from FERigidIC
	REGISTER_FECORE_CLASS(FERigidBodyVelocity       , "initial_rigid_velocity"        );
	REGISTER_FECORE_CLASS(FERigidBodyAngularVelocity, "initial_rigid_angular_velocity");

	//-----------------------------------------------------------------------------
	// classes derived directly from FERigidLoad
	REGISTER_FECORE_CLASS(FERigidAxialForce    , "rigid_axial_force"    );
	REGISTER_FECORE_CLASS(FERigidBodyForce     , "rigid_force"          );
	REGISTER_FECORE_CLASS(FERigidBodyMoment    , "rigid_moment"         );
    REGISTER_FECORE_CLASS(FERigidFollowerForce , "rigid_follower_force" );
    REGISTER_FECORE_CLASS(FERigidFollowerMoment, "rigid_follower_moment");
	REGISTER_FECORE_CLASS(FERigidCable         , "rigid_cable");
	REGISTER_FECORE_CLASS(FERigidCablePoint	   , "rigid_cable_point");

	//-----------------------------------------------------------------------------
	// classes derived from FEPlotData
	REGISTER_FECORE_CLASS(FEPlotElementVelocity, "velocity");
	REGISTER_FECORE_CLASS(FEPlotElementAcceleration, "acceleration");
	REGISTER_FECORE_CLASS(FEPlotDensity, "density");
	REGISTER_FECORE_CLASS(FEPlotElementStress, "stress");
	REGISTER_FECORE_CLASS(FEPlotElementPK2Stress, "PK2 stress");
	REGISTER_FECORE_CLASS(FEPlotElementPK1Stress, "PK1 stress");
	REGISTER_FECORE_CLASS(FEPlotElementMixtureStress, "mixture stress");
	REGISTER_FECORE_CLASS(FEPlotElementUncoupledPressure, "uncoupled pressure");
	REGISTER_FECORE_CLASS(FEPlotElementElasticity, "elasticity");
    REGISTER_FECORE_CLASS(FEPlotElementDevElasticity, "deviatoric elasticity");
	REGISTER_FECORE_CLASS(FEPlotRelativeVolume, "relative volume");
	REGISTER_FECORE_CLASS(FEPlotShellRelativeVolume, "shell relative volume");// , FECORE_SPEC(3, 0)); // NOTE: deprecated
	REGISTER_FECORE_CLASS(FEPlotFiberVector, "fiber vector");
	REGISTER_FECORE_CLASS(FEPlotFiberStretch, "fiber stretch");
	REGISTER_FECORE_CLASS(FEPlotDevFiberStretch, "deviatoric fiber stretch");
	REGISTER_FECORE_CLASS(FEPlotMaterialAxes, "material axes");
	REGISTER_FECORE_CLASS(FEPlotShellThickness, "shell thickness");
	REGISTER_FECORE_CLASS(FEPlotShellDirector, "shell director");
	REGISTER_FECORE_CLASS(FEPlotDamage, "damage");
	REGISTER_FECORE_CLASS(FEPlotIntactBondFraction, "intact bond fraction");
    REGISTER_FECORE_CLASS(FEPlotFatigueBondFraction, "fatigue bond fraction");
	REGISTER_FECORE_CLASS(FEPlotYieldedBondFraction, "yielded bond fraction");
	REGISTER_FECORE_CLASS(FEPlotOctahedralPlasticStrain, "octahedral plastic strain");
	REGISTER_FECORE_CLASS(FEPlotReactivePlasticityHeatSupply, "plasticity heat supply density");
	REGISTER_FECORE_CLASS(FEPlotMixtureVolumeFraction, "volume fraction");
	REGISTER_FECORE_CLASS(FEPlotUT4NodalStresses, "ut4 nodal stress");
	REGISTER_FECORE_CLASS(FEPlotContactGap, "contact gap");
	REGISTER_FECORE_CLASS(FEPlotNodalContactGap, "nodal contact gap");
	REGISTER_FECORE_CLASS(FEPlotVectorGap, "vector gap");
	REGISTER_FECORE_CLASS(FEPlotNodalVectorGap, "nodal vector gap");
	REGISTER_FECORE_CLASS(FEPlotContactPressure, "contact pressure");
	REGISTER_FECORE_CLASS(FEPlotNodalContactPressure, "nodal contact pressure");
	REGISTER_FECORE_CLASS(FEPlotContactTraction, "contact traction");
	REGISTER_FECORE_CLASS(FEPlotNodalContactTraction, "nodal contact traction");
	REGISTER_FECORE_CLASS(FEPlotStickStatus, "contact stick");
	REGISTER_FECORE_CLASS(FEPlotContactForce, "contact force");
	REGISTER_FECORE_CLASS(FEPlotContactArea, "contact area");
	REGISTER_FECORE_CLASS(FEPlotContactPenalty, "contact penalty");
	REGISTER_FECORE_CLASS(FEPlotContactStatus, "contact status");
	REGISTER_FECORE_CLASS(FEPlotSPRStresses, "SPR stress");
	REGISTER_FECORE_CLASS(FEPlotSPRLinearStresses, "SPR-P1 stress");
	REGISTER_FECORE_CLASS(FEPlotSPRPrincStresses, "SPR principal stress");
	REGISTER_FECORE_CLASS(FEPlotNodalStresses, "nodal stress");
	REGISTER_FECORE_CLASS(FEPlotShellStrain, "shell strain");
	REGISTER_FECORE_CLASS(FEPlotDeformationGradient, "deformation gradient");
	REGISTER_FECORE_CLASS(FEPlotLagrangeStrain, "Lagrange strain");
	REGISTER_FECORE_CLASS(FEPlotInfStrain, "infinitesimal strain");
	REGISTER_FECORE_CLASS(FEPlotSPRLagrangeStrain, "SPR Lagrange strain");
    REGISTER_FECORE_CLASS(FEPlotRightStretch, "right stretch");
    REGISTER_FECORE_CLASS(FEPlotLeftStretch, "left stretch");
    REGISTER_FECORE_CLASS(FEPlotRightHencky, "right Hencky");
    REGISTER_FECORE_CLASS(FEPlotLeftHencky, "left Hencky");
    REGISTER_FECORE_CLASS(FEPlotRateOfDeformation, "rate of deformation");
	REGISTER_FECORE_CLASS(FEPlotMortarContactGap, "mortar-gap");
	REGISTER_FECORE_CLASS(FEPlotSurfaceTraction, "surface traction");
	REGISTER_FECORE_CLASS(FEPlotNodalSurfaceTraction, "nodal surface traction");
	REGISTER_FECORE_CLASS(FEPlotEnclosedVolume, "enclosed volume");
	REGISTER_FECORE_CLASS(FEPlotSurfaceArea, "surface area");
	REGISTER_FECORE_CLASS(FEPlotFacetArea, "facet area");
	REGISTER_FECORE_CLASS(FEPlotStrainEnergyDensity, "strain energy density");
	REGISTER_FECORE_CLASS(FEPlotDevStrainEnergyDensity, "deviatoric strain energy density");
	REGISTER_FECORE_CLASS(FEPlotSpecificStrainEnergy, "specific strain energy");
	REGISTER_FECORE_CLASS(FEPlotKineticEnergyDensity, "kinetic energy density");
	REGISTER_FECORE_CLASS(FEPlotElementStrainEnergy, "element strain energy");
	REGISTER_FECORE_CLASS(FEPlotElementKineticEnergy, "element kinetic energy");
	REGISTER_FECORE_CLASS(FEPlotElementCenterOfMass, "element center of mass");
	REGISTER_FECORE_CLASS(FEPlotElementLinearMomentum, "element linear momentum");
	REGISTER_FECORE_CLASS(FEPlotElementAngularMomentum, "element angular momentum");
	REGISTER_FECORE_CLASS(FEPlotElementStressPower, "element stress power");
	REGISTER_FECORE_CLASS(FEPlotCurrentElementStrainEnergy, "current element strain energy");
	REGISTER_FECORE_CLASS(FEPlotCurrentElementKineticEnergy, "current element kinetic energy");
	REGISTER_FECORE_CLASS(FEPlotCurrentElementCenterOfMass, "current element center of mass");
	REGISTER_FECORE_CLASS(FEPlotCurrentElementLinearMomentum, "current element linear momentum");
	REGISTER_FECORE_CLASS(FEPlotCurrentElementAngularMomentum, "current element angular momentum");
	REGISTER_FECORE_CLASS(FEPlotNodeDisplacement, "displacement");
	REGISTER_FECORE_CLASS(FEPlotNodeRotation, "rotation");
	REGISTER_FECORE_CLASS(FEPlotNodeVelocity, "nodal velocity");
	REGISTER_FECORE_CLASS(FEPlotNodeAcceleration, "nodal acceleration");
	REGISTER_FECORE_CLASS(FEPlotNodeReactionForces, "reaction forces");
	REGISTER_FECORE_CLASS(FEPlotRigidReactionForce, "rigid force");
	REGISTER_FECORE_CLASS(FEPlotRigidReactionTorque, "rigid torque");
	REGISTER_FECORE_CLASS(FEPlotRigidDisplacement, "rigid position");
	REGISTER_FECORE_CLASS(FEPlotRigidVelocity, "rigid velocity");
	REGISTER_FECORE_CLASS(FEPlotRigidAcceleration, "rigid acceleration");
	REGISTER_FECORE_CLASS(FEPlotRigidRotation, "rigid angular position");
	REGISTER_FECORE_CLASS(FEPlotRigidAngularVelocity, "rigid angular velocity");
	REGISTER_FECORE_CLASS(FEPlotRigidAngularAcceleration, "rigid angular acceleration");
	REGISTER_FECORE_CLASS(FEPlotRigidLinearMomentum, "rigid linear momentum");
	REGISTER_FECORE_CLASS(FEPlotRigidAngularMomentum, "rigid angular momentum");
	REGISTER_FECORE_CLASS(FEPlotRigidKineticEnergy, "rigid kinetic energy");
	REGISTER_FECORE_CLASS(FEPlotRigidEuler, "Euler angle");
	REGISTER_FECORE_CLASS(FEPlotRigidRotationVector, "rigid rotation vector");
	REGISTER_FECORE_CLASS(FEPlotScalarSurfaceLoad, "scalar surface load");
	REGISTER_FECORE_CLASS(FEPlotNetSurfaceReactionForce, "surface reaction force");
	REGISTER_FECORE_CLASS(FEPlotNetSurfaceReactionMoment, "surface reaction moment");
	REGISTER_FECORE_CLASS(FEPlotStressError, "stress error");
	REGISTER_FECORE_CLASS(FEPlotFiberTargetStretch, "in-situ target stretch");
	REGISTER_FECORE_CLASS(FEPlotPreStrainStretch, "prestrain stretch");
	REGISTER_FECORE_CLASS(FEPlotPreStrainStretchError, "prestrain stretch error");
	REGISTER_FECORE_CLASS(FEPlotPreStrainCorrection, "prestrain correction");
	REGISTER_FECORE_CLASS(FEPlotSPRPreStrainCorrection, "SPR prestrain correction");
	REGISTER_FECORE_CLASS(FEPlotPreStrainCompatibility, "prestrain compatibility");
	REGISTER_FECORE_CLASS(FEPlotDiscreteElementStretch, "discrete element stretch");
	REGISTER_FECORE_CLASS(FEPlotDiscreteElementElongation, "discrete element elongation");
	REGISTER_FECORE_CLASS(FEPlotDiscreteElementPercentElongation, "discrete element percent elongation");
	REGISTER_FECORE_CLASS(FEPlotDiscreteElementForce, "discrete element force");
	REGISTER_FECORE_CLASS(FEPlotDiscreteElementSignedForce, "discrete element signed force");
	REGISTER_FECORE_CLASS(FEPlotDiscreteElementStrainEnergy, "discrete element strain energy");
	REGISTER_FECORE_CLASS(FEPlotContinuousDamage_D    , "continuous damage");
	REGISTER_FECORE_CLASS(FEPlotContinuousDamage_D1   , "continuous damage D1");
	REGISTER_FECORE_CLASS(FEPlotContinuousDamage_Ds   , "continuous damage Ds");
	REGISTER_FECORE_CLASS(FEPlotContinuousDamage_D2   , "continuous damage D2");
	REGISTER_FECORE_CLASS(FEPlotContinuousDamage_D3   , "continuous damage D3");
	REGISTER_FECORE_CLASS(FEPlotContinuousDamage_P    , "continuous damage P");
	REGISTER_FECORE_CLASS(FEPlotContinuousDamage_Psi0 , "continuous damage Psi0");
	REGISTER_FECORE_CLASS(FEPlotContinuousDamage_beta , "continuous damage beta");
	REGISTER_FECORE_CLASS(FEPlotContinuousDamage_gamma, "continuous damage gamma");
	REGISTER_FECORE_CLASS(FEPlotContinuousDamage_D2beta, "continuous damage D2beta");
    REGISTER_FECORE_CLASS(FEPlotRVEgenerations, "RVE generations");
    REGISTER_FECORE_CLASS(FEPlotRVEbonds, "RVE reforming bonds");
    REGISTER_FECORE_CLASS(FEPlotRVErecruitment, "RVE recruitment");
    REGISTER_FECORE_CLASS(FEPlotRVEstrain, "RVE strain");
    REGISTER_FECORE_CLASS(FEPlotStrongBondSED, "strong bond SED");
    REGISTER_FECORE_CLASS(FEPlotWeakBondSED, "weak bond SED");
    REGISTER_FECORE_CLASS(FEPlotStrongBondDevSED, "deviatoric strong bond SED");
    REGISTER_FECORE_CLASS(FEPlotWeakBondDevSED, "deviatoric weak bond SED");
    REGISTER_FECORE_CLASS(FEPlotTrussStretch  , "truss stretch");
    REGISTER_FECORE_CLASS(FEPlotGrowthLagrangeStrain, "growth Lagrange strain");
    REGISTER_FECORE_CLASS(FEPlotGrowthInfStrain, "growth infinitesimal strain");
    REGISTER_FECORE_CLASS(FEPlotGrowthRightStretch, "growth right stretch");
    REGISTER_FECORE_CLASS(FEPlotGrowthLeftStretch, "growth left stretch");
    REGISTER_FECORE_CLASS(FEPlotGrowthRightHencky, "growth right Hencky");
    REGISTER_FECORE_CLASS(FEPlotGrowthLeftHencky, "growth left Hencky");
    REGISTER_FECORE_CLASS(FEPlotGrowthRelativeVolume, "growth relative volume");

	// beam variables
	REGISTER_FECORE_CLASS(FEPlotBeamStress      , "beam stress");
	REGISTER_FECORE_CLASS(FEPlotBeamStressCouple, "beam stress couple");
	REGISTER_FECORE_CLASS(FEPlotBeamStrain      , "beam strain");
	REGISTER_FECORE_CLASS(FEPlotBeamCurvature   , "beam curvature");

	REGISTER_FECORE_CLASS(FEPlotBeamReferenceStress      , "beam reference stress");
	REGISTER_FECORE_CLASS(FEPlotBeamReferenceStressCouple, "beam reference stress couple");

	// 2O continuum fields
	REGISTER_FECORE_CLASS(FEPlotElementsnorm, "s norm");

	//-----------------------------------------------------------------------------
	// Derived from FELogNodeData
	REGISTER_FECORE_CLASS(FENodeXPos, "x");
	REGISTER_FECORE_CLASS(FENodeYPos, "y");
	REGISTER_FECORE_CLASS(FENodeZPos, "z");
	REGISTER_FECORE_CLASS(FENodeXDisp, "ux");
	REGISTER_FECORE_CLASS(FENodeYDisp, "uy");
	REGISTER_FECORE_CLASS(FENodeZDisp, "uz");
	REGISTER_FECORE_CLASS(FENodeXVel, "vx");
	REGISTER_FECORE_CLASS(FENodeYVel, "vy");
	REGISTER_FECORE_CLASS(FENodeZVel, "vz");
	REGISTER_FECORE_CLASS(FENodeXAcc, "ax");
	REGISTER_FECORE_CLASS(FENodeYAcc, "ay");
	REGISTER_FECORE_CLASS(FENodeZAcc, "az");
	REGISTER_FECORE_CLASS(FENodeForceX, "Rx");
	REGISTER_FECORE_CLASS(FENodeForceY, "Ry");
	REGISTER_FECORE_CLASS(FENodeForceZ, "Rz");

	//-----------------------------------------------------------------------------
	// Derived from FELogFaceData
	REGISTER_FECORE_CLASS(FELogContactGap     , "contact gap");
	REGISTER_FECORE_CLASS(FELogContactPressure, "contact pressure");
	REGISTER_FECORE_CLASS(FELogContactTractionX, "contact_traction.x");
	REGISTER_FECORE_CLASS(FELogContactTractionY, "contact_traction.y");
	REGISTER_FECORE_CLASS(FELogContactTractionZ, "contact_traction.z");

	//-----------------------------------------------------------------------------
	// Derived from FELogElemData
	REGISTER_FECORE_CLASS(FELogElemPosX, "x");
	REGISTER_FECORE_CLASS(FELogElemPosY, "y");
	REGISTER_FECORE_CLASS(FELogElemPosZ, "z");
	REGISTER_FECORE_CLASS(FELogElemJacobian, "J");
	REGISTER_FECORE_CLASS(FELogElemStrainX, "Ex");
	REGISTER_FECORE_CLASS(FELogElemStrainY, "Ey");
	REGISTER_FECORE_CLASS(FELogElemStrainZ, "Ez");
	REGISTER_FECORE_CLASS(FELogElemStrainXY, "Exy");
	REGISTER_FECORE_CLASS(FELogElemStrainYZ, "Eyz");
	REGISTER_FECORE_CLASS(FELogElemStrainXZ, "Exz");
	REGISTER_FECORE_CLASS(FELogElemStrainEffective, "effective strain");
	REGISTER_FECORE_CLASS(FELogElemStrain1, "E1");
	REGISTER_FECORE_CLASS(FELogElemStrain2, "E2");
	REGISTER_FECORE_CLASS(FELogElemStrain3, "E3");
	REGISTER_FECORE_CLASS(FELogElemInfStrainX, "ex");
	REGISTER_FECORE_CLASS(FELogElemInfStrainY, "ey");
	REGISTER_FECORE_CLASS(FELogElemInfStrainZ, "ez");
	REGISTER_FECORE_CLASS(FELogElemInfStrainXY, "exy");
	REGISTER_FECORE_CLASS(FELogElemInfStrainYZ, "eyz");
	REGISTER_FECORE_CLASS(FELogElemInfStrainXZ, "exz");
    REGISTER_FECORE_CLASS(FELogElemRightStretchX, "Ux");
    REGISTER_FECORE_CLASS(FELogElemRightStretchY, "Uy");
    REGISTER_FECORE_CLASS(FELogElemRightStretchZ, "Uz");
    REGISTER_FECORE_CLASS(FELogElemRightStretchXY, "Uxy");
    REGISTER_FECORE_CLASS(FELogElemRightStretchYZ, "Uyz");
    REGISTER_FECORE_CLASS(FELogElemRightStretchXZ, "Uxz");
    REGISTER_FECORE_CLASS(FELogElemRightStretchEffective, "effective right stretch");
    REGISTER_FECORE_CLASS(FELogElemRightStretch1, "U1");
    REGISTER_FECORE_CLASS(FELogElemRightStretch2, "U2");
    REGISTER_FECORE_CLASS(FELogElemRightStretch3, "U3");
    REGISTER_FECORE_CLASS(FELogElemLeftStretchX, "Vx");
    REGISTER_FECORE_CLASS(FELogElemLeftStretchY, "Vy");
    REGISTER_FECORE_CLASS(FELogElemLeftStretchZ, "Vz");
    REGISTER_FECORE_CLASS(FELogElemLeftStretchXY, "Vxy");
    REGISTER_FECORE_CLASS(FELogElemLeftStretchYZ, "Vyz");
    REGISTER_FECORE_CLASS(FELogElemLeftStretchXZ, "Vxz");
    REGISTER_FECORE_CLASS(FELogElemLeftStretchEffective, "effective left stretch");
    REGISTER_FECORE_CLASS(FELogElemLeftStretch1, "V1");
    REGISTER_FECORE_CLASS(FELogElemLeftStretch2, "V2");
    REGISTER_FECORE_CLASS(FELogElemLeftStretch3, "V3");
    REGISTER_FECORE_CLASS(FELogElemRightHenckyX, "Hx");
    REGISTER_FECORE_CLASS(FELogElemRightHenckyY, "Hy");
    REGISTER_FECORE_CLASS(FELogElemRightHenckyZ, "Hz");
    REGISTER_FECORE_CLASS(FELogElemRightHenckyXY, "Hxy");
    REGISTER_FECORE_CLASS(FELogElemRightHenckyYZ, "Hyz");
    REGISTER_FECORE_CLASS(FELogElemRightHenckyXZ, "Hxz");
    REGISTER_FECORE_CLASS(FELogElemRightHenckyEffective, "effective right Hencky");
    REGISTER_FECORE_CLASS(FELogElemRightHencky1, "H1");
    REGISTER_FECORE_CLASS(FELogElemRightHencky2, "H2");
    REGISTER_FECORE_CLASS(FELogElemRightHencky3, "H3");
    REGISTER_FECORE_CLASS(FELogElemLeftHenckyX, "hx");
    REGISTER_FECORE_CLASS(FELogElemLeftHenckyY, "hy");
    REGISTER_FECORE_CLASS(FELogElemLeftHenckyZ, "hz");
    REGISTER_FECORE_CLASS(FELogElemLeftHenckyXY, "hxy");
    REGISTER_FECORE_CLASS(FELogElemLeftHenckyYZ, "hyz");
    REGISTER_FECORE_CLASS(FELogElemLeftHenckyXZ, "hxz");
    REGISTER_FECORE_CLASS(FELogElemLeftHenckyEffective, "effective left Hencky");
    REGISTER_FECORE_CLASS(FELogElemLeftHencky1, "h1");
    REGISTER_FECORE_CLASS(FELogElemLeftHencky2, "h2");
    REGISTER_FECORE_CLASS(FELogElemLeftHencky3, "h3");
	REGISTER_FECORE_CLASS(FELogElemStressX, "sx");
	REGISTER_FECORE_CLASS(FELogElemStressY, "sy");
	REGISTER_FECORE_CLASS(FELogElemStressZ, "sz");
	REGISTER_FECORE_CLASS(FELogElemStressXY, "sxy");
	REGISTER_FECORE_CLASS(FELogElemStressYZ, "syz");
	REGISTER_FECORE_CLASS(FELogElemStressXZ, "sxz");
	REGISTER_FECORE_CLASS(FELogElemStressEffective, "effective stress");
	REGISTER_FECORE_CLASS(FELogElemStress1, "s1");
	REGISTER_FECORE_CLASS(FELogElemStress2, "s2");
	REGISTER_FECORE_CLASS(FELogElemStress3, "s3");
	REGISTER_FECORE_CLASS(FELogElemPK2StressX , "Sx");
	REGISTER_FECORE_CLASS(FELogElemPK2StressY , "Sy");
	REGISTER_FECORE_CLASS(FELogElemPK2StressZ , "Sz");
	REGISTER_FECORE_CLASS(FELogElemPK2StressXY, "Sxy");
	REGISTER_FECORE_CLASS(FELogElemPK2StressYZ, "Syz");
	REGISTER_FECORE_CLASS(FELogElemPK2StressXZ, "Sxz");
	REGISTER_FECORE_CLASS_T2(FELogElemStressEigenVector_T, 0, 0, "s1x");
	REGISTER_FECORE_CLASS_T2(FELogElemStressEigenVector_T, 0, 1, "s1y");
	REGISTER_FECORE_CLASS_T2(FELogElemStressEigenVector_T, 0, 2, "s1z");
	REGISTER_FECORE_CLASS_T2(FELogElemStressEigenVector_T, 1, 0, "s2x");
	REGISTER_FECORE_CLASS_T2(FELogElemStressEigenVector_T, 1, 1, "s2y");
	REGISTER_FECORE_CLASS_T2(FELogElemStressEigenVector_T, 1, 2, "s2z");
	REGISTER_FECORE_CLASS_T2(FELogElemStressEigenVector_T, 2, 0, "s3x");
	REGISTER_FECORE_CLASS_T2(FELogElemStressEigenVector_T, 2, 1, "s3y");
	REGISTER_FECORE_CLASS_T2(FELogElemStressEigenVector_T, 2, 2, "s3z");
	REGISTER_FECORE_CLASS(FELogElemDeformationGradientXX, "Fxx");
	REGISTER_FECORE_CLASS(FELogElemDeformationGradientXY, "Fxy");
	REGISTER_FECORE_CLASS(FELogElemDeformationGradientXZ, "Fxz");
	REGISTER_FECORE_CLASS(FELogElemDeformationGradientYX, "Fyx");
	REGISTER_FECORE_CLASS(FELogElemDeformationGradientYY, "Fyy");
	REGISTER_FECORE_CLASS(FELogElemDeformationGradientYZ, "Fyz");
	REGISTER_FECORE_CLASS(FELogElemDeformationGradientZX, "Fzx");
	REGISTER_FECORE_CLASS(FELogElemDeformationGradientZY, "Fzy");
	REGISTER_FECORE_CLASS(FELogElemDeformationGradientZZ, "Fzz");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 0, "cxxxx");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 1, "cxxyy");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 2, "cyyyy");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 3, "cxxzz");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 4, "cyyzz");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 5, "czzzz");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 6, "cxxxy");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 7, "cyyxy");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 8, "czzxy");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 9, "cxyxy");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 10, "cxxyz");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 11, "cyyyz");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 12, "czzyz");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 13, "cxyyz");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 14, "cyzyz");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 15, "cxxxz");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 16, "cyyxz");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 17, "czzxz");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 18, "cxyxz");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 19, "cyzxz");
	REGISTER_FECORE_CLASS_T(FELogElemElasticity, 20, "cxzxz");
	REGISTER_FECORE_CLASS(FELogElemStrainEnergyDensity, "sed");
	REGISTER_FECORE_CLASS(FELogElemDevStrainEnergyDensity, "devsed");
	REGISTER_FECORE_CLASS(FELogElemFiberStretch, "fiber_stretch");
	REGISTER_FECORE_CLASS(FELogElemFiberVectorX, "fiber_x");
	REGISTER_FECORE_CLASS(FELogElemFiberVectorY, "fiber_y");
	REGISTER_FECORE_CLASS(FELogElemFiberVectorZ, "fiber_z");
	REGISTER_FECORE_CLASS(FELogDamage, "D");
    REGISTER_FECORE_CLASS(FELogIntactBonds, "wi");
    REGISTER_FECORE_CLASS(FELogYieldedBonds, "wy");
    REGISTER_FECORE_CLASS(FELogFatigueBonds, "wf");
    REGISTER_FECORE_CLASS(FELogOctahedralPlasticStrain, "ops");
    REGISTER_FECORE_CLASS(FELogDiscreteElementStretch   , "discrete element stretch");
    REGISTER_FECORE_CLASS(FELogDiscreteElementElongation, "discrete element elongation");
    REGISTER_FECORE_CLASS(FELogDiscreteElementForce     , "discrete element force"  );
    REGISTER_FECORE_CLASS(FELogDiscreteElementForceX    , "Fde.x");
    REGISTER_FECORE_CLASS(FELogDiscreteElementForceY    , "Fde.y");
    REGISTER_FECORE_CLASS(FELogDiscreteElementForceZ    , "Fde.z");
	REGISTER_FECORE_CLASS(FELogContactArea, "contact area");
	REGISTER_FECORE_CLASS_T2(FELogElementMixtureStress_T, 0, 0, "mixture_stress[0].xx");
	REGISTER_FECORE_CLASS_T2(FELogElementMixtureStress_T, 0, 1, "mixture_stress[0].xy");
	REGISTER_FECORE_CLASS_T2(FELogElementMixtureStress_T, 0, 2, "mixture_stress[0].yy");
	REGISTER_FECORE_CLASS_T2(FELogElementMixtureStress_T, 0, 3, "mixture_stress[0].xz");
	REGISTER_FECORE_CLASS_T2(FELogElementMixtureStress_T, 0, 4, "mixture_stress[0].yz");
	REGISTER_FECORE_CLASS_T2(FELogElementMixtureStress_T, 0, 5, "mixture_stress[0].zz");
	REGISTER_FECORE_CLASS_T2(FELogElementMixtureStress_T, 1, 0, "mixture_stress[1].xx");
	REGISTER_FECORE_CLASS_T2(FELogElementMixtureStress_T, 1, 1, "mixture_stress[1].xy");
	REGISTER_FECORE_CLASS_T2(FELogElementMixtureStress_T, 1, 2, "mixture_stress[1].yy");
	REGISTER_FECORE_CLASS_T2(FELogElementMixtureStress_T, 1, 3, "mixture_stress[1].xz");
	REGISTER_FECORE_CLASS_T2(FELogElementMixtureStress_T, 1, 4, "mixture_stress[1].yz");
	REGISTER_FECORE_CLASS_T2(FELogElementMixtureStress_T, 1, 5, "mixture_stress[1].zz");
	REGISTER_FECORE_CLASS_T2(FELogElementMixtureStress_T, 2, 0, "mixture_stress[2].xx");
	REGISTER_FECORE_CLASS_T2(FELogElementMixtureStress_T, 2, 1, "mixture_stress[2].xy");
	REGISTER_FECORE_CLASS_T2(FELogElementMixtureStress_T, 2, 2, "mixture_stress[2].yy");
	REGISTER_FECORE_CLASS_T2(FELogElementMixtureStress_T, 2, 3, "mixture_stress[2].xz");
	REGISTER_FECORE_CLASS_T2(FELogElementMixtureStress_T, 2, 4, "mixture_stress[2].yz");
	REGISTER_FECORE_CLASS_T2(FELogElementMixtureStress_T, 2, 5, "mixture_stress[2].zz");

	//-----------------------------------------------------------------------------
	// Derived from FELogObjectData
	REGISTER_FECORE_CLASS(FELogRigidBodyPosX, "x");
	REGISTER_FECORE_CLASS(FELogRigidBodyPosY, "y");
	REGISTER_FECORE_CLASS(FELogRigidBodyPosZ, "z");
	REGISTER_FECORE_CLASS(FELogRigidBodyVelX, "vx");
	REGISTER_FECORE_CLASS(FELogRigidBodyVelY, "vy");
	REGISTER_FECORE_CLASS(FELogRigidBodyVelZ, "vz");
	REGISTER_FECORE_CLASS(FELogRigidBodyAccX, "ax");
	REGISTER_FECORE_CLASS(FELogRigidBodyAccY, "ay");
	REGISTER_FECORE_CLASS(FELogRigidBodyAccZ, "az");
	REGISTER_FECORE_CLASS(FELogRigidBodyAngPosX, "thx");
	REGISTER_FECORE_CLASS(FELogRigidBodyAngPosY, "thy");
	REGISTER_FECORE_CLASS(FELogRigidBodyAngPosZ, "thz");
	REGISTER_FECORE_CLASS(FELogRigidBodyAngVelX, "omx");
	REGISTER_FECORE_CLASS(FELogRigidBodyAngVelY, "omy");
	REGISTER_FECORE_CLASS(FELogRigidBodyAngVelZ, "omz");
	REGISTER_FECORE_CLASS(FELogRigidBodyAngAccX, "alx");
	REGISTER_FECORE_CLASS(FELogRigidBodyAngAccY, "aly");
	REGISTER_FECORE_CLASS(FELogRigidBodyAngAccZ, "alz");
	REGISTER_FECORE_CLASS(FELogRigidBodyQuatX, "qx");
	REGISTER_FECORE_CLASS(FELogRigidBodyQuatY, "qy");
	REGISTER_FECORE_CLASS(FELogRigidBodyQuatZ, "qz");
	REGISTER_FECORE_CLASS(FELogRigidBodyQuatW, "qw");
	REGISTER_FECORE_CLASS(FELogRigidBodyR11, "R11");
	REGISTER_FECORE_CLASS(FELogRigidBodyR12, "R12");
	REGISTER_FECORE_CLASS(FELogRigidBodyR13, "R13");
	REGISTER_FECORE_CLASS(FELogRigidBodyR21, "R21");
	REGISTER_FECORE_CLASS(FELogRigidBodyR22, "R22");
	REGISTER_FECORE_CLASS(FELogRigidBodyR23, "R23");
	REGISTER_FECORE_CLASS(FELogRigidBodyR31, "R31");
	REGISTER_FECORE_CLASS(FELogRigidBodyR32, "R32");
	REGISTER_FECORE_CLASS(FELogRigidBodyR33, "R33");
	REGISTER_FECORE_CLASS(FELogRigidBodyEulerX, "EulerX");
	REGISTER_FECORE_CLASS(FELogRigidBodyEulerY, "EulerY");
	REGISTER_FECORE_CLASS(FELogRigidBodyEulerZ, "EulerZ");
	REGISTER_FECORE_CLASS(FELogRigidBodyForceX, "Fx");
	REGISTER_FECORE_CLASS(FELogRigidBodyForceY, "Fy");
	REGISTER_FECORE_CLASS(FELogRigidBodyForceZ, "Fz");
	REGISTER_FECORE_CLASS(FELogRigidBodyTorqueX, "Mx");
	REGISTER_FECORE_CLASS(FELogRigidBodyTorqueY, "My");
	REGISTER_FECORE_CLASS(FELogRigidBodyTorqueZ, "Mz");
	REGISTER_FECORE_CLASS(FELogRigidBodyKineticEnergy, "KE");

	//-----------------------------------------------------------------------------
	// Derived from FELogConnectorData
	REGISTER_FECORE_CLASS(FELogRigidConnectorForceX, "RCFx");
	REGISTER_FECORE_CLASS(FELogRigidConnectorForceY, "RCFy");
	REGISTER_FECORE_CLASS(FELogRigidConnectorForceZ, "RCFz");
	REGISTER_FECORE_CLASS(FELogRigidConnectorMomentX, "RCMx");
	REGISTER_FECORE_CLASS(FELogRigidConnectorMomentY, "RCMy");
	REGISTER_FECORE_CLASS(FELogRigidConnectorMomentZ, "RCMz");
	REGISTER_FECORE_CLASS(FELogRigidConnectorTranslationX, "RCx");
	REGISTER_FECORE_CLASS(FELogRigidConnectorTranslationY, "RCy");
	REGISTER_FECORE_CLASS(FELogRigidConnectorTranslationZ, "RCz");
	REGISTER_FECORE_CLASS(FELogRigidConnectorRotationX, "RCthx");
	REGISTER_FECORE_CLASS(FELogRigidConnectorRotationY, "RCthy");
	REGISTER_FECORE_CLASS(FELogRigidConnectorRotationZ, "RCthz");

	//-----------------------------------------------------------------------------
	// Derived from FELogNLConstraintData
	REGISTER_FECORE_CLASS(FELogVolumeConstraint, "constrained volume");
	REGISTER_FECORE_CLASS(FELogVolumePressure, "volume pressure");

	//-----------------------------------------------------------------------------
	// Derived from DataRecord
	REGISTER_FECORE_CLASS(ObjectDataRecord, "rigid_body_data");

	//-----------------------------------------------------------------------------
	// Derived from FEMeshAdaptorCriterion
	REGISTER_FECORE_CLASS(FEStressCriterion, "stress");
	REGISTER_FECORE_CLASS(FEDamageAdaptorCriterion, "damage");
	REGISTER_FECORE_CLASS(FESpringForceCriterion, "spring force");
	REGISTER_FECORE_CLASS(FESpringStretchCriterion, "spring stretch");

	//-----------------------------------------------------------------------------
	// Derived from FEElemDataGenerator
	REGISTER_FECORE_CLASS(FEDeformationMapGenerator, "defgrad");

	//-----------------------------------------------------------------------------
	// Model update requests
	febio.OnCreateEvent(UpdateModelWhenCreating<FESolidAnalysis>([](FEModelUpdate& fem) {
			fem.AddPlotVariable("displacement");
			fem.AddPlotVariable("stress");
		}));

	febio.OnCreateEvent(AddPlotVariableWhenCreating<FEContactInterface>("contact pressure"));
	febio.OnCreateEvent(AddPlotVariableWhenCreating<FEContactInterface>("contact gap"));

	febio.SetActiveModule(0);
}
