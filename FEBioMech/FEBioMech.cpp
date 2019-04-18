/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FE2DFiberNeoHookean.h"
#include "FE2DTransIsoMooneyRivlin.h"
#include "FE2DTransIsoVerondaWestmann.h"
#include "FEArrudaBoyce.h"
#include "FECarterHayesOld.h"
#include "FECellGrowth.h"
#include "FECubicCLE.h"
#include "FEDamageMooneyRivlin.h"
#include "FEDamageNeoHookean.h"
#include "FEDamageTransIsoMooneyRivlin.h"
#include "FESpringMaterial.h"
#include "FEDonnanEquilibrium.h"
#include "FEEFDDonnanEquilibrium.h"
#include "FEEFDMooneyRivlin.h"
#include "FEEFDNeoHookean.h"
#include "FEEFDUncoupled.h"
#include "FEEFDVerondaWestmann.h"
#include "FEElasticMixture.h"
#include "FEElasticMultigeneration.h"
#include "FEEllipsoidalFiberDistribution.h"
#include "FEFatigueMaterial.h"
#include "FEFiberExpPow.h"
#include "FEFiberExpPowUncoupled.h"
#include "FEFiberNeoHookean.h"
#include "FEFiberPowLinear.h"
#include "FEFiberPowLinearUncoupled.h"
#include "FEFiberEFDNeoHookean.h"
#include "FEFiberExponentialPowerUC.h"
#include "FEFiberNHUC.h"
#include "FEFungOrthoCompressible.h"
#include "FEFungOrthotropic.h"
#include "FEGasserOgdenHolzapfel.h"
#include "FEGasserOgdenHolzapfelUC.h"
#include "FEHolmesMow.h"
#include "FEHuiskesSupply.h"
#include "FEIncompNeoHookean.h"
#include "FEIsotropicElastic.h"
#include "FEMooneyRivlin.h"
#include "FEMRVonMisesFibers.h"
#include "FEMuscleMaterial.h"
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
#include "FERemodelingElasticMaterial.h"
#include "FERigidMaterial.h"
#include "FESphericalFiberDistribution.h"
#include "FEStVenantKirchhoff.h"
#include "FETCNonlinearOrthotropic.h"
#include "FETendonMaterial.h"
#include "FETransIsoMooneyRivlin.h"
#include "FETransIsoVerondaWestmann.h"
#include "FETrussMaterial.h"
#include "FEUncoupledActiveContraction.h"
#include "FEUncoupledElasticMixture.h"
#include "FEUncoupledViscoElasticMaterial.h"
#include "FEVerondaWestmann.h"
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
#include "FEMicroMaterial.h"
#include "FEMicroMaterial2O.h"
#include "FESpringMaterial.h"
#include "FEElasticMultigeneration.h"
#include "FEPRLig.h"
#include "FECoupledMooneyRivlin.h"
#include "FECoupledVerondaWestmann.h"
#include "FEReactiveViscoelastic.h"
#include "FEUncoupledReactiveViscoelastic.h"
#include "FEBondRelaxation.h"
#include "FEDamageMaterial.h"
#include "FEDamageMaterialUC.h"
#include "FEDamageCDF.h"
#include "FEDamageCriterion.h"
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
#include "FEMindlinElastic2O.h"

#include "FEPressureLoad.h"
#include "FETractionLoad.h"
#include "FEGenericBodyForce.h"
#include "FECentrifugalBodyForce.h"
#include "FEPointBodyForce.h"

#include "FEFacet2FacetSliding.h"
#include "FEPeriodicBoundary.h"
#include "FEPeriodicBoundary2O.h"
#include "FERigidWallInterface.h"
#include "FERigidSlidingContact.h"
#include "FESlidingInterface.h"
#include "FESlidingInterfaceBW.h"
#include "FEPeriodicSurfaceConstraint.h"
#include "FETiedInterface.h"
#include "FETiedElasticInterface.h"
#include "FEStickyInterface.h"
#include "FEPointConstraint.h"
#include "FEFacet2FacetTied.h"
#include "FEVolumeConstraint.h"
#include "FEDistanceConstraint.h"
#include "FE2OMicroConstraint.h"
#include "FEMortarSlidingContact.h"
#include "FEMortarTiedContact.h"

#include "FEAugLagLinearConstraint.h"
#include "FESymmetryPlane.h"
#include "FERigidJoint.h"
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
#include "FERigidSolidDomain.h"
#include "FERigidShellDomain.h"
#include "FERemodelingElasticDomain.h"
#include "FEUDGHexDomain.h"
#include "FEUT4Domain.h"
#include "FEElasticMultiscaleDomain1O.h"
#include "FEElasticMultiscaleDomain2O.h"
#include "FESRIElasticSolidDomain.h"
#include "FE3FieldElasticSolidDomain.h"
#include "FE3FieldElasticShellDomain.h"
#include "FEDiscreteSpringDomain.h"
#include "FEDeformableSpringDomain.h"
#include "RigidBC.h"
#include "FEBCPrescribedDeformation.h"
#include "FEPrescribedNormalDisplacement.h"
#include "FEMaxStressCriterion.h"
#include "FEMaxDamageCriterion.h"
#include "FEInitialVelocity.h"

//-----------------------------------------------------------------------------
const char* FEBioMech::GetVariableName(FEBioMech::MECH_VARIABLE var)
{
	switch (var)
	{
	case DISPLACEMENT           : return "displacement"               ; break;
	case ROTATION               : return "rotation"                   ; break;
	case RIGID_ROTATION         : return "rigid rotation"             ; break;
	case SHELL_DISPLACEMENT     : return "shell displacement"         ; break;
	case VELOCTIY               : return "velocity"                   ; break;
	case PREV_ROTATION          : return "previous rotation"          ; break;
	case SHELL_VELOCITY         : return "shell velocity"             ; break;
	case PREV_SHELL_DISPLACEMENT: return "previous shell displacement"; break;
	case SHELL_ACCELERATION     : return "shell acceleration"         ; break;
	case PREV_SHELL_VELOCITY    : return "previous shell velocity"    ; break;
	case PREV_SHELL_ACCELERATION: return "previous_shell_acceleration"; break;
	}

	assert(false);
	return nullptr;
}


//-----------------------------------------------------------------------------
//! Register all the classes of the FEBioMech module with the FEBio framework.
void FEBioMech::InitModule()
{
//-----------------------------------------------------------------------------
// Domain factory
	FECoreKernel& febio = FECoreKernel::GetInstance();
	febio.RegisterDomain(new FESolidDomainFactory);

//-----------------------------------------------------------------------------
// Solver classes
REGISTER_FECORE_CLASS(FESolidSolver        , "solid_old"     );
REGISTER_FECORE_CLASS(FESolidSolver2       , "solid"         );
REGISTER_FECORE_CLASS(FEExplicitSolidSolver, "explicit-solid");
REGISTER_FECORE_CLASS(FECGSolidSolver      , "CG-solid"      );

//-----------------------------------------------------------------------------
// material classes

// elastic materials (derived from FEElasticMaterial)
REGISTER_FECORE_CLASS(FE2DFiberNeoHookean              , "2D fiber neo-Hookean"                );
REGISTER_FECORE_CLASS(FECellGrowth                     , "cell growth"                         );
REGISTER_FECORE_CLASS(FECubicCLE                       , "cubic CLE"                           );
REGISTER_FECORE_CLASS(FEDamageNeoHookean               , "damage neo-Hookean"                  );
REGISTER_FECORE_CLASS(FEDonnanEquilibrium              , "Donnan equilibrium"                  );
REGISTER_FECORE_CLASS(FEEFDDonnanEquilibrium           , "EFD Donnan equilibrium"              );
REGISTER_FECORE_CLASS(FEEFDNeoHookean                  , "EFD neo-Hookean (new)"               );
REGISTER_FECORE_CLASS(FEEFDNeoHookeanOld               , "EFD neo-Hookean"                     );
REGISTER_FECORE_CLASS(FEEllipsoidalFiberDistribution   , "ellipsoidal fiber distribution"      );
REGISTER_FECORE_CLASS(FEEllipsoidalFiberDistributionOld, "ellipsoidal fiber distribution (old)");
REGISTER_FECORE_CLASS(FEFungOrthoCompressible          , "Fung-ortho-compressible"             );
REGISTER_FECORE_CLASS(FEGasserOgdenHolzapfel           , "Gasser-Ogden-Holzapfel"              );
REGISTER_FECORE_CLASS(FECompressibleGentMaterial       , "compressible Gent"                   );
REGISTER_FECORE_CLASS(FEHolmesMow                      , "Holmes-Mow"                          );
REGISTER_FECORE_CLASS(FEIsotropicElastic               , "isotropic elastic"                   );
REGISTER_FECORE_CLASS(FECoupledMooneyRivlin            , "coupled Mooney-Rivlin"               );
REGISTER_FECORE_CLASS(FECoupledVerondaWestmann         , "coupled Veronda-Westmann"            );
REGISTER_FECORE_CLASS(FENeoHookean                     , "neo-Hookean"                         );
REGISTER_FECORE_CLASS(FENeoHookeanTransIso             , "neo-Hookean transiso"                );
REGISTER_FECORE_CLASS(FENewtonianViscousSolid          , "Newtonian viscous solid"             );
REGISTER_FECORE_CLASS(FEOgdenUnconstrained             , "Ogden unconstrained"                 );
REGISTER_FECORE_CLASS(FEOrthoElastic                   , "orthotropic elastic"                 );
REGISTER_FECORE_CLASS(FEOrthotropicCLE                 , "orthotropic CLE"                     );
REGISTER_FECORE_CLASS(FEOsmoticVirialExpansion         , "osmotic virial expansion"            );
REGISTER_FECORE_CLASS(FEPerfectOsmometer               , "perfect osmometer"                   );
REGISTER_FECORE_CLASS(FESphericalFiberDistribution     , "spherical fiber distribution"        );
REGISTER_FECORE_CLASS(FEStVenantKirchhoff              , "St.Venant-Kirchhoff"                 );
REGISTER_FECORE_CLASS(FEViscoElasticMaterial           , "viscoelastic"                        );
REGISTER_FECORE_CLASS(FEElasticMultigeneration         , "multigeneration"                     );
REGISTER_FECORE_CLASS(FERemodelingElasticMaterial      , "remodeling solid"                    );
REGISTER_FECORE_CLASS(FECarterHayesOld                 , "Carter-Hayes (old)"                  );
REGISTER_FECORE_CLASS(FEContinuousFiberDistribution    , "continuous fiber distribution"       );
REGISTER_FECORE_CLASS(FECoupledTransIsoVerondaWestmann , "coupled trans-iso Veronda-Westmann"  );
REGISTER_FECORE_CLASS(FECoupledTransIsoMooneyRivlin    , "coupled trans-iso Mooney-Rivlin"     );
REGISTER_FECORE_CLASS(FEMicroMaterial                  , "micro-material"                      );
REGISTER_FECORE_CLASS(FEMicroMaterial2O                , "micro-material2O"                    );
REGISTER_FECORE_CLASS(FEMindlinElastic2O               , "mindlin elastic"                     );

// These materials are derived from FEElasticMaterial and use FEElasticMaterials
REGISTER_FECORE_CLASS(FEElasticMixture                 , "solid mixture"                       );
REGISTER_FECORE_CLASS(FEGenerationMaterial             , "generation"                          );
REGISTER_FECORE_CLASS(FEReactiveViscoelasticMaterial   , "reactive viscoelastic"               );
REGISTER_FECORE_CLASS(FEDamageMaterial                 , "elastic damage"                      );
REGISTER_FECORE_CLASS(FEFatigueMaterial                , "reactive fatigue"                    );

// Uncoupled elastic materials (derived from FEUncoupledMaterial)
REGISTER_FECORE_CLASS(FEArrudaBoyce                          , "Arruda-Boyce"                           );
REGISTER_FECORE_CLASS(FE2DTransIsoMooneyRivlin               , "2D trans iso Mooney-Rivlin"             );
REGISTER_FECORE_CLASS(FE2DTransIsoVerondaWestmann            , "2D trans iso Veronda-Westmann"          );
REGISTER_FECORE_CLASS(FEDamageMooneyRivlin                   , "damage Mooney-Rivlin"                   );
REGISTER_FECORE_CLASS(FEDamageTransIsoMooneyRivlin           , "damage trans iso Mooney-Rivlin"         );
REGISTER_FECORE_CLASS(FEEFDMooneyRivlin                      , "EFD Mooney-Rivlin"                      );
REGISTER_FECORE_CLASS(FEEFDUncoupled                         , "EFD uncoupled"                          );
REGISTER_FECORE_CLASS(FEEFDVerondaWestmann                   , "EFD Veronda-Westmann"                   );
REGISTER_FECORE_CLASS(FEFungOrthotropic                      , "Fung orthotropic"                       );
REGISTER_FECORE_CLASS(FEGasserOgdenHolzapfelUC               , "Gasser-Ogden-Holzapfel-uncoupled"       );
REGISTER_FECORE_CLASS(FEGentMaterial                         , "Gent"                                   );
REGISTER_FECORE_CLASS(FEIncompNeoHookean                     , "incomp neo-Hookean"                     );
REGISTER_FECORE_CLASS(FEMooneyRivlin                         , "Mooney-Rivlin"                          );
REGISTER_FECORE_CLASS(FEMuscleMaterial                       , "muscle material"                        );
REGISTER_FECORE_CLASS(FENewtonianViscousSolidUC              , "Newtonian viscous solid uncoupled"      );
REGISTER_FECORE_CLASS(FEOgdenMaterial                        , "Ogden"                                  );
REGISTER_FECORE_CLASS(FETCNonlinearOrthotropic               , "TC nonlinear orthotropic"               );
REGISTER_FECORE_CLASS(FETendonMaterial                       , "tendon material"                        );
REGISTER_FECORE_CLASS(FETransIsoMooneyRivlin                 , "trans iso Mooney-Rivlin"                );
REGISTER_FECORE_CLASS(FETransIsoVerondaWestmann              , "trans iso Veronda-Westmann"             );
REGISTER_FECORE_CLASS(FEUncoupledElasticMixture              , "uncoupled solid mixture"                );
REGISTER_FECORE_CLASS(FEVerondaWestmann                      , "Veronda-Westmann"                       );
REGISTER_FECORE_CLASS(FEUncoupledViscoElasticMaterial        , "uncoupled viscoelastic"                 );
REGISTER_FECORE_CLASS(FEMRVonMisesFibers                     , "Mooney-Rivlin von Mises Fibers"         );
REGISTER_FECORE_CLASS(FEUncoupledActiveContraction           , "uncoupled active contraction"           );
REGISTER_FECORE_CLASS(FEContinuousFiberDistributionUC        , "continuous fiber distribution uncoupled");
REGISTER_FECORE_CLASS(FEPRLig					             , "PRLig"                                  );
REGISTER_FECORE_CLASS(FEUncoupledReactiveViscoelasticMaterial, "uncoupled reactive viscoelastic"        );
REGISTER_FECORE_CLASS(FEDamageMaterialUC                     , "uncoupled elastic damage"               );

// Fiber materials
REGISTER_FECORE_CLASS(FEFiberExpPow                  , "fiber-exp-pow"                        );
REGISTER_FECORE_CLASS(FEFiberExpPowUncoupled         , "fiber-exp-pow-uncoupled"              );
REGISTER_FECORE_CLASS(FEFiberEFDNeoHookean           , "fiber neo-Hookean"                    );
REGISTER_FECORE_CLASS(FEFiberPowLinear               , "fiber-pow-linear"                     );
REGISTER_FECORE_CLASS(FEFiberPowLinearUncoupled      , "fiber-pow-linear-uncoupled"           );
REGISTER_FECORE_CLASS(FEFiberExponentialPower        , "fiber-exponential-power-law"          );
REGISTER_FECORE_CLASS(FEFiberExponentialPowerUC      , "fiber-exponential-power-law-uncoupled");
REGISTER_FECORE_CLASS(FEFiberNH                      , "fiber-NH"                             );
REGISTER_FECORE_CLASS(FEFiberNHUC                    , "fiber-NH-uncoupled"                   );
REGISTER_FECORE_CLASS(FEFiberPowerLinear             , "fiber-power-linear"                   );
REGISTER_FECORE_CLASS(FEFiberExpLinear				 , "fiber-exp-linear"                     );
REGISTER_FECORE_CLASS(FEUncoupledFiberExpLinear      , "uncoupled fiber-exp-linear"           );

// solid materials (derived from FESolidMaterial)
REGISTER_FECORE_CLASS(FERigidMaterial                , "rigid body"          );
REGISTER_FECORE_CLASS(FEVonMisesPlasticity           , "von-Mises plasticity");

// Fiber density distributions for CFD materials
REGISTER_FECORE_CLASS(FESphericalFiberDensityDistribution  , "spherical"              );
REGISTER_FECORE_CLASS(FEEllipsodialFiberDensityDistribution, "ellipsoidal"            );
REGISTER_FECORE_CLASS(FEVonMises3DFiberDensityDistribution , "von-Mises-3d"           );
REGISTER_FECORE_CLASS(FEVonMises3DTwoFDDAxisymmetric       , "von-Mises-3d-two-axisym");
REGISTER_FECORE_CLASS(FECircularFiberDensityDistribution   , "circular"               );
REGISTER_FECORE_CLASS(FEEllipticalFiberDensityDistribution , "elliptical"             );
REGISTER_FECORE_CLASS(FEVonMises2DFiberDensityDistribution , "von-Mises-2d"           );

// Fiber distribution integration schemes for CFD materials
REGISTER_FECORE_CLASS(FEFiberIntegrationGauss              , "fibers-3d-gauss"      );
REGISTER_FECORE_CLASS(FEFiberIntegrationGeodesic           , "fibers-3d-geodesic"   );
REGISTER_FECORE_CLASS(FEFiberIntegrationGaussKronrod       , "fibers-3d-gkt"        );
REGISTER_FECORE_CLASS(FEFiberIntegrationTriangle           , "fibers-3d-fei"        );
REGISTER_FECORE_CLASS(FEFiberIntegrationTrapezoidal        , "fibers-2d-trapezoidal");

// Other materials 
REGISTER_FECORE_CLASS(FETrussMaterial         , "linear truss"      );
REGISTER_FECORE_CLASS(FEHuiskesSupply         , "Huiskes-supply"    );
REGISTER_FECORE_CLASS(FEActiveFiberContraction, "active_contraction");
REGISTER_FECORE_CLASS(FEMicroProbe            , "probe"             );
REGISTER_FECORE_CLASS(FEWrinkleOgdenMaterial  , "wrinkle Ogden"     );
REGISTER_FECORE_CLASS(FEElasticMembrane       , "elastic membrane"  );

// active contraction materials
REGISTER_FECORE_CLASS(FEPrescribedActiveContractionUniaxial   , "prescribed uniaxial active contraction");
REGISTER_FECORE_CLASS(FEPrescribedActiveContractionUniaxialUC , "uncoupled prescribed uniaxial active contraction");
REGISTER_FECORE_CLASS(FEPrescribedActiveContractionTransIso   , "prescribed trans iso active contraction");
REGISTER_FECORE_CLASS(FEPrescribedActiveContractionTransIsoUC , "uncoupled prescribed trans iso active contraction");
REGISTER_FECORE_CLASS(FEPrescribedActiveContractionIsotropic  , "prescribed isotropic active contraction");
REGISTER_FECORE_CLASS(FEPrescribedActiveContractionIsotropicUC, "uncoupled prescribed isotropic active contraction");

// spring materials
REGISTER_FECORE_CLASS(FELinearSpring           , "linear spring"             );
REGISTER_FECORE_CLASS(FETensionOnlyLinearSpring, "tension-only linear spring");
REGISTER_FECORE_CLASS(FENonLinearSpring        , "nonlinear spring"          );
REGISTER_FECORE_CLASS(FEExperimentalSpring     , "experimental spring"       );

// bond relaxation materials (used by reactive visco-elastic materials)
REGISTER_FECORE_CLASS(FEBondRelaxationExponential    , "relaxation-exponential"     );
REGISTER_FECORE_CLASS(FEBondRelaxationExpDistortion  , "relaxation-exp-distortion"  );
REGISTER_FECORE_CLASS(FEBondRelaxationFung           , "relaxation-Fung"            );
REGISTER_FECORE_CLASS(FEBondRelaxationPark           , "relaxation-Park"            );
REGISTER_FECORE_CLASS(FEBondRelaxationParkDistortion , "relaxation-Park-distortion" );
REGISTER_FECORE_CLASS(FEBondRelaxationPower          , "relaxation-power"           );
REGISTER_FECORE_CLASS(FEBondRelaxationPowerDistortion, "relaxation-power-distortion");
REGISTER_FECORE_CLASS(FEBondRelaxationCarreau        , "relaxation-Carreau"         );
    
// damage cumulative distribution functions (used by damage materials)
REGISTER_FECORE_CLASS(FEDamageCDFSimo      , "CDF Simo"      );
REGISTER_FECORE_CLASS(FEDamageCDFLogNormal , "CDF log-normal");
REGISTER_FECORE_CLASS(FEDamageCDFWeibull   , "CDF Weibull"   );
REGISTER_FECORE_CLASS(FEDamageCDFStep      , "CDF step"      );
REGISTER_FECORE_CLASS(FEDamageCDFPQP       , "CDF quintic"   );
REGISTER_FECORE_CLASS(FEDamageCriterionSimo, "DC Simo"       );

// damage criterion (used by damage materials)
REGISTER_FECORE_CLASS(FEDamageCriterionSED , "DC strain energy density"     );
REGISTER_FECORE_CLASS(FEDamageCriterionSSE , "DC specific strain energy"    );
REGISTER_FECORE_CLASS(FEDamageCriterionVMS , "DC von Mises stress"          );
REGISTER_FECORE_CLASS(FEDamageCriterionMSS , "DC max shear stress"          );
REGISTER_FECORE_CLASS(FEDamageCriterionMNS , "DC max normal stress"         );
REGISTER_FECORE_CLASS(FEDamageCriterionMNLS, "DC max normal Lagrange strain");

//-----------------------------------------------------------------------------
// domain classes
REGISTER_FECORE_CLASS(FERigidSolidDomain         , "rigid-solid"       );
REGISTER_FECORE_CLASS(FERigidShellDomain         , "rigid-shell"       );
REGISTER_FECORE_CLASS(FERigidShellDomainOld      , "rigid-shell-old"   );
REGISTER_FECORE_CLASS(FERemodelingElasticDomain  , "remodeling-solid"  );
REGISTER_FECORE_CLASS(FEElasticMultiscaleDomain1O, "elastic-mm-solid"  );
REGISTER_FECORE_CLASS(FEElasticMultiscaleDomain2O, "elastic-mm-solid2O");
REGISTER_FECORE_CLASS(FEElasticSolidDomain2O     , "elastic-solid2O"   );
REGISTER_FECORE_CLASS(FE3FieldElasticSolidDomain , "three-field-solid" );
REGISTER_FECORE_CLASS(FE3FieldElasticShellDomain , "three-field-shell" );
REGISTER_FECORE_CLASS(FEUDGHexDomain             , "udg-hex"           );
REGISTER_FECORE_CLASS(FESRIElasticSolidDomain    , "sri-solid"         );
REGISTER_FECORE_CLASS(FEUT4Domain                , "ut4-solid"         );
REGISTER_FECORE_CLASS(FEElasticSolidDomain       , "elastic-solid"     );
REGISTER_FECORE_CLASS(FEElasticShellDomain       , "elastic-shell"     );
REGISTER_FECORE_CLASS(FEElasticShellDomainOld    , "elastic-shell-old" );
REGISTER_FECORE_CLASS(FEElasticEASShellDomain    , "elastic-shell-eas" );
REGISTER_FECORE_CLASS(FEElasticANSShellDomain    , "elastic-shell-ans" );
REGISTER_FECORE_CLASS(FEElasticTrussDomain       , "elastic-truss"     );
REGISTER_FECORE_CLASS(FEDiscreteSpringDomain     , "discrete-spring"   );
REGISTER_FECORE_CLASS(FEDeformableSpringDomain   , "deformable-spring" );
REGISTER_FECORE_CLASS(FEDeformableSpringDomain2  , "deformable-spring2");

//-----------------------------------------------------------------------------
// classes derived from FEBoundaryCondition
REGISTER_FECORE_CLASS(FEBCPrescribedDeformation     , "prescribed deformation"   );
REGISTER_FECORE_CLASS(FEBCPrescribedDeformation2O   , "prescribed deformation 2O");
REGISTER_FECORE_CLASS(FEPrescribedNormalDisplacement, "normal displacement"      );

//-----------------------------------------------------------------------------
// classes derived from FEInitialCondition
REGISTER_FECORE_CLASS(FEInitialVelocity, "velocity");

//-----------------------------------------------------------------------------
// classes derived from FESurfaceLoad
REGISTER_FECORE_CLASS(FEPressureLoad, "pressure");
REGISTER_FECORE_CLASS(FETractionLoad, "traction");

//-----------------------------------------------------------------------------
// classes derived from FEBodyForce
REGISTER_FECORE_CLASS(FEConstBodyForceOld   , "const"    , 0x0205);	// obsolete in 3.0
REGISTER_FECORE_CLASS(FENonConstBodyForceOld, "non-const", 0x0205);	// obsolete in 3.0

REGISTER_FECORE_CLASS(FEGenericBodyForce    , "body force" );
REGISTER_FECORE_CLASS(FECentrifugalBodyForce, "centrifugal");
REGISTER_FECORE_CLASS(FEPointBodyForce      , "point"      );

//-----------------------------------------------------------------------------
// constraint classes
REGISTER_FECORE_CLASS(FEPointConstraint      , "point"                  );
REGISTER_FECORE_CLASS(FELinearConstraintSet  , "linear constraint"      );
REGISTER_FECORE_CLASS(FESymmetryPlane        , "symmetry plane"         );
REGISTER_FECORE_CLASS(FERigidJoint           , "rigid joint"            );
REGISTER_FECORE_CLASS(FERigidSphericalJoint  , "rigid spherical joint"  );
REGISTER_FECORE_CLASS(FERigidRevoluteJoint   , "rigid revolute joint"   );
REGISTER_FECORE_CLASS(FERigidPrismaticJoint  , "rigid prismatic joint"  );
REGISTER_FECORE_CLASS(FERigidCylindricalJoint, "rigid cylindrical joint");
REGISTER_FECORE_CLASS(FERigidPlanarJoint     , "rigid planar joint"     );
REGISTER_FECORE_CLASS(FERigidLock            , "rigid lock"             );
REGISTER_FECORE_CLASS(FERigidSpring          , "rigid spring"           );
REGISTER_FECORE_CLASS(FERigidDamper          , "rigid damper"           );
REGISTER_FECORE_CLASS(FERigidAngularDamper   , "rigid angular damper"   );
REGISTER_FECORE_CLASS(FERigidContractileForce, "rigid contractile force");
REGISTER_FECORE_CLASS(FEVolumeConstraint     , "volume"                 );
REGISTER_FECORE_CLASS(FEDiscreteContact      , "discrete contact"       );
REGISTER_FECORE_CLASS(FEDiscreteContact2     , "discrete contact2"      );
REGISTER_FECORE_CLASS(FEDistanceConstraint   , "node distance"          );
REGISTER_FECORE_CLASS(FE2OMicroConstraint    , "2O microfluc"           );

//-----------------------------------------------------------------------------
// classes derived from FEContactInterface

REGISTER_FECORE_CLASS(FESlidingInterface         , "sliding-node-on-facet"      );
REGISTER_FECORE_CLASS(FEFacet2FacetSliding       , "sliding-facet-on-facet"     );
REGISTER_FECORE_CLASS(FESlidingInterfaceBW       , "sliding-elastic"            );
REGISTER_FECORE_CLASS(FETiedInterface            , "tied-node-on-facet"         );
REGISTER_FECORE_CLASS(FEFacet2FacetTied          , "tied-facet-on-facet"        );
REGISTER_FECORE_CLASS(FETiedElasticInterface     , "tied-elastic"               );
REGISTER_FECORE_CLASS(FEPeriodicBoundary         , "periodic boundary"          );
REGISTER_FECORE_CLASS(FEPeriodicBoundary1O       , "periodic boundary1O"        );
REGISTER_FECORE_CLASS(FEPeriodicBoundary2O       , "periodic boundary2O"        );
REGISTER_FECORE_CLASS(FERigidWallInterface       , "rigid_wall"                 );
REGISTER_FECORE_CLASS(FERigidSlidingContact      , "rigid sliding"              );
REGISTER_FECORE_CLASS(FEPeriodicSurfaceConstraint, "surface constraint"         );
REGISTER_FECORE_CLASS(FEStickyInterface          , "sticky"                     );
REGISTER_FECORE_CLASS(FEMortarSlidingContact     , "mortar-sliding"             );
REGISTER_FECORE_CLASS(FEMortarTiedContact        , "mortar-tied"                );

//-----------------------------------------------------------------------------
// classes derived from FERigidSurface
REGISTER_FECORE_CLASS(FERigidPlane    , "plane"   );
REGISTER_FECORE_CLASS(FERigidSphere   , "sphere"  );
REGISTER_FECORE_CLASS(FERigidCylinder , "cylinder");
REGISTER_FECORE_CLASS(FERigidEllipsoid, "ellipsoid");

//-----------------------------------------------------------------------------
// classes derived directly from FEModelLoad
// TODO: define another SUPER_CLASS_ID for this
REGISTER_FECORE_CLASS_EXPLICIT(FERigidAxialForce      , FEBC_ID, "rigid_axial_force");
REGISTER_FECORE_CLASS_EXPLICIT(FERigidBodyForce       , FEBC_ID, "rigid_force"      );
REGISTER_FECORE_CLASS_EXPLICIT(FERigidBodyFixedBC     , FERIGIDBC_ID, "rigid_fixed"      );
REGISTER_FECORE_CLASS_EXPLICIT(FERigidBodyDisplacement, FERIGIDBC_ID, "rigid_prescribed" );
REGISTER_FECORE_CLASS_EXPLICIT(FERigidCable           , FEBC_ID, "rigid_cable"      );

//-----------------------------------------------------------------------------
// classes derived from FEPlotData
REGISTER_FECORE_CLASS(FEPlotElementVelocity              , "velocity"                        );
REGISTER_FECORE_CLASS(FEPlotElementAcceleration          , "acceleration"                    );
REGISTER_FECORE_CLASS(FEPlotDensity                      , "density"                         );
REGISTER_FECORE_CLASS(FEPlotElementStress                , "stress"                          );
REGISTER_FECORE_CLASS(FEPlotElementUncoupledPressure     , "uncoupled pressure"              );
REGISTER_FECORE_CLASS(FEPlotElementElasticity            , "elasticity"                      );
REGISTER_FECORE_CLASS(FEPlotRelativeVolume               , "relative volume"                 );
REGISTER_FECORE_CLASS(FEPlotFiberVector                  , "fiber vector"                    );
REGISTER_FECORE_CLASS(FEPlotFiberStretch                 , "fiber stretch"                   );
REGISTER_FECORE_CLASS(FEPlotDevFiberStretch              , "deviatoric fiber stretch"        );
REGISTER_FECORE_CLASS(FEPlotMaterialAxes		         , "material axes"                   );
REGISTER_FECORE_CLASS(FEPlotShellThickness               , "shell thickness"                 );
REGISTER_FECORE_CLASS(FEPlotShellDirector                , "shell director"                  );
REGISTER_FECORE_CLASS(FEPlotDamage                       , "damage"                          );
REGISTER_FECORE_CLASS(FEPlotNestedDamage                 , "nested damage"                   );
REGISTER_FECORE_CLASS(FEPlotMixtureVolumeFraction        , "volume fraction"                 );
REGISTER_FECORE_CLASS(FEPlotUT4NodalStresses             , "ut4 nodal stress"                );
REGISTER_FECORE_CLASS(FEPlotShellStrain                  , "shell strain"                    );
REGISTER_FECORE_CLASS(FEPlotShellRelativeVolume          , "shell relative volume"           );
REGISTER_FECORE_CLASS(FEPlotContactGap			         , "contact gap"                     );
REGISTER_FECORE_CLASS(FEPlotNodalContactGap              , "nodal contact gap"               );
REGISTER_FECORE_CLASS(FEPlotVectorGap			         , "vector gap"                      );
REGISTER_FECORE_CLASS(FEPlotNodalVectorGap               , "nodal vector gap"                );
REGISTER_FECORE_CLASS(FEPlotContactPressure		         , "contact pressure"                );
REGISTER_FECORE_CLASS(FEPlotNodalContactPressure         , "nodal contact pressure"          );
REGISTER_FECORE_CLASS(FEPlotContactTraction		         , "contact traction"                );
REGISTER_FECORE_CLASS(FEPlotNodalContactTraction         , "nodal contact traction"          );
REGISTER_FECORE_CLASS(FEPlotStickStatus			         , "contact stick"                   );
REGISTER_FECORE_CLASS(FEPlotContactForce 		         , "contact force"                   );
REGISTER_FECORE_CLASS(FEPlotContactArea 		         , "contact area"                    );
REGISTER_FECORE_CLASS(FEPlotContactPenalty 		         , "contact penalty"                 );
REGISTER_FECORE_CLASS(FEPlotSPRStresses                  , "SPR stress"                      );
REGISTER_FECORE_CLASS(FEPlotSPRLinearStresses            , "SPR-P1 stress"                   );
REGISTER_FECORE_CLASS(FEPlotSPRPrincStresses             , "SPR principal stress"            );
REGISTER_FECORE_CLASS(FEPlotNodalStresses		         , "nodal stress"		            );
REGISTER_FECORE_CLASS(FEPlotLagrangeStrain               , "Lagrange strain"                 );
REGISTER_FECORE_CLASS(FEPlotSPRLagrangeStrain            , "SPR Lagrange strain"             );
REGISTER_FECORE_CLASS(FEPlotMortarContactGap             , "mortar-gap"                      );
REGISTER_FECORE_CLASS(FEPlotSurfaceTraction		         , "surface traction"                );
REGISTER_FECORE_CLASS(FEPlotNodalSurfaceTraction         , "nodal surface traction"          );
REGISTER_FECORE_CLASS(FEPlotEnclosedVolume               , "enclosed volume"                 );
REGISTER_FECORE_CLASS(FEPlotStrainEnergyDensity          , "strain energy density"           );
REGISTER_FECORE_CLASS(FEPlotDevStrainEnergyDensity       , "deviatoric strain energy density");
REGISTER_FECORE_CLASS(FEPlotSpecificStrainEnergy         , "specific strain energy"          );
REGISTER_FECORE_CLASS(FEPlotKineticEnergyDensity         , "kinetic energy density"          );
REGISTER_FECORE_CLASS(FEPlotElementStrainEnergy          , "element strain energy"           );
REGISTER_FECORE_CLASS(FEPlotElementKineticEnergy         , "element kinetic energy"          );
REGISTER_FECORE_CLASS(FEPlotElementCenterOfMass          , "element center of mass"          );
REGISTER_FECORE_CLASS(FEPlotElementLinearMomentum        , "element linear momentum"         );
REGISTER_FECORE_CLASS(FEPlotElementAngularMomentum       , "element angular momentum"        );
REGISTER_FECORE_CLASS(FEPlotElementStressPower           , "element stress power"            );
REGISTER_FECORE_CLASS(FEPlotCurrentElementStrainEnergy   , "current element strain energy"   );
REGISTER_FECORE_CLASS(FEPlotCurrentElementKineticEnergy  , "current element kinetic energy"  );
REGISTER_FECORE_CLASS(FEPlotCurrentElementCenterOfMass   , "current element center of mass"  );
REGISTER_FECORE_CLASS(FEPlotCurrentElementLinearMomentum , "current element linear momentum" );
REGISTER_FECORE_CLASS(FEPlotCurrentElementAngularMomentum, "current element angular momentum");
REGISTER_FECORE_CLASS(FEPlotNodeVelocity                 , "nodal velocity"                  );
REGISTER_FECORE_CLASS(FEPlotNodeAcceleration             , "nodal acceleration"              );
REGISTER_FECORE_CLASS(FEPlotNodeReactionForces           , "reaction forces"                 );
REGISTER_FECORE_CLASS(FEPlotRigidReactionForce           , "rigid force"                     );
REGISTER_FECORE_CLASS(FEPlotRigidReactionTorque          , "rigid torque"                    );
REGISTER_FECORE_CLASS(FEPlotRigidDisplacement            , "rigid position"                  );
REGISTER_FECORE_CLASS(FEPlotRigidVelocity                , "rigid velocity"                  );
REGISTER_FECORE_CLASS(FEPlotRigidAcceleration            , "rigid acceleration"              );
REGISTER_FECORE_CLASS(FEPlotRigidRotation                , "rigid angular position"          );
REGISTER_FECORE_CLASS(FEPlotRigidAngularVelocity         , "rigid angular velocity"          );
REGISTER_FECORE_CLASS(FEPlotRigidAngularAcceleration     , "rigid angular acceleration"      );
REGISTER_FECORE_CLASS(FEPlotRigidLinearMomentum          , "rigid linear momentum"           );
REGISTER_FECORE_CLASS(FEPlotRigidAngularMomentum         , "rigid angular momentum"          );
REGISTER_FECORE_CLASS(FEPlotRigidKineticEnergy           , "rigid kinetic energy"            );
REGISTER_FECORE_CLASS(FEPlotRigidEuler                   , "Euler angle"                     );
REGISTER_FECORE_CLASS(FEPlotRigidRotationVector          , "rigid rotation vector"           );

// 2O continuum fields
REGISTER_FECORE_CLASS(FEPlotElementGnorm      , "G norm"      );
REGISTER_FECORE_CLASS(FEPlotElementsnorm      , "s norm"      );
REGISTER_FECORE_CLASS(FEPlotElementPK1norm    , "PK1 norm"    );
REGISTER_FECORE_CLASS(FEPlotElementQK1norm    , "QK1 norm"    );
REGISTER_FECORE_CLASS(FEPlotElementMicroEnergy, "micro energy");

//-----------------------------------------------------------------------------
// Derived from FENodeLogData
REGISTER_FECORE_CLASS(FENodeXPos  , "x");
REGISTER_FECORE_CLASS(FENodeYPos  , "y");
REGISTER_FECORE_CLASS(FENodeZPos  , "z");
REGISTER_FECORE_CLASS(FENodeXDisp , "ux");
REGISTER_FECORE_CLASS(FENodeYDisp , "uy");
REGISTER_FECORE_CLASS(FENodeZDisp , "uz");
REGISTER_FECORE_CLASS(FENodeXVel  , "vx");
REGISTER_FECORE_CLASS(FENodeYVel  , "vy");
REGISTER_FECORE_CLASS(FENodeZVel  , "vz");
REGISTER_FECORE_CLASS(FENodeXAcc  , "ax");
REGISTER_FECORE_CLASS(FENodeYAcc  , "ay");
REGISTER_FECORE_CLASS(FENodeZAcc  , "az");
REGISTER_FECORE_CLASS(FENodeForceX, "Rx");
REGISTER_FECORE_CLASS(FENodeForceY, "Ry");
REGISTER_FECORE_CLASS(FENodeForceZ, "Rz");

//-----------------------------------------------------------------------------
// Derived from FELogElemData
REGISTER_FECORE_CLASS(FELogElemPosX                 , "x");
REGISTER_FECORE_CLASS(FELogElemPosY                 , "y");
REGISTER_FECORE_CLASS(FELogElemPosZ                 , "z");
REGISTER_FECORE_CLASS(FELogElemJacobian             , "J");
REGISTER_FECORE_CLASS(FELogElemStrainX              , "Ex");
REGISTER_FECORE_CLASS(FELogElemStrainY              , "Ey");
REGISTER_FECORE_CLASS(FELogElemStrainZ              , "Ez");
REGISTER_FECORE_CLASS(FELogElemStrainXY             , "Exy");
REGISTER_FECORE_CLASS(FELogElemStrainYZ             , "Eyz");
REGISTER_FECORE_CLASS(FELogElemStrainXZ             , "Exz");
REGISTER_FECORE_CLASS(FELogElemStrain1              , "E1");
REGISTER_FECORE_CLASS(FELogElemStrain2              , "E2");
REGISTER_FECORE_CLASS(FELogElemStrain3              , "E3");
REGISTER_FECORE_CLASS(FELogElemInfStrainX           , "ex");
REGISTER_FECORE_CLASS(FELogElemInfStrainY           , "ey");
REGISTER_FECORE_CLASS(FELogElemInfStrainZ           , "ez");
REGISTER_FECORE_CLASS(FELogElemInfStrainXY          , "exy");
REGISTER_FECORE_CLASS(FELogElemInfStrainYZ          , "eyz");
REGISTER_FECORE_CLASS(FELogElemInfStrainXZ          , "exz");
REGISTER_FECORE_CLASS(FELogElemStressX              , "sx");
REGISTER_FECORE_CLASS(FELogElemStressY              , "sy");
REGISTER_FECORE_CLASS(FELogElemStressZ              , "sz");
REGISTER_FECORE_CLASS(FELogElemStressXY             , "sxy");
REGISTER_FECORE_CLASS(FELogElemStressYZ             , "syz");
REGISTER_FECORE_CLASS(FELogElemStressXZ             , "sxz");
REGISTER_FECORE_CLASS(FELogElemStress1              , "s1");
REGISTER_FECORE_CLASS(FELogElemStress2              , "s2");
REGISTER_FECORE_CLASS(FELogElemStress3              , "s3");
REGISTER_FECORE_CLASS(FELogElemDeformationGradientXX, "Fxx");
REGISTER_FECORE_CLASS(FELogElemDeformationGradientXY, "Fxy");
REGISTER_FECORE_CLASS(FELogElemDeformationGradientXZ, "Fxz");
REGISTER_FECORE_CLASS(FELogElemDeformationGradientYX, "Fyx");
REGISTER_FECORE_CLASS(FELogElemDeformationGradientYY, "Fyy");
REGISTER_FECORE_CLASS(FELogElemDeformationGradientYZ, "Fyz");
REGISTER_FECORE_CLASS(FELogElemDeformationGradientZX, "Fzx");
REGISTER_FECORE_CLASS(FELogElemDeformationGradientZY, "Fzy");
REGISTER_FECORE_CLASS(FELogElemDeformationGradientZZ, "Fzz");
REGISTER_FECORE_CLASS_T(FELogElemElasticity,  0, "cxxxx");
REGISTER_FECORE_CLASS_T(FELogElemElasticity,  1, "cxxyy");
REGISTER_FECORE_CLASS_T(FELogElemElasticity,  2, "cyyyy");
REGISTER_FECORE_CLASS_T(FELogElemElasticity,  3, "cxxzz");
REGISTER_FECORE_CLASS_T(FELogElemElasticity,  4, "cyyzz");
REGISTER_FECORE_CLASS_T(FELogElemElasticity,  5, "czzzz");
REGISTER_FECORE_CLASS_T(FELogElemElasticity,  6, "cxxxy");
REGISTER_FECORE_CLASS_T(FELogElemElasticity,  7, "cyyxy");
REGISTER_FECORE_CLASS_T(FELogElemElasticity,  8, "czzxy");
REGISTER_FECORE_CLASS_T(FELogElemElasticity,  9, "cxyxy");
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
REGISTER_FECORE_CLASS(FELogElemStrainEnergyDensity    , "sed");
REGISTER_FECORE_CLASS(FELogElemDevStrainEnergyDensity , "devsed");
REGISTER_FECORE_CLASS(FELogElemFiberStretch           , "fiber_stretch");
REGISTER_FECORE_CLASS(FELogElemFiberVectorX           , "fiber_x");
REGISTER_FECORE_CLASS(FELogElemFiberVectorY           , "fiber_y");
REGISTER_FECORE_CLASS(FELogElemFiberVectorZ           , "fiber_z");
REGISTER_FECORE_CLASS(FELogDamage                     , "D");

//-----------------------------------------------------------------------------
// Derived from FELogObjectData
REGISTER_FECORE_CLASS(FELogRigidBodyPosX         , "x");
REGISTER_FECORE_CLASS(FELogRigidBodyPosY         , "y");
REGISTER_FECORE_CLASS(FELogRigidBodyPosZ         , "z");
REGISTER_FECORE_CLASS(FELogRigidBodyVelX         , "vx");
REGISTER_FECORE_CLASS(FELogRigidBodyVelY         , "vy");
REGISTER_FECORE_CLASS(FELogRigidBodyVelZ         , "vz");
REGISTER_FECORE_CLASS(FELogRigidBodyAccX         , "ax");
REGISTER_FECORE_CLASS(FELogRigidBodyAccY         , "ay");
REGISTER_FECORE_CLASS(FELogRigidBodyAccZ         , "az");
REGISTER_FECORE_CLASS(FELogRigidBodyAngPosX      , "thx");
REGISTER_FECORE_CLASS(FELogRigidBodyAngPosY      , "thy");
REGISTER_FECORE_CLASS(FELogRigidBodyAngPosZ      , "thz");
REGISTER_FECORE_CLASS(FELogRigidBodyAngVelX      , "omx");
REGISTER_FECORE_CLASS(FELogRigidBodyAngVelY      , "omy");
REGISTER_FECORE_CLASS(FELogRigidBodyAngVelZ      , "omz");
REGISTER_FECORE_CLASS(FELogRigidBodyAngAccX      , "alx");
REGISTER_FECORE_CLASS(FELogRigidBodyAngAccY      , "aly");
REGISTER_FECORE_CLASS(FELogRigidBodyAngAccZ      , "alz");
REGISTER_FECORE_CLASS(FELogRigidBodyQuatX        , "qx");
REGISTER_FECORE_CLASS(FELogRigidBodyQuatY        , "qy");
REGISTER_FECORE_CLASS(FELogRigidBodyQuatZ        , "qz");
REGISTER_FECORE_CLASS(FELogRigidBodyQuatW        , "qw");
REGISTER_FECORE_CLASS(FELogRigidBodyR11          , "R11");
REGISTER_FECORE_CLASS(FELogRigidBodyR12          , "R12");
REGISTER_FECORE_CLASS(FELogRigidBodyR13          , "R13");
REGISTER_FECORE_CLASS(FELogRigidBodyR21          , "R21");
REGISTER_FECORE_CLASS(FELogRigidBodyR22          , "R22");
REGISTER_FECORE_CLASS(FELogRigidBodyR23          , "R23");
REGISTER_FECORE_CLASS(FELogRigidBodyR31          , "R31");
REGISTER_FECORE_CLASS(FELogRigidBodyR32          , "R32");
REGISTER_FECORE_CLASS(FELogRigidBodyR33          , "R33");
REGISTER_FECORE_CLASS(FELogRigidBodyForceX       , "Fx");
REGISTER_FECORE_CLASS(FELogRigidBodyForceY       , "Fy");
REGISTER_FECORE_CLASS(FELogRigidBodyForceZ       , "Fz");
REGISTER_FECORE_CLASS(FELogRigidBodyTorqueX      , "Mx");
REGISTER_FECORE_CLASS(FELogRigidBodyTorqueY      , "My");
REGISTER_FECORE_CLASS(FELogRigidBodyTorqueZ      , "Mz");
REGISTER_FECORE_CLASS(FELogRigidBodyKineticEnergy, "KE");

//-----------------------------------------------------------------------------
    // Derived from FELogConnectorData
REGISTER_FECORE_CLASS(FELogRigidConnectorForceX , "RCFx");
REGISTER_FECORE_CLASS(FELogRigidConnectorForceY , "RCFy");
REGISTER_FECORE_CLASS(FELogRigidConnectorForceZ , "RCFz");
REGISTER_FECORE_CLASS(FELogRigidConnectorMomentX, "RCMx");
REGISTER_FECORE_CLASS(FELogRigidConnectorMomentY, "RCMy");
REGISTER_FECORE_CLASS(FELogRigidConnectorMomentZ, "RCMz");

//-----------------------------------------------------------------------------
// Derived from FEMeshAdaptorCriterion
REGISTER_FECORE_CLASS(FEMaxStressCriterion, "max_stress");
REGISTER_FECORE_CLASS(FEMaxDamageCriterion, "max_damage");
}

//-----------------------------------------------------------------------------
// Derived from FELogNLConstraintData
REGISTER_FECORE_CLASS(FELogVolumeConstraint , "constrained volume");
REGISTER_FECORE_CLASS(FELogVolumePressure   , "volume pressure"   );

