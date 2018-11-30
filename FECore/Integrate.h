#pragma once
#include "matrix.h"
#include <functional>

//-----------------------------------------------------------------------------
// The purpose of this file is to explore mechanisms for evaluating the integrals 
// that commonly appear in FE formulations. The goal is that this would simplify 
// the implementation of new FE features.

//-----------------------------------------------------------------------------
class FESolidDomain;
class FESolidElement;
class FEMaterialPoint;

//-----------------------------------------------------------------------------
// Integrator function for BDB forms
// where B is the shape function gradients
FECORE_API void IntegrateBDB(FESolidDomain& dom, FESolidElement& el, double D, matrix& ke);
FECORE_API void IntegrateBDB(FESolidDomain& dom, FESolidElement& el, const mat3ds& D, matrix& ke);
FECORE_API void IntegrateBDB(FESolidDomain& dom, FESolidElement& el, std::function<mat3ds (const FEMaterialPoint& mp)> f, matrix& ke);

//-----------------------------------------------------------------------------
// Integrator function for NCN forms
// where N are the shape functions
FECORE_API void IntegrateNCN(FESolidDomain& dom, FESolidElement& el, double C, matrix& ke);
