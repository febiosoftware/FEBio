#pragma once
#include "matrix.h"

//-----------------------------------------------------------------------------
// The purpose of this file is to explore mechanisms for evaluating the integrals 
// that commonly appear in FE formulations. The goal is that this would simplify 
// the implementation of new FE features.

//-----------------------------------------------------------------------------
class FESolidDomain;
class FESolidElement;
class FEMaterialPoint;

//-----------------------------------------------------------------------------
// This class can be used to evaluate quantities that depend on the material point.
// TODO: I wonder if it might be possible for a material to store these data classes
// directly. Maybe I can integrate this with the FEProperty class and create a
// FEMaterialProperty class that offers the functionality presented here. 
template <typename T>
class FEMaterialPointValue
{
public:
	FEMaterialPointValue(){}
	virtual ~FEMaterialPointValue(){}

	// overload this function and implement
	virtual T operator () (FEMaterialPoint& mp) = 0;
};

//-----------------------------------------------------------------------------
// Integrator function for BDB forms
// where B is the shape function gradients
void IntegrateBDB(FESolidDomain& dom, FESolidElement& el, double D, matrix& ke);
void IntegrateBDB(FESolidDomain& dom, FESolidElement& el, const mat3ds& D, matrix& ke);
void IntegrateBDB(FESolidDomain& dom, FESolidElement& el, FEMaterialPointValue<mat3ds>& D, matrix& ke);

//-----------------------------------------------------------------------------
// Integrator function for NCN forms
// where N are the shape functions
void IntegrateNCN(FESolidDomain& dom, FESolidElement& el, double C, matrix& ke);
