#pragma once
#include "matrix.h"

class FESolidDomain;
class FESolidElement;

//-----------------------------------------------------------------------------
// Integrator function for BDB forms
// where B is the shape function gradients
void IntegrateBDB(FESolidDomain& dom, FESolidElement& el, double D, matrix& ke);

//-----------------------------------------------------------------------------
// Integrator function for NCN forms
// where N are the shape functions
void IntegrateNCN(FESolidDomain& dom, FESolidElement& el, double C, matrix& ke);
