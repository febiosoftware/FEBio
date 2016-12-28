#pragma once
#include "matrix.h"

class FESolidDomain;
class FESolidElement;

void IntegrateBDB(FESolidDomain& dom, FESolidElement& el, double D, matrix& ke);
