// Harwell-Boeing matrix print diagnostic class
//
///////////////////////////////////////////////

#pragma once
#include "FEDiagnostic.h"

//! Harwell-Boeing Matrix Print Diagnostic

//! Class to run a diagnostic to print the initial matrix in
//! Harwell-Boeing matrix format

class FEPrintHBMatrixDiagnostic :	public FEDiagnostic
{
public:
	FEPrintHBMatrixDiagnostic(FEM& fem);
	~FEPrintHBMatrixDiagnostic(void);

	bool ParseSection(XMLTag& tag);

	bool Run();

};
