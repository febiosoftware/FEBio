#pragma once
#include "FEDiagnostic.h"

class FEPrintHBMatrixDiagnostic :	public FEDiagnostic
{
public:
	FEPrintHBMatrixDiagnostic(FEM& fem);
	~FEPrintHBMatrixDiagnostic(void);

	bool ParseSection(XMLTag& tag);

	bool Run();

};
