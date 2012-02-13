#pragma once
#include "FEDiagnostic.h"

class FEPrintMatrixDiagnostic :	public FEDiagnostic
{
public:
	FEPrintMatrixDiagnostic(FEModel& fem);
	~FEPrintMatrixDiagnostic(void);

	bool ParseSection(XMLTag& tag);

	bool Run();

protected:
	char	m_szout[1024];
	int		m_rng[4];
};
