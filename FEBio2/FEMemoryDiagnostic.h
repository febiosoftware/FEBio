#pragma once
#include "FEDiagnostic.h"

class FEMemoryDiagnostic :	public FEDiagnostic
{
public:
	FEMemoryDiagnostic(FEModel& fem);
	~FEMemoryDiagnostic();

	bool Init();

	bool Run();

	bool ParseSection(XMLTag& tag);

protected:
	char	m_szfile[512];
	int		m_iters;
};
