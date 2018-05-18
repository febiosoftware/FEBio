#pragma once
#include "FEDiagnostic.h"
#include <FECore/FECoreTask.h>

class FEBioDiagnostic : public FECoreTask
{
public:
	FEBioDiagnostic(FEModel* pfem) : FECoreTask(pfem){ m_pdia = 0; }

	//! initialization
	bool Init(const char* szfile);

	//! Run the FE model
	virtual bool Run();

private:
	FEDiagnostic*	m_pdia;
};
