#pragma once
#include <FECore/FECoreTask.h>

//-----------------------------------------------------------------------------
// This diagnostics tests the running and cold restart features
class FERestartDiagnostic : public FECoreTask
{
public:
	// constructor
	FERestartDiagnostic(FEModel* pfem);

	// initialize the diagnostic
	bool Init(const char* sz);

	// run the diagnostic
	bool Run();

public:
	bool	m_bok;
	char	m_szdmp[256];	// restart file name
};
