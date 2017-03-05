#pragma once
#include <FECore/FECoreTask.h>
#include <FECore/DumpMemStream.h>

//-----------------------------------------------------------------------------
// This diagnostics tests the running and cold restart features
class FECORE_EXPORT FERestartDiagnostic : public FECoreTask
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
	bool	m_bfile;		// file or memory stream?
	char	m_szdmp[256];	// restart file name
	DumpMemStream	m_dmp;
};
