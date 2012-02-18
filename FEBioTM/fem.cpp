#include "stdafx.h"
#include "fem.h"
#include "FEBioLib/FEBioImport.h"

//-----------------------------------------------------------------------------
FEM::FEM()
{
}

//-----------------------------------------------------------------------------
bool FEM::Input(const char* szfile)
{
	// create file reader
	FEFEBioImport fim;

	// Load the file
	if (fim.Load(*this, szfile) == false)
	{
		char szerr[256];
		fim.GetErrorMessage(szerr);
		fprintf(stderr, szerr);

		return false;
	}

	// see if user redefined output filenames
	if (fim.m_szdmp[0]) SetDumpFilename(fim.m_szdmp);
	if (fim.m_szlog[0]) SetLogFilename (fim.m_szlog);
	if (fim.m_szplt[0]) SetPlotFilename(fim.m_szplt);

	// we're done reading
	return true;
}
