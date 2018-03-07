#include "stdafx.h"
#include "FEOptimize.h"
#include "FEOptimizeData.h"
#include "FECore/FECoreKernel.h"
#include "FECore/log.h"

//-----------------------------------------------------------------------------
#define VERSION 2
#define SUB_VERSION	0

//-----------------------------------------------------------------------------
FEOptimize::FEOptimize(FEModel* pfem) : FECoreTask(pfem), m_opt(*m_pfem)
{

}

//-----------------------------------------------------------------------------
bool FEOptimize::Init(const char* szfile)
{
	// read the data from the xml input file
	if (m_opt.Input(szfile) == false) return false;

	// do initialization
	felog.SetMode(Logfile::LOG_NEVER);
	if (m_opt.Init() == false) return false;

	felog.SetMode(Logfile::LOG_FILE_AND_SCREEN);
	char szversion[32] = {0};
	sprintf(szversion, "version %d.%d", VERSION, SUB_VERSION);
	felog.printbox("P A R A M E T E R   O P T I M I Z A T I O N   M O D U L E", szversion);

	return true;
}

//-----------------------------------------------------------------------------
bool FEOptimize::Run()
{
	// solve the problem
	bool bret = m_opt.Solve();

	if (bret)
		felog.printf("\n\n N O R M A L   T E R M I N A T I O N\n\n");
	else 
		felog.printf("\n\n E R R O R   T E R M I N A T I O N\n\n");

	return bret;
}
