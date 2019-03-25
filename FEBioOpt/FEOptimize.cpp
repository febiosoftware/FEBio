#include "stdafx.h"
#include "FEOptimize.h"
#include "FEOptimizeData.h"
#include "FECore/FECoreKernel.h"
#include "FECore/log.h"

//-----------------------------------------------------------------------------
#define VERSION 2
#define SUB_VERSION	0

//-----------------------------------------------------------------------------
FEOptimize::FEOptimize(FEModel* pfem) : FECoreTask(pfem), m_opt(pfem)
{

}

//-----------------------------------------------------------------------------
bool FEOptimize::Init(const char* szfile)
{
	char szversion[32] = { 0 };
	sprintf(szversion, "version %d.%d", VERSION, SUB_VERSION);
	feLog("P A R A M E T E R   O P T I M I Z A T I O N   M O D U L E\n%s\n\n", szversion);

	// read the data from the xml input file
	if (m_opt.Input(szfile) == false) return false;

	// do initialization
	bool ret = m_opt.Init();

	if (ret == false)
	{
		feLogErrorEx(m_opt.GetFEModel(), "Failed to initialize the optimization data.");
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEOptimize::Run()
{
	// solve the problem
	bool bret = m_opt.Solve();

	if (bret)
		feLog("\n\n N O R M A L   T E R M I N A T I O N\n\n");
	else 
		feLog("\n\n E R R O R   T E R M I N A T I O N\n\n");

	return bret;
}
