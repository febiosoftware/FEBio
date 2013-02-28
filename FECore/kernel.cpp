#include "stdafx.h"
#include "febio.h"
#include "Logfile.h"
using namespace std;

//-----------------------------------------------------------------------------
FEBioKernel* FEBioKernel::m_pKernel = 0;

//-----------------------------------------------------------------------------
FEBioKernel& FEBioKernel::GetInstance()
{
	if (m_pKernel == 0) m_pKernel = new FEBioKernel;
	return *m_pKernel;
}

//-----------------------------------------------------------------------------
Logfile& FEBioKernel::GetLogfile()
{
	return *m_plog;
}

//-----------------------------------------------------------------------------
FEBioKernel::FEBioKernel()
{
	m_plog = Logfile::GetInstance();
}
