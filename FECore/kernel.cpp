#include "stdafx.h"
#include "febio.h"
using namespace std;

FEBioKernel* FEBioKernel::m_pKernel = 0;

FEBioKernel& FEBioKernel::GetInstance()
{
	if (m_pKernel == 0) m_pKernel = new FEBioKernel;
	return *m_pKernel;
}
