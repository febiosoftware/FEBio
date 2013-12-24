#include "stdafx.h"
#include "FECoreKernel.h"
#include "Logfile.h"
using namespace std;

//-----------------------------------------------------------------------------
FECoreKernel* FECoreKernel::m_pKernel = 0;

//-----------------------------------------------------------------------------
FECoreKernel& FECoreKernel::GetInstance()
{
	if (m_pKernel == 0) m_pKernel = new FECoreKernel;
	return *m_pKernel;
}

//-----------------------------------------------------------------------------
Logfile& FECoreKernel::GetLogfile()
{
	return *m_plog;
}

//-----------------------------------------------------------------------------
//! Create an object. An object is created by specifying the super-class id
//! and the type-string. 
void* FECoreKernel::Create(SUPER_CLASS_ID id, const char* sztype, FEModel* pfem)
{
	std::vector<FECoreFactory*>::iterator pf;
	for (pf=m_Fac.begin(); pf!= m_Fac.end(); ++pf)
	{
		FECoreFactory* pfac = *pf;
		if (pfac->GetSuperClassID() == id) {
			if (strcmp(pfac->GetTypeStr(), sztype) == 0) return pfac->CreateInstance(pfem);
		}
	}

#ifdef _DEBUG
	fprintf(stderr, "Unable to create class\n. These are the possible values:\n");
	for (pf=m_Fac.begin(); pf!=m_Fac.end(); ++pf)
	  {
	    FECoreFactory* pfac = *pf;
	    if (pfac->GetSuperClassID() == id) fprintf(stderr, "%s\n", pfac->GetTypeStr());
	  }
#endif

	return 0;
}

//-----------------------------------------------------------------------------
int FECoreKernel::Count(SUPER_CLASS_ID sid)
{
	int N = 0;
	std::vector<FECoreFactory*>::iterator pf;
	for (pf=m_Fac.begin(); pf!= m_Fac.end(); ++pf)
	{
		FECoreFactory* pfac = *pf;
		if (pfac->GetSuperClassID() == sid) N++;
	}
	return N;
}

//-----------------------------------------------------------------------------
void FECoreKernel::List(SUPER_CLASS_ID sid)
{
  std::vector<FECoreFactory*>::iterator pf;
  for (pf = m_Fac.begin(); pf != m_Fac.end(); ++pf)
    {
      FECoreFactory* pfac = *pf;
      if (pfac->GetSuperClassID() == sid) fprintf(stdout, "%s\n", pfac->GetTypeStr());
    }
}

//-----------------------------------------------------------------------------
FECoreKernel::FECoreKernel()
{
	m_plog = Logfile::GetInstance();
}

//-----------------------------------------------------------------------------
int FECoreKernel::GetDomainType(const FE_Element_Spec& spec, FEMaterial* pmat)
{
	for (int i=0; i<(int)m_Dom.size(); ++i)
	{
		int ndom = m_Dom[i]->GetDomainType(spec, pmat);
		if (ndom != 0) return ndom;
	}
	return 0;
}

//-----------------------------------------------------------------------------
FEDomain* FECoreKernel::CreateDomain(int dtype, FEMesh* pm, FEMaterial* pmat)
{
	for (int i=0; i<(int)m_Dom.size(); ++i)
	{
		FEDomain* pdom = m_Dom[i]->CreateDomain(dtype, pm, pmat);
		if (pdom != 0) return pdom;
	}
	return 0;
}
