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
//! Create an object. An object is created by specifying the super-class id
//! and the type-string. 
void* FEBioKernel::Create(SUPER_CLASS_ID id, const char* sztype, FEModel* pfem)
{
	std::vector<FEBioFactory*>::iterator pf;
	for (pf=m_Fac.begin(); pf!= m_Fac.end(); ++pf)
	{
		FEBioFactory* pfac = *pf;
		if (pfac->GetSuperClassID() == id) {
			if (strcmp(pfac->GetTypeStr(), sztype) == 0) return pfac->CreateInstance(pfem);
		}
	}

#ifdef _DEBUG
	fprintf(stderr, "Unable to create class\n. These are the possible values:\n");
	for (pf=m_Fac.begin(); pf!=m_Fac.end(); ++pf)
	  {
	    FEBioFactory* pfac = *pf;
	    if (pfac->GetSuperClassID() == id) fprintf(stderr, "%s\n", pfac->GetTypeStr());
	  }
#endif

	return 0;
}

//-----------------------------------------------------------------------------
int FEBioKernel::Count(SUPER_CLASS_ID sid)
{
	int N = 0;
	std::vector<FEBioFactory*>::iterator pf;
	for (pf=m_Fac.begin(); pf!= m_Fac.end(); ++pf)
	{
		FEBioFactory* pfac = *pf;
		if (pfac->GetSuperClassID() == sid) N++;
	}
	return N;
}

//-----------------------------------------------------------------------------
void FEBioKernel::List(SUPER_CLASS_ID sid)
{
  std::vector<FEBioFactory*>::iterator pf;
  for (pf = m_Fac.begin(); pf != m_Fac.end(); ++pf)
    {
      FEBioFactory* pfac = *pf;
      if (pfac->GetSuperClassID() == sid) fprintf(stdout, "%s\n", pfac->GetTypeStr());
    }
}

//-----------------------------------------------------------------------------
FEBioKernel::FEBioKernel()
{
	m_plog = Logfile::GetInstance();
}

//-----------------------------------------------------------------------------
int FEBioKernel::GetDomainType(const FE_Element_Spec& spec, FEMaterial* pmat)
{
	for (int i=0; i<(int)m_Dom.size(); ++i)
	{
		int ndom = m_Dom[i]->GetDomainType(spec, pmat);
		if (ndom != 0) return ndom;
	}
	return 0;
}

//-----------------------------------------------------------------------------
FEDomain* FEBioKernel::CreateDomain(int dtype, FEMesh* pm, FEMaterial* pmat)
{
	for (int i=0; i<(int)m_Dom.size(); ++i)
	{
		FEDomain* pdom = m_Dom[i]->CreateDomain(dtype, pm, pmat);
		if (pdom != 0) return pdom;
	}
	return 0;
}
