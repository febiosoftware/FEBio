#include "stdafx.h"
#include "FECoreKernel.h"
#include "Logfile.h"
#include <stdarg.h>
using namespace std;

// set the default linear solver (0 is equivalent to skyline solver)
int FECoreKernel::m_ndefault_solver = 0;

//-----------------------------------------------------------------------------
//! Helper function for reporting errors
bool fecore_error(const char* sz, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// make the message
	char szerr[512] = {0};
	va_start(args, sz);
	vsprintf(szerr, sz, args);
	va_end(args);

	// TODO: Perhaps I should report it to the logfile?
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	fecore.SetErrorString(szerr);

	return false;
}

//-----------------------------------------------------------------------------
const char* fecore_get_error_string()
{
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	return fecore.GetErrorString();
}

//-----------------------------------------------------------------------------
FECoreKernel* FECoreKernel::m_pKernel = 0;

//-----------------------------------------------------------------------------
FECoreKernel& FECoreKernel::GetInstance()
{
	if (m_pKernel == 0) m_pKernel = new FECoreKernel;
	return *m_pKernel;
}

//-----------------------------------------------------------------------------
// This function is used by plugins to make sure that the plugin and the executable
// are using the same kernel
void FECoreKernel::SetInstance(FECoreKernel* pkernel)
{
	m_pKernel = pkernel;
}

//-----------------------------------------------------------------------------
Logfile& FECoreKernel::GetLogfile()
{
	return *m_pKernel->m_plog;
}

//-----------------------------------------------------------------------------
FECoreKernel::FECoreKernel()
{
	m_plog = Logfile::GetInstance();
	m_szerr = 0;
	m_szmod = 0;
}

//-----------------------------------------------------------------------------
// Sets the error string.
// Calling SetErrorString(null) can be used to clear the error string.
void FECoreKernel::SetErrorString(const char* sz)
{
	// always clear the current error first
	if (m_szerr) delete [] m_szerr; m_szerr = 0;

	// make sure there is a new error string
	if (sz == 0) return;
	int l = strlen(sz);
	m_szerr = new char[l+1];
	strncpy(m_szerr, sz, l);
	m_szerr[l] = 0;
}

//-----------------------------------------------------------------------------
const char* FECoreKernel::GetErrorString()
{
	return m_szerr;
}

//-----------------------------------------------------------------------------
void FECoreKernel::RegisterClass(FECoreFactory* ptf)
{
	const char* szmod = m_szmod;
	if (szmod == 0) szmod = "";

	// see if the name already exists
	for (int i=0; i<m_Fac.size(); ++i)
	{
		FECoreFactory* pfi = m_Fac[i];

		if (pfi->GetSuperClassID() == ptf->GetSuperClassID())
		{
			const char* szmod_i = pfi->GetModuleName();
			if (szmod_i == 0) szmod_i = "";

			if ((strcmp(szmod, szmod_i) == 0) && (strcmp(pfi->GetTypeStr(), ptf->GetTypeStr()) == 0))
			{
#ifdef _DEBUG
				fprintf(stderr, "WARNING: %s feature is redefined\n", ptf->GetTypeStr());
#endif
				m_Fac[i] = ptf;
				return;
			}
		}
	}

	ptf->SetModuleName(m_szmod);
	m_Fac.push_back(ptf);
}

//-----------------------------------------------------------------------------
//! Create an object. An object is created by specifying the super-class id
//! and the type-string. 
void* FECoreKernel::Create(SUPER_CLASS_ID id, const char* sztype, FEModel* pfem)
{
	if (sztype == 0) return 0;

	const char* szmod = m_szmod;
	if (szmod == 0) szmod = "";

	// first find by module name
	std::vector<FECoreFactory*>::iterator pf;
	for (pf=m_Fac.begin(); pf!= m_Fac.end(); ++pf)
	{
		FECoreFactory* pfac = *pf;
		if (pfac->GetSuperClassID() == id) {

			const char* szmod_i = pfac->GetModuleName();
			if (szmod_i == 0) szmod_i = "";

			if ((strcmp(szmod_i, szmod) == 0) && (strcmp(pfac->GetTypeStr(), sztype) == 0))
			{
				return pfac->CreateInstance(pfem);
			}
		}
	}

	// we didn't find it.
	// Let's ignore module name
	// TODO: This is mostly for backward compatibility, but eventually should be removed
	for (pf = m_Fac.begin(); pf != m_Fac.end(); ++pf)
	{
		FECoreFactory* pfac = *pf;
		if (pfac->GetSuperClassID() == id) {
			if (strcmp(pfac->GetTypeStr(), sztype) == 0)
			{
				return pfac->CreateInstance(pfem);
			}
		}
	}

/*
#ifdef _DEBUG
	fprintf(stderr, "Unable to create class\n. These are the possible values:\n");
	for (pf=m_Fac.begin(); pf!=m_Fac.end(); ++pf)
	  {
	    FECoreFactory* pfac = *pf;
	    if (pfac->GetSuperClassID() == id) fprintf(stderr, "%s\n", pfac->GetTypeStr());
	  }
#endif
*/
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
int FECoreKernel::FactoryClasses()
{
	return (int) m_Fac.size();
}

//-----------------------------------------------------------------------------
const FECoreFactory* FECoreKernel::GetFactoryClass(int i)
{
	return m_Fac[i];
}

//-----------------------------------------------------------------------------
FECoreFactory* FECoreKernel::FindFactoryClass(int classID, const char* sztype)
{
	for (size_t i=0; i<m_Fac.size(); ++i)
	{
		FECoreFactory* fac = m_Fac[i];
		if ((fac->GetSuperClassID() == classID) &&
			(strcmp(fac->GetTypeStr(), sztype) == 0)) return fac;
	}
	return 0;
}

//-----------------------------------------------------------------------------
//! set the active module
void FECoreKernel::SetActiveModule(const char* szmod)
{
	m_szmod = szmod;
}

//-----------------------------------------------------------------------------
//! Get the active module
const char* FECoreKernel::GetActiveModule() const 
{ 
	return m_szmod; 
}

//-----------------------------------------------------------------------------
//! Register a new domain class
void FECoreKernel::RegisterDomain(FEDomainFactory* pf)
{
	m_Dom.push_back(pf); 
}

//-----------------------------------------------------------------------------
FEDomain* FECoreKernel::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	for (int i=0; i<(int)m_Dom.size(); ++i)
	{
		FEDomain* pdom = m_Dom[i]->CreateDomain(spec, pm, pmat);
		if (pdom != 0) return pdom;
	}
	return 0;
}

//-----------------------------------------------------------------------------
void FECoreKernel::RegisterLinearSolver(FELinearSolverFactory* pf)
{
	m_LS.push_back(pf); 
}

//-----------------------------------------------------------------------------
LinearSolver* FECoreKernel::CreateLinearSolver(int nsolver)
{
	for (int i=0; i<(int)m_LS.size(); ++i)
	{
		FELinearSolverFactory* pls = m_LS[i];
		if (pls->GetID() == nsolver) return pls->Create();
	}
	return 0;
}

//-----------------------------------------------------------------------------
FELinearSolverFactory* FECoreKernel::FindLinearSolverFactory(int nsolver)
{
	for (int i = 0; i<(int)m_LS.size(); ++i)
	{
		FELinearSolverFactory* pls = m_LS[i];
		if (pls->GetID() == nsolver) return pls;
	}
	return 0;
}
