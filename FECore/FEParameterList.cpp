#include "stdafx.h"
#include "FEParameterList.h"
#include <cstring>
#include <assert.h>

//-----------------------------------------------------------------------------
// This function adds a parameter to the parameter list
void FEParameterList::AddParameter(void *pv, FEParamType itype, int ndim, const char *sz)
{
	// create a new parameter object
	FEParam p;

	// set the pointer to the value
	assert(pv);
	p.m_pv = pv;

	// set the type
	p.m_itype = itype;

	// set the dimension
	// make sure that the type is a vector type if the dimension > 1
	assert((ndim >= 1) && ( (ndim > 1 ? (itype >= 100) : true)));
	p.m_ndim = ndim;

	// set the name
	// note that we just copy the pointer, not the actual string
	// this is okay as long as the name strings are defined
	// as literal strings
	assert(sz);
	p.m_szname = sz;

	// add the parameter to the list
	m_pl.push_back(p);
}

//-----------------------------------------------------------------------------
// This function searches the parameters in the list for a parameter
// with the name given by the input argument
// \param sz name of parameter to find
FEParam* FEParameterList::Find(const char* sz)
{
	FEParam* pp = 0;
	if (m_pl.size() > 0)
	{
		list<FEParam>::iterator it;
		for (it = m_pl.begin(); it != m_pl.end(); ++it)
		{
			if (strcmp(it->m_szname, sz) == 0)
			{
				pp = &(*it);
				break;
			}
		}
	}

	return pp;
}

//-----------------------------------------------------------------------------
FEParamContainer::FEParamContainer()
{
	m_pParam = 0;
}

//-----------------------------------------------------------------------------
FEParamContainer::~FEParamContainer()
{
	delete m_pParam; 
	m_pParam = 0;
}

//-----------------------------------------------------------------------------
FEParameterList& FEParamContainer::GetParameterList()
{
	if (m_pParam == 0) 
	{
		m_pParam = new FEParameterList;
		BuildParamList();
	}
	return *m_pParam;
}

//-----------------------------------------------------------------------------
// Find a parameter from its name
FEParam* FEParamContainer::GetParameter(const char* sz)
{
	FEParameterList& pl = GetParameterList();
	return pl.Find(sz);
}

//-----------------------------------------------------------------------------
// Add a parameter to the parameter list
void FEParamContainer::AddParameter(void* pv, FEParamType itype, int ndim, const char* sz)
{
	assert(m_pParam);
	m_pParam->AddParameter(pv, itype, ndim, sz);
}

//-----------------------------------------------------------------------------
// Serialize parameters to archive
void FEParamContainer::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		int NP = m_pParam->Parameters();
		ar << NP;
		list<FEParam>::iterator it = m_pParam->first();
		for (int i=0; i<NP; ++i)
		{
			FEParam& p = *it++;
			ar << p.m_nlc;
			ar << p.m_scl;
			ar << (int) p.m_itype;
			switch (p.m_itype)
			{
			case FE_PARAM_INT   : ar << p.value<int   >(); break;
			case FE_PARAM_BOOL  : ar << p.value<bool  >(); break;
			case FE_PARAM_DOUBLE: ar << p.value<double>(); break;
			case FE_PARAM_VEC3D : ar << p.value<vec3d >(); break;
			case FE_PARAM_STRING: ar << (const char*) p.m_pv; break;
			case FE_PARAM_INTV:
				{
					int* pi = (int*) p.m_pv;
					ar << p.m_ndim;
					for (int i=0; i<p.m_ndim; ++i) ar << pi[i];
				}
				break;
			case FE_PARAM_DOUBLEV:
				{
					double* pv = (double*) p.m_pv;
					ar << p.m_ndim;
					for (int i=0; i<p.m_ndim; ++i) ar << pv[i];
				}
				break;
			default:
				assert(false);
			}
		}
	}
	else
	{
		FEParameterList& pl = GetParameterList();
		int NP = 0;
		ar >> NP;
		if (NP != pl.Parameters()) throw DumpFile::ReadError();
		list<FEParam>::iterator it = pl.first();
		for (int i=0; i<NP; ++i)
		{
			FEParam& p = *it++;
			ar >> p.m_nlc;
			ar >> p.m_scl;
			int ntype;
			ar >> ntype;
			if (ntype != p.m_itype) throw DumpFile::ReadError();
			switch (p.m_itype)
			{
			case FE_PARAM_INT   : ar >> p.value<int   >(); break;
			case FE_PARAM_BOOL  : ar >> p.value<bool  >(); break;
			case FE_PARAM_DOUBLE: ar >> p.value<double>(); break;
			case FE_PARAM_VEC3D : ar >> p.value<vec3d >(); break;
			case FE_PARAM_STRING: ar >> (char*) p.m_pv; break;
			case FE_PARAM_INTV:
				{
					int* pi = (int*) p.m_pv;
					int ndim;
					ar >> ndim;
					if (ndim != p.m_ndim) throw DumpFile::ReadError();
					for (int i=0; i<p.m_ndim; ++i) ar >> pi[i];
				}
				break;
			case FE_PARAM_DOUBLEV:
				{
					double* pv = (double*) p.m_pv;
					int ndim;
					ar >> ndim;
					if (ndim != p.m_ndim) throw DumpFile::ReadError();
					for (int i=0; i<p.m_ndim; ++i) ar >> pv[i];
				}
				break;
			default:
				assert(false);
			}
		}
	}
}
