#include "stdafx.h"
#include "FEParameterList.h"
#include <cstring>
#include <assert.h>

//-----------------------------------------------------------------------------
//! This function copies the parameter data from the passed parameter list.
//! This assumes that the two parameter lists are identical.
void FEParameterList::operator = (FEParameterList& l)
{
	if (m_pl.size() != l.m_pl.size()) { assert(false); return; }

	list<FEParam>::iterator ps, pd;
	ps = l.m_pl.begin();
	pd = m_pl.begin();
	for (; pd != m_pl.end(); ++pd, ++ps)
	{
		FEParam& s = *ps;
		FEParam& d = *pd;

		if (s.m_itype != d.m_itype) { assert(false); return; }
		if (s.m_ndim != d.m_ndim) { assert(false); return; }
		switch (s.m_itype)
		{
		case FE_PARAM_INT   : d.value<int   >() = s.value<int   >(); break;
		case FE_PARAM_BOOL  : d.value<bool  >() = s.value<bool  >(); break;
		case FE_PARAM_DOUBLE: d.value<double>() = s.value<double>(); break;
		case FE_PARAM_VEC3D : d.value<vec3d >() = s.value<vec3d >(); break;
		case FE_PARAM_MAT3D : d.value<mat3d >() = s.value<mat3d >(); break;
		case FE_PARAM_MAT3DS: d.value<mat3ds>() = s.value<mat3ds>(); break;
		case FE_PARAM_INTV  :
			{
				for (int i=0; i<s.m_ndim; ++i) d.pvalue<int>()[i] = s.pvalue<int>()[i];
			}
			break;
		case FE_PARAM_DOUBLEV:
			{
				for (int i=0; i<s.m_ndim; ++i) d.pvalue<double>()[i] = s.pvalue<double>()[i];
			}
			break;
		default:
			assert(false);
		}
	}
}

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

//=============================================================================

void cstrncpy(char* pd, char* ps, int l)
{
	for (int i=0; i<=l; ++i) pd[i] = ps[i];
}

ParamString::ParamString(const char* sz)
{
	m_sz = 0;
	m_nc = 0;
	m_nl = strlen(sz);
	m_sz = new char[m_nl+1];
	strcpy(m_sz, sz);
	char* ch = m_sz;
	while (ch = strchr(ch, '.')) { m_nc++; *ch = 0; ch = ch+1; }
	m_nc++;
}

//-----------------------------------------------------------------------------
ParamString::ParamString(const ParamString& p)
{
	m_nc = p.m_nc;
	m_nl = p.m_nl;
	m_sz = new char[m_nl+1];
	cstrncpy(m_sz, p.m_sz, m_nl);
}

//-----------------------------------------------------------------------------
void ParamString::operator=(const ParamString& p)
{
	m_nc = p.m_nc;
	m_nl = p.m_nl;
	delete [] m_sz;
	m_sz = new char[m_nl+1];
	cstrncpy(m_sz, p.m_sz, m_nl);
}

//-----------------------------------------------------------------------------
ParamString::~ParamString()
{
	delete [] m_sz;
	m_sz = 0;
	m_nc = 0;
	m_nl = 0;
}

//-----------------------------------------------------------------------------
int ParamString::count() const
{
	return m_nc;
}

//-----------------------------------------------------------------------------
ParamString ParamString::next() const
{
	int nc = m_nc - 1;
	int l = strlen(m_sz);
	char* sz = strchr(m_sz, '\0');
	char* sznew = new char[m_nl-l+1];
	cstrncpy(sznew, sz+1, m_nl-l);
	return ParamString(sznew, nc, m_nl-l);
}

//-----------------------------------------------------------------------------
bool ParamString::operator==(const char* sz) const
{
	return strcmp(m_sz, sz) == 0;
}

//=============================================================================
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
		m_pParam = new FEParameterList(this);
		BuildParamList();
	}
	return *m_pParam;
}

//-----------------------------------------------------------------------------
// Find a parameter from its name
FEParam* FEParamContainer::GetParameter(const ParamString& s)
{
	FEParameterList& pl = GetParameterList();
	return pl.Find(s.c_str());
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
		int NP = 0;
		list<FEParam>::iterator it;
		if (m_pParam) // If the input file doesn't set any parameters, the parameter list won't be created.
		{
			NP = m_pParam->Parameters();
			it = m_pParam->first();
		}
		ar << NP;
		for (int i=0; i<NP; ++i)
		{
			FEParam& p = *it++;
			ar << p.m_nlc;
			ar << p.m_scl;
			ar << p.m_vscl;
			ar << (int) p.m_itype;
			switch (p.m_itype)
			{
			case FE_PARAM_INT   : ar << p.value<int   >(); break;
			case FE_PARAM_BOOL  : ar << p.value<bool  >(); break;
			case FE_PARAM_DOUBLE: ar << p.value<double>(); break;
			case FE_PARAM_VEC3D : ar << p.value<vec3d >(); break;
			case FE_PARAM_MAT3D : ar << p.value<mat3d >(); break;
			case FE_PARAM_MAT3DS: ar << p.value<mat3ds>(); break;
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
		int NP = 0;
		ar >> NP;
		if (NP)
		{
			FEParameterList& pl = GetParameterList();
			if (NP != pl.Parameters()) throw DumpFile::ReadError();
			list<FEParam>::iterator it = pl.first();
			for (int i=0; i<NP; ++i)
			{
				FEParam& p = *it++;
				ar >> p.m_nlc;
				ar >> p.m_scl;
				ar >> p.m_vscl;
				int ntype;
				ar >> ntype;
				if (ntype != p.m_itype) throw DumpFile::ReadError();
				switch (p.m_itype)
				{
				case FE_PARAM_INT   : ar >> p.value<int   >(); break;
				case FE_PARAM_BOOL  : ar >> p.value<bool  >(); break;
				case FE_PARAM_DOUBLE: ar >> p.value<double>(); break;
				case FE_PARAM_VEC3D : ar >> p.value<vec3d >(); break;
				case FE_PARAM_MAT3D : ar >> p.value<mat3d >(); break;
				case FE_PARAM_MAT3DS: ar >> p.value<mat3ds>(); break;
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
}
