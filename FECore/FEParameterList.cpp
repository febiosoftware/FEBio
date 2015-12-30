#include "stdafx.h"
#include "FEParameterList.h"
#include <cstring>
#include <assert.h>

//-----------------------------------------------------------------------------
bool is_inside_range_int(int ival, FEParamRange rng, int imin, int imax)
{
	switch (rng)
	{
	case FE_GREATER         : return (ival >  imin); break;
	case FE_GREATER_OR_EQUAL: return (ival >= imin); break;
	case FE_LESS            : return (ival <  imin); break;
	case FE_LESS_OR_EQUAL   : return (ival <= imin); break;
	case FE_OPEN            : return ((ival >  imin) && (ival <  imax)); break;
	case FE_CLOSED          : return ((ival >= imin) && (ival <= imax)); break;
	case FE_LEFT_OPEN       : return ((ival >  imin) && (ival <= imax)); break;
	case FE_RIGHT_OPEN      : return ((ival >= imin) && (ival <  imax)); break;
	case FE_NOT_EQUAL       : return (ival != imin); break;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool is_inside_range_double(double val, FEParamRange rng, double dmin, double dmax)
{
	switch (rng)
	{
	case FE_GREATER         : return (val >  dmin); break;
	case FE_GREATER_OR_EQUAL: return (val >= dmin); break;
	case FE_LESS            : return (val <  dmin); break;
	case FE_LESS_OR_EQUAL   : return (val <= dmin); break;
	case FE_OPEN            : return ((val >  dmin) && (val <  dmax)); break;
	case FE_CLOSED          : return ((val >= dmin) && (val <= dmax)); break;
	case FE_LEFT_OPEN       : return ((val >  dmin) && (val <= dmax)); break;
	case FE_RIGHT_OPEN      : return ((val >= dmin) && (val <  dmax)); break;
	case FE_NOT_EQUAL       : return (val != dmin); break;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEParam::is_inside_range()
{
	if (m_irange == FE_DONT_CARE) return true;

	if (m_ndim == 1)
	{
		if (m_itype == FE_PARAM_INT)
		{
			int ival = value<int>();
			return is_inside_range_int(ival, m_irange, m_imin, m_imax);
		}
		else if (m_itype == FE_PARAM_DOUBLE)
		{
			double val = value<double>();
			return is_inside_range_double(val, m_irange, m_dmin, m_dmax);
		}
	}
	else
	{
		if (m_itype == FE_PARAM_INT)
		{
			for (int i=0; i<m_ndim; ++i)
			{
				int val = *(pvalue<int>(i));
				bool b = is_inside_range_int(val, m_irange, m_imin, m_imax);
				if (b == false) return false;
			}
			return true;
		}
		else if (m_itype == FE_PARAM_DOUBLE)
		{
			for (int i=0; i<m_ndim; ++i)
			{
				double val = *(pvalue<double>(i));
				bool b = is_inside_range_double(val, m_irange, m_dmin, m_dmax);;
				if (b == false) return false;
			}
			return true;
		}
	}

	// we can get here if the user specified a range type
	// for a parameter that is not an int or double.
	return false;
}

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
		if (s.m_ndim == 1)
		{
			switch (s.m_itype)
			{
			case FE_PARAM_INT   : d.value<int   >() = s.value<int   >(); break;
			case FE_PARAM_BOOL  : d.value<bool  >() = s.value<bool  >(); break;
			case FE_PARAM_DOUBLE: d.value<double>() = s.value<double>(); break;
			case FE_PARAM_VEC3D : d.value<vec3d >() = s.value<vec3d >(); break;
			case FE_PARAM_MAT3D : d.value<mat3d >() = s.value<mat3d >(); break;
			case FE_PARAM_MAT3DS: d.value<mat3ds>() = s.value<mat3ds>(); break;
			default:
				assert(false);
			}
		}
		else
		{
			switch (s.m_itype)
			{
			case FE_PARAM_INT:
				{
					for (int i=0; i<s.m_ndim; ++i) d.pvalue<int>()[i] = s.pvalue<int>()[i];
				}
				break;
			case FE_PARAM_DOUBLE:
				{
					for (int i=0; i<s.m_ndim; ++i) d.pvalue<double>()[i] = s.pvalue<double>()[i];
				}
				break;
			default:
				assert(false);
			}
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
// This function adds a parameter to the parameter list
void FEParameterList::AddParameter(void *pv, FEParamType itype, int ndim, FEParamRange rng, double fmin, double fmax, const char *sz)
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

	// set the range
	// (range checking is only supported for int and double params)
	p.m_irange = rng;
	if (itype == FE_PARAM_INT)
	{
		p.m_imax = (int) fmax;
		p.m_imin = (int) fmin;
	}
	else if (itype == FE_PARAM_DOUBLE)
	{
		p.m_dmax = fmax;
		p.m_dmin = fmin;
	}
	else p.m_irange = FE_DONT_CARE;

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
// Add a parameter to the parameter list
void FEParamContainer::AddParameter(void* pv, FEParamType itype, int ndim, RANGE rng, const char* sz)
{
	assert(m_pParam);
	m_pParam->AddParameter(pv, itype, ndim, rng.m_rt, rng.m_fmin, rng.m_fmax, sz);
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
			ar << p.m_ndim;
			if (p.m_ndim == 1)
			{
				switch (p.m_itype)
				{
				case FE_PARAM_INT   : ar << p.value<int   >(); break;
				case FE_PARAM_BOOL  : ar << p.value<bool  >(); break;
				case FE_PARAM_DOUBLE: ar << p.value<double>(); break;
				case FE_PARAM_VEC3D : ar << p.value<vec3d >(); break;
				case FE_PARAM_MAT3D : ar << p.value<mat3d >(); break;
				case FE_PARAM_MAT3DS: ar << p.value<mat3ds>(); break;
				case FE_PARAM_STRING: ar << (const char*) p.m_pv; break;
				default:
					assert(false);
				}
			}
			else
			{
				switch (p.m_itype)
				{
				case FE_PARAM_INT:
					{
						int* pi = (int*) p.m_pv;
						for (int i=0; i<p.m_ndim; ++i) ar << pi[i];
					}
					break;
				case FE_PARAM_DOUBLE:
					{
						double* pv = (double*) p.m_pv;
						for (int i=0; i<p.m_ndim; ++i) ar << pv[i];
					}
					break;
				default:
					assert(false);
				}
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
				int ntype, ndim;
				ar >> ntype;
				ar >> ndim;
				if (ndim != p.m_ndim) throw DumpFile::ReadError();
				if (ntype != p.m_itype) throw DumpFile::ReadError();
				if (p.m_ndim == 1)
				{
					switch (p.m_itype)
					{
					case FE_PARAM_INT   : ar >> p.value<int   >(); break;
					case FE_PARAM_BOOL  : ar >> p.value<bool  >(); break;
					case FE_PARAM_DOUBLE: ar >> p.value<double>(); break;
					case FE_PARAM_VEC3D : ar >> p.value<vec3d >(); break;
					case FE_PARAM_MAT3D : ar >> p.value<mat3d >(); break;
					case FE_PARAM_MAT3DS: ar >> p.value<mat3ds>(); break;
					case FE_PARAM_STRING: ar >> (char*) p.m_pv; break;
					default:
						assert(false);
					}
				}
				else
				{
					switch (p.m_itype)
					{
					case FE_PARAM_INT:
						{
							int* pi = (int*) p.m_pv;
							for (int i=0; i<p.m_ndim; ++i) ar >> pi[i];
						}
						break;
					case FE_PARAM_DOUBLE:
						{
							double* pv = (double*) p.m_pv;
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
}
