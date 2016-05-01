#include "stdafx.h"
#include "FEParam.h"
#include "FEParamValidator.h"
#include "DumpStream.h"
#include "FEFunction1D.h"

//-----------------------------------------------------------------------------
FEParam::FEParam(void* pdata, FEParamType itype, int ndim, const char* szname)
{
	m_pv = pdata;
	m_itype = itype;
	m_ndim = ndim;
	m_nlc = -1;
	m_scl = 1.0;
	m_vscl = vec3d(0, 0, 0);

	// set the name
	// note that we just copy the pointer, not the actual string
	// this is okay as long as the name strings are defined
	// as literal strings
	m_szname = szname;

	m_pvalid = 0;	// no default validator
}

//-----------------------------------------------------------------------------
FEParam::FEParam(const FEParam& p)
{
	m_pv = p.m_pv;
	m_itype = p.m_itype;
	m_ndim = p.m_ndim;
	m_nlc = p.m_nlc;
	m_scl = p.m_scl;
	m_vscl = p.m_vscl;
	m_szname = p.m_szname;

	m_pvalid = (p.m_pvalid ? p.m_pvalid->copy() : 0);
}

//-----------------------------------------------------------------------------
FEParam& FEParam::operator=(const FEParam& p)
{
	m_pv = p.m_pv;
	m_itype = p.m_itype;
	m_ndim = p.m_ndim;
	m_nlc = p.m_nlc;
	m_scl = p.m_scl;
	m_vscl = p.m_vscl;
	m_szname = p.m_szname;

	if (m_pvalid) delete m_pvalid;
	m_pvalid = (p.m_pvalid ? p.m_pvalid->copy() : 0);

	return *this;
}

//-----------------------------------------------------------------------------
bool FEParam::is_valid() const
{
	if (m_pvalid) return m_pvalid->is_valid(*this);
	return true;
}

//-----------------------------------------------------------------------------
//! This function deletes the existing validator and replaces it with the parameter
//! passed in the function. 
//! The pvalid can be null in which case the parameter will no longer be validated.
//! (i.e. is_valid() will always return true.
//! TODO: Should I delete the validator here? What if it was allocated in a plugin?
//!       Perhaps I should just return the old validator?
void FEParam::SetValidator(FEParamValidator* pvalid)
{
	if (m_pvalid) delete m_pvalid;
	m_pvalid = pvalid;
}

//-----------------------------------------------------------------------------
//! Sets the load curve ID and scale factor
void FEParam::SetLoadCurve(int lc)
{
	m_nlc = lc;
}

//-----------------------------------------------------------------------------
//! Sets the load curve ID and scale factor
void FEParam::SetLoadCurve(int lc, double s)
{
	assert(m_itype == FE_PARAM_DOUBLE);
	m_nlc = lc;
	m_scl = s;
}

//-----------------------------------------------------------------------------
//! Sets the load curve ID and scale factor
void FEParam::SetLoadCurve(int lc, const vec3d& v)
{
	assert(m_itype == FE_PARAM_VEC3D);
	m_nlc = lc;
	m_vscl = v;
}

//-----------------------------------------------------------------------------
void FEParam::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_nlc;
		ar << m_scl;
		ar << m_vscl;
		ar << (int) m_itype;
		ar << m_ndim;
		if (m_ndim == 1)
		{
			switch (m_itype)
			{
			case FE_PARAM_INT   : ar << value<int   >(); break;
			case FE_PARAM_BOOL  : ar << value<bool  >(); break;
			case FE_PARAM_DOUBLE: ar << value<double>(); break;
			case FE_PARAM_VEC3D : ar << value<vec3d >(); break;
			case FE_PARAM_MAT3D : ar << value<mat3d >(); break;
			case FE_PARAM_MAT3DS: ar << value<mat3ds>(); break;
			case FE_PARAM_STRING: ar << (const char*) data_ptr(); break;
			case FE_PARAM_FUNC1D: 
				{
					FEFunction1D& f = value<FEFunction1D>();
					f.Serialize(ar);
				}
				break;
			default:
				assert(false);
			}
		}
		else
		{
			switch (m_itype)
			{
			case FE_PARAM_INT:
				{
					int* pi = (int*) m_pv;
					for (int i=0; i<m_ndim; ++i) ar << pi[i];
				}
				break;
			case FE_PARAM_DOUBLE:
				{
					double* pv = (double*) m_pv;
					for (int i=0; i<m_ndim; ++i) ar << pv[i];
				}
				break;
			default:
				assert(false);
			}
		}	
	}
	else
	{
		ar >> m_nlc;
		ar >> m_scl;
		ar >> m_vscl;
		int ntype, ndim;
		ar >> ntype;
		ar >> ndim;
		if (ndim != m_ndim) throw DumpStream::ReadError();
		if (ntype != m_itype) throw DumpStream::ReadError();
		if (m_ndim == 1)
		{
			switch (m_itype)
			{
			case FE_PARAM_INT   : ar >> value<int   >(); break;
			case FE_PARAM_BOOL  : ar >> value<bool  >(); break;
			case FE_PARAM_DOUBLE: ar >> value<double>(); break;
			case FE_PARAM_VEC3D : ar >> value<vec3d >(); break;
			case FE_PARAM_MAT3D : ar >> value<mat3d >(); break;
			case FE_PARAM_MAT3DS: ar >> value<mat3ds>(); break;
			case FE_PARAM_STRING: ar >> (char*) data_ptr(); break;
			case FE_PARAM_FUNC1D: 
				{
					FEFunction1D& f = value<FEFunction1D>();
					f.Serialize(ar);
				}
				break;
			default:
				assert(false);
			}
		}
		else
		{
			switch (m_itype)
			{
			case FE_PARAM_INT:
				{
					int* pi = (int*) data_ptr();
					for (int i=0; i<m_ndim; ++i) ar >> pi[i];
				}
				break;
			case FE_PARAM_DOUBLE:
				{
					double* pv = (double*) data_ptr();
					for (int i=0; i<m_ndim; ++i) ar >> pv[i];
				}
				break;
			default:
				assert(false);
			}
		}
	}

	// serialize the validator
	if (m_pvalid) m_pvalid->Serialize(ar);
}
