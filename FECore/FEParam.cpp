#include "stdafx.h"
#include "FEParam.h"
#include "FEParamValidator.h"

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

