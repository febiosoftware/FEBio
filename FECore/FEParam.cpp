#include "stdafx.h"
#include "FEParam.h"
#include "FEParamValidator.h"
#include "DumpStream.h"
#include "FEDataArray.h"
#include "tens3d.h"
#include "FEModelParam.h"

//-----------------------------------------------------------------------------
FEParam::FEParam(void* pdata, FEParamType itype, int ndim, const char* szname)
{
	m_pv = pdata;
	m_type = itype;
	m_dim = ndim;

	m_flag = 0;

	m_nlc = -1;
	m_scl = 1.0;
	m_vscl = vec3d(0, 0, 0);

	// set the name
	// note that we just copy the pointer, not the actual string
	// this is okay as long as the name strings are defined
	// as literal strings
	m_szname = szname;

	m_szenum = 0;

	m_pvalid = 0;	// no default validator

	m_parent = 0;
}

//-----------------------------------------------------------------------------
FEParam::FEParam(const FEParam& p)
{
	m_pv = p.m_pv;
	m_type = p.m_type;
	m_dim = p.m_dim;

	m_flag = p.m_flag;

	m_nlc = p.m_nlc;
	m_scl = p.m_scl;
	m_vscl = p.m_vscl;
	m_szname = p.m_szname;
	m_szenum = 0;
	m_parent = p.m_parent;

	m_pvalid = (p.m_pvalid ? p.m_pvalid->copy() : 0);
}

//-----------------------------------------------------------------------------
FEParam& FEParam::operator=(const FEParam& p)
{
	m_pv = p.m_pv;
	m_type = p.m_type;
	m_dim = p.m_dim;

	m_flag = p.m_flag;

	m_nlc = p.m_nlc;
	m_scl = p.m_scl;
	m_vscl = p.m_vscl;
	m_szname = p.m_szname;
	m_szenum = 0;
	m_parent = p.m_parent;

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
FEParamValue FEParam::paramValue(int i)
{
	if (m_dim == 1)
	{
		assert(i == -1);
		if (m_type == FE_PARAM_DOUBLE)
		{
			int lc = GetLoadCurve();
			if (lc == -1) return FEParamValue(this, &value<double>(), FE_PARAM_DOUBLE, -1);
			else return FEParamValue(this, &m_scl, FE_PARAM_DOUBLE, -1);
		}
		else if (m_type == FE_PARAM_DOUBLE_MAPPED)
		{
			FEParamDouble& p = value<FEParamDouble>();
			if (p.isConst()) return FEParamValue(this, &p.constValue(), FE_PARAM_DOUBLE, -1);
			else return FEParamValue(this, m_pv, m_type, -1);
		}
		else if (m_type == FE_PARAM_VEC3D_MAPPED)
		{
			FEParamVec3& p = value<FEParamVec3>();
			if (p.isConst()) return FEParamValue(this, &p.constValue(), FE_PARAM_VEC3D, -1);
			else return FEParamValue(this, m_pv, m_type, -1);
		}
		else if (m_type == FE_PARAM_MAT3D_MAPPED)
		{
			FEParamMat3d& p = value<FEParamMat3d>();
			if (p.isConst()) return FEParamValue(this, &p.constValue(), FE_PARAM_MAT3D, -1);
			else return FEParamValue(this, m_pv, m_type, -1);
		}
		else return FEParamValue(this, m_pv, m_type, -1);
	}
	else
	{
		switch (m_type)
		{
		case FE_PARAM_DOUBLE: return FEParamValue(this, &value<double>(i), FE_PARAM_DOUBLE, i); break;
		case FE_PARAM_VEC3D: return FEParamValue(this, &value<vec3d>(i), FE_PARAM_VEC3D, i); break;
		case FE_PARAM_STD_VECTOR_DOUBLE:
		{
			vector<double>& data = value< vector<double> >();
			if ((i >= 0) && (i < (int)data.size()))
				return FEParamValue(this, &data[i], FE_PARAM_DOUBLE, i);
		}
		break;
		case FE_PARAM_STD_VECTOR_VEC2D:
		{
			vector<vec2d>& data = value< vector<vec2d> >();
			if ((i >= 0) && (i < (int)data.size()))
				return FEParamValue(this, &data[i], FE_PARAM_VEC2D, i);
		}
		break;
		case FE_PARAM_STD_VECTOR_STRING:
		{
			vector<string>& data = value< vector<string> >();
			if ((i >= 0) && (i < (int)data.size()))
				return FEParamValue(this, &data[i], FE_PARAM_STD_STRING, i);
		}
		break;
		case FE_PARAM_DOUBLE_MAPPED:
		{
			FEParamDouble& data = value<FEParamDouble>(i);
			return FEParamValue(this, &data, FE_PARAM_DOUBLE_MAPPED, i);
		}
		break;
		case FE_PARAM_VEC3D_MAPPED:
		{
			FEParamVec3& data = value<FEParamVec3>(i);
			return FEParamValue(this, &data, FE_PARAM_VEC3D_MAPPED, i);
		}
		break;
		case FE_PARAM_MAT3D_MAPPED:
		{
			FEParamMat3d& data = value<FEParamMat3d>(i);
			return FEParamValue(this, &data, FE_PARAM_MAT3D_MAPPED, i);
		}
		break;
		}
	}

	return FEParamValue();
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
	assert(m_type == FE_PARAM_DOUBLE);
	m_nlc = lc;
	m_scl = s;
}

//-----------------------------------------------------------------------------
//! Sets the load curve ID and scale factor
void FEParam::SetLoadCurve(int lc, const vec3d& v)
{
	assert(m_type == FE_PARAM_VEC3D);
	m_nlc = lc;
	m_vscl = v;
}

//-----------------------------------------------------------------------------
void FEParam::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << (int) m_type;
		if (m_dim == 1)
		{
			switch (m_type)
			{
			case FE_PARAM_INT       : ar << value<int>(); break;
			case FE_PARAM_BOOL      : ar << value<bool>(); break;
			case FE_PARAM_DOUBLE    : ar << value<double>(); break;
			case FE_PARAM_VEC3D     : ar << value<vec3d>(); break;
			case FE_PARAM_MAT3D     : ar << value<mat3d>(); break;
			case FE_PARAM_MAT3DS    : ar << value<mat3ds>(); break;
			case FE_PARAM_TENS3DRS  : ar << value<tens3ds>(); break;
			case FE_PARAM_DATA_ARRAY:
			{
				FEDataArray& m = value<FEDataArray>();
				m.Serialize(ar);
			}
			break;
			case FE_PARAM_STRING: ar << (const char*)data_ptr(); break;
			case FE_PARAM_STD_STRING: ar << value<string>(); break;
			default:
				assert(false);
			}
		}
		else
		{
			switch (m_type)
			{
			case FE_PARAM_INT:
			{
				int* pi = (int*) m_pv;
				for (int i = 0; i<m_dim; ++i) ar << pi[i];
			}
			break;
			case FE_PARAM_DOUBLE:
			{
				double* pv = (double*) m_pv;
				for (int i = 0; i<m_dim; ++i) ar << pv[i];
			}
			break;
			default:
				assert(false);
			}
		}
	}
	else
	{
		int ntype, ndim;
		ar >> ntype;
		ar >> ndim;
		if (ndim != m_dim) throw DumpStream::ReadError();
		if (ntype != (int) m_type) throw DumpStream::ReadError();
		if (m_dim == 1)
		{
			switch (m_type)
			{
			case FE_PARAM_INT       : ar >> value<int         >(); break;
			case FE_PARAM_BOOL      : ar >> value<bool        >(); break;
			case FE_PARAM_DOUBLE    : ar >> value<double      >(); break;
			case FE_PARAM_VEC3D     : ar >> value<vec3d       >(); break;
			case FE_PARAM_MAT3D     : ar >> value<mat3d       >(); break;
			case FE_PARAM_MAT3DS    : ar >> value<mat3ds      >(); break;
			case FE_PARAM_TENS3DRS  : ar >> value<tens3drs>(); break;
			case FE_PARAM_DATA_ARRAY:
			{
				FEDataArray& m = value<FEDataArray>();
				m.Serialize(ar);
			}
			break;
			case FE_PARAM_STRING: ar >> (char*)data_ptr(); break;
			case FE_PARAM_STD_STRING: ar >> value<string>(); break;
			default:
				assert(false);
			}
		}
		else
		{
			switch (m_type)
			{
			case FE_PARAM_INT:
			{
				int* pi = (int*)data_ptr();
				for (int i = 0; i<m_dim; ++i) ar >> pi[i];
			}
			break;
			case FE_PARAM_DOUBLE:
			{
				double* pv = (double*)data_ptr();
				for (int i = 0; i<m_dim; ++i) ar >> pv[i];
			}
			break;
			default:
				assert(false);
			}
		}
	}
	// serialize the parameter 
	if (ar.IsSaving())
	{
		ar << m_nlc;
		ar << m_scl;
		ar << m_vscl;
	}
	else
	{
		ar >> m_nlc;
		ar >> m_scl;
		ar >> m_vscl;
	}

	// serialize the validator
	if (m_pvalid) m_pvalid->Serialize(ar);
}

//-----------------------------------------------------------------------------
//! This function copies the state of a parameter to this parameter.
//! This assumes that the parameters are compatible (i.e. have the same type)
//! This is used in FEParamContainer::CopyParameterListState()
bool FEParam::CopyState(const FEParam& p)
{
	if (p.type() != type()) return false;

	m_nlc = p.m_nlc;
	m_scl = p.m_scl;
	m_vscl = p.m_vscl;

	return true;
}

//-----------------------------------------------------------------------------
// helper functions for accessing components of parameters via parameter strings
FEParamValue GetParameterComponent(const ParamString& paramName, FEParam* param)
{
	// make sure we have something to do
	if (param == 0) return FEParamValue();

	if (param->type() == FE_PARAM_DOUBLE)
	{
		return param->paramValue(paramName.Index());
	}
	else if (param->type() == FE_PARAM_VEC3D)
	{
		vec3d* v = param->pvalue<vec3d>(0);
		assert(v);
		if (v)
		{
			if      (paramName == "x") return FEParamValue(param, &v->x, FE_PARAM_DOUBLE, 0);
			else if (paramName == "y") return FEParamValue(param, &v->y, FE_PARAM_DOUBLE, 1);
			else if (paramName == "z") return FEParamValue(param, &v->z, FE_PARAM_DOUBLE, 2);
			else return FEParamValue();
		}
		else return FEParamValue();
	}
	else if (param->type() == FE_PARAM_STD_VECTOR_DOUBLE)
	{
		int index = paramName.Index();
		return param->paramValue(index);
	}
	else if (param->type() == FE_PARAM_STD_VECTOR_VEC2D)
	{
		int index = paramName.Index();
		return param->paramValue(index);
	}
	else if (param->type() == FE_PARAM_STD_VECTOR_STRING)
	{
		int index = paramName.Index();
		return param->paramValue(index);
	}
	else if (param->type() == FE_PARAM_DOUBLE_MAPPED)
	{
		return param->paramValue(paramName.Index());
	}
	else if (param->type() == FE_PARAM_VEC3D_MAPPED)
	{
		return param->paramValue(paramName.Index());
	}
	else if (param->type() == FE_PARAM_MAT3D_MAPPED)
	{
		return param->paramValue(paramName.Index());
	}
	else if (param->type() == FE_PARAM_MATERIALPOINT)
	{
		return param->paramValue(paramName.Index());
	}

	return FEParamValue();
}
