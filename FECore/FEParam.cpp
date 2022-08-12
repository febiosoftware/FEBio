/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FEParam.h"
#include "FEParamValidator.h"
#include "DumpStream.h"
#include "FEDataArray.h"
#include "tens3d.h"
#include "FEModelParam.h"
using namespace std;

FEParamValue FEParamValue::component(int n)
{
	switch (m_itype)
	{
	case FE_PARAM_VEC2D:
	{
		vec2d& r = value<vec2d>();
		if (n == 0) return FEParamValue(m_param, &r.x(), FE_PARAM_DOUBLE);
		if (n == 1) return FEParamValue(m_param, &r.y(), FE_PARAM_DOUBLE);
	}
	break;
	case FE_PARAM_VEC3D:
	{
		vec3d& r = value<vec3d>();
		if (n == 0) return FEParamValue(m_param, &r.x, FE_PARAM_DOUBLE);
		if (n == 1) return FEParamValue(m_param, &r.y, FE_PARAM_DOUBLE);
		if (n == 2) return FEParamValue(m_param, &r.z, FE_PARAM_DOUBLE);
	}
	break;
	case FE_PARAM_MAT3DS:
	{
		mat3ds& a = value<mat3ds>();
		if (n == 0) return FEParamValue(m_param, &a.xx(), FE_PARAM_DOUBLE);
		if (n == 1) return FEParamValue(m_param, &a.xy(), FE_PARAM_DOUBLE);
		if (n == 2) return FEParamValue(m_param, &a.yy(), FE_PARAM_DOUBLE);
		if (n == 3) return FEParamValue(m_param, &a.xz(), FE_PARAM_DOUBLE);
		if (n == 4) return FEParamValue(m_param, &a.yz(), FE_PARAM_DOUBLE);
		if (n == 5) return FEParamValue(m_param, &a.zz(), FE_PARAM_DOUBLE);
	}
	break;
	}

	assert(false);
	return FEParamValue();
}


//-----------------------------------------------------------------------------
FEParam::FEParam(void* pdata, FEParamType itype, int ndim, const char* szname, bool* watch)
{
	m_pv = pdata;
	m_type = itype;
	m_dim = ndim;
	m_watch = watch;
	if (m_watch) *m_watch = false;

	// default flags depend on type
	// (see also FEModel::EvaluateLoadParameters())
	m_flag = 0;
	if (ndim == 1)
	{
		switch (itype)
		{
		case FE_PARAM_DOUBLE:
		case FE_PARAM_VEC3D:
		case FE_PARAM_DOUBLE_MAPPED:
		case FE_PARAM_VEC3D_MAPPED:
			// all these types can be modified via a load curve
			m_flag = FE_PARAM_VOLATILE;
			break;
		}
	}

	m_group = -1;

	// set the name
	// note that we just copy the pointer, not the actual string
	// this is okay as long as the name strings are defined
	// as literal strings
	m_szname = szname;
	m_szlongname = szname;

	m_szenum = 0;

	m_pvalid = 0;	// no default validator

	m_parent = 0;

	m_szunit = nullptr;
}

//-----------------------------------------------------------------------------
FEParam::FEParam(const FEParam& p)
{
	m_pv = p.m_pv;
	m_type = p.m_type;
	m_dim = p.m_dim;
	m_watch = p.m_watch;

	m_flag = p.m_flag;
	m_group = p.m_group;

	m_szname = p.m_szname;
	m_szlongname = p.m_szlongname;

	m_szenum = 0;
	m_parent = p.m_parent;

	m_szunit = p.m_szunit;

	m_pvalid = (p.m_pvalid ? p.m_pvalid->copy() : 0);
}

//-----------------------------------------------------------------------------
int FEParam::GetParamGroup() const
{
	return m_group;
}

//-----------------------------------------------------------------------------
void FEParam::SetParamGroup(int i)
{
	m_group = i;
}

//-----------------------------------------------------------------------------
FEParam::~FEParam()
{
	if (m_flag & FEParamFlag::FE_PARAM_USER)
	{
		free((void*)m_szname);
		assert(m_type == FE_PARAM_DOUBLE);
		delete (double*)m_pv;
	}
}

//-----------------------------------------------------------------------------
FEParam& FEParam::operator=(const FEParam& p)
{
	m_pv = p.m_pv;
	m_type = p.m_type;
	m_dim = p.m_dim;
	m_watch = p.m_watch;

	m_flag = p.m_flag;

	m_szname = p.m_szname;
	m_szenum = 0;
	m_parent = p.m_parent;

	m_szunit = p.m_szunit;

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
// return the name of the parameter
const char* FEParam::name() const 
{ 
	return m_szname; 
}

//-----------------------------------------------------------------------------
const char* FEParam::longName() const
{
	return m_szlongname;
}

//-----------------------------------------------------------------------------
// return the enum values
const char* FEParam::enums() const 
{ 
	return m_szenum; 
}

//-----------------------------------------------------------------------------
// get the current enum value (or nullptr)
const char* FEParam::enumKey() const
{
	const char* sz = enums();
	if (sz == nullptr) return nullptr;
	if (sz[0] == '$') return nullptr;

	int n = value<int>();
	if (n < 0) return nullptr;
	for (int i = 0; i < n; ++i)
	{
		sz += strlen(sz) + 1;
		if ((sz == nullptr) || (*sz == 0)) return nullptr;
	}
	return sz;
}

//-----------------------------------------------------------------------------
const char* FEParam::units() const
{
	return m_szunit;
}

//-----------------------------------------------------------------------------
FEParam* FEParam::setUnits(const char* szunit) { m_szunit = szunit; return this; }

//-----------------------------------------------------------------------------
// set the enum values (\0 separated. Make sure the end of the string has two \0's)
FEParam* FEParam::setEnums(const char* sz) { m_szenum = sz; return this; }

//-----------------------------------------------------------------------------
FEParam* FEParam::setLongName(const char* sz)
{
	m_szlongname = sz; 
	return this;
}

//-----------------------------------------------------------------------------
// parameter dimension
int FEParam::dim() const { return m_dim; }

//-----------------------------------------------------------------------------
// parameter type
FEParamType FEParam::type() const { return m_type; }

//-----------------------------------------------------------------------------
// data pointer
void* FEParam::data_ptr() const { return m_pv; }

//-----------------------------------------------------------------------------
//! override the template for char pointers
char* FEParam::cvalue() { return (char*)data_ptr(); }

//-----------------------------------------------------------------------------
FEParamValue FEParam::paramValue(int i)
{
	switch (m_type)
	{
	case FE_PARAM_DOUBLE:
	{
		if (i == -1)
		{
			assert(m_dim == 1);
			return FEParamValue(this, &value<double>(), FE_PARAM_DOUBLE);
		}
		else return FEParamValue(this, &value<double>(i), FE_PARAM_DOUBLE);
	}
	break;
	case FE_PARAM_VEC3D:
	{
		if (i == -1)
		{
			assert(m_dim == 1);
			return FEParamValue(this, &value<vec3d>(), FE_PARAM_VEC3D);
		}
		else return FEParamValue(this, &value<vec3d>(i), FE_PARAM_VEC3D);
	}
	break;
	case FE_PARAM_DOUBLE_MAPPED:
	{
		if (i == -1)
		{
			assert(m_dim == 1);
			FEParamDouble& p = value<FEParamDouble>();
			if (p.isConst()) return FEParamValue(this, &p.constValue(), FE_PARAM_DOUBLE);
			else return FEParamValue(this, m_pv, m_type);
		}
		else
		{
			FEParamDouble& data = value<FEParamDouble>(i);
			return FEParamValue(this, &data, FE_PARAM_DOUBLE_MAPPED);
		}
	}
	break;
	case FE_PARAM_VEC3D_MAPPED:
	{
		if (i == -1)
		{
			assert(m_dim == 1);
			FEParamVec3& p = value<FEParamVec3>();
			if (p.isConst()) return FEParamValue(this, &p.constValue(), FE_PARAM_VEC3D);
			else return FEParamValue(this, m_pv, m_type);
		}
		else
		{
			FEParamVec3& data = value<FEParamVec3>(i);
			return FEParamValue(this, &data, FE_PARAM_VEC3D_MAPPED);
		}
	}
	break;
	case FE_PARAM_MAT3D_MAPPED:
	{
		if (i == -1)
		{
			assert(m_dim == 1);
			FEParamMat3d& p = value<FEParamMat3d>();
			if (p.isConst()) return FEParamValue(this, &p.constValue(), FE_PARAM_MAT3D);
			else return FEParamValue(this, m_pv, m_type);
		}
		else
		{
			FEParamMat3d& data = value<FEParamMat3d>(i);
			return FEParamValue(this, &data, FE_PARAM_MAT3D_MAPPED);
		}
	}
	break;
	case FE_PARAM_STD_VECTOR_DOUBLE:
	{
		vector<double>& data = value< vector<double> >();
		if ((i >= 0) && (i < (int)data.size()))
			return FEParamValue(this, &data[i], FE_PARAM_DOUBLE);
		else assert(false);
	}
	break;
	case FE_PARAM_STD_VECTOR_VEC2D:
	{
		vector<vec2d>& data = value< vector<vec2d> >();
		if ((i >= 0) && (i < (int)data.size()))
			return FEParamValue(this, &data[i], FE_PARAM_VEC2D);
		else assert(false);
	}
	break;
	case FE_PARAM_STD_VECTOR_STRING:
	{
		vector<string>& data = value< vector<string> >();
		if ((i >= 0) && (i < (int)data.size()))
			return FEParamValue(this, &data[i], FE_PARAM_STD_STRING);
	}
	break;
	default:
	{
		if (i == -1)
		{
			assert(m_dim == 1);
			return FEParamValue(this, m_pv, m_type);
		}
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
void FEParam::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << (int) m_type << m_flag;

		bool b = (m_watch ? *m_watch : false);
		ar << (b ? 1 : 0);

		if ((ar.IsShallow() == false) && (m_flag & FEParamFlag::FE_PARAM_USER))
		{
			assert(m_type == FE_PARAM_DOUBLE);
			ar << m_szname;
		}

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
			case FE_PARAM_STD_VECTOR_VEC2D: ar << value< std::vector<vec2d> >(); break;
			case FE_PARAM_DOUBLE_MAPPED:
			{
				FEParamDouble& p = value<FEParamDouble>();
				p.Serialize(ar);
			}
			break;
			case FE_PARAM_VEC3D_MAPPED:
			{
				FEParamVec3& p = value<FEParamVec3>();
				p.Serialize(ar);
			}
			break;
			case FE_PARAM_MAT3D_MAPPED:
			{
				FEParamMat3d& p = value<FEParamMat3d>();
				p.Serialize(ar);
			}
			break;
			case FE_PARAM_STD_VECTOR_INT:
			{
				vector<int>& p = value< vector<int> >();
				ar & p;
			}
			break;
			case FE_PARAM_STD_VECTOR_DOUBLE:
			{
				vector<double>& p = value< vector<double> >();
				ar & p;
			}
			break;
			default:
				assert(false);
			}
		}
		else
		{
			ar << m_dim;
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
			case FE_PARAM_DOUBLE_MAPPED:
			{
				FEParamDouble* p = (FEParamDouble*)(m_pv);
				for (int i = 0; i < m_dim; ++i)
				{
					p[i].Serialize(ar);
				}
			}
			break;
			default:
				assert(false);
			}
		}
	}
	else
	{
		int ntype;
		ar >> ntype >> m_flag;
		if (ntype != (int) m_type) throw DumpStream::ReadError();

		int watch = 0;
		ar >> watch;
		if (m_watch) *m_watch = (watch == 1);

		if ((ar.IsShallow() == false) && (m_flag & FEParamFlag::FE_PARAM_USER))
		{
			assert(m_type == FE_PARAM_DOUBLE);
			char sz[512] = { 0 };
			ar >> sz;
			m_szname = strdup(sz);
			m_pv = new double(0);
		}

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
			case FE_PARAM_STD_VECTOR_VEC2D: ar >> value< std::vector<vec2d> >(); break;
			case FE_PARAM_DOUBLE_MAPPED:
			{
				FEParamDouble& p = value<FEParamDouble>();
				p.Serialize(ar);
			}
			break;
			case FE_PARAM_VEC3D_MAPPED:
			{
				FEParamVec3& p = value<FEParamVec3>();
				p.Serialize(ar);
			}
			break;
			case FE_PARAM_MAT3D_MAPPED:
			{
				FEParamMat3d& p = value<FEParamMat3d>();
				p.Serialize(ar);
			}
			break;
			case FE_PARAM_STD_VECTOR_INT:
			{
				vector<int>& p = value< vector<int> >();
				ar & p;
			}
			break;
			case FE_PARAM_STD_VECTOR_DOUBLE:
			{
				vector<double>& p = value< vector<double> >();
				ar & p;
			}
			break;			default:
				assert(false);
			}
		}
		else
		{
			int ndim = 0;
			ar >> ndim;
			if (ndim != (int) m_dim) throw DumpStream::ReadError();
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
			case FE_PARAM_DOUBLE_MAPPED:
			{
				FEParamDouble* p = (FEParamDouble*)(m_pv);
				for (int i = 0; i < m_dim; ++i)
				{
					p[i].Serialize(ar);
				}
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

void FEParam::SaveClass(DumpStream& ar, FEParam* p)
{

}

FEParam* FEParam::LoadClass(DumpStream& ar, FEParam* p)
{
	return p;
}

//-----------------------------------------------------------------------------
//! This function copies the state of a parameter to this parameter.
//! This assumes that the parameters are compatible (i.e. have the same type)
//! This is used in FEParamContainer::CopyParameterListState()
bool FEParam::CopyState(const FEParam& p)
{
	if (p.type() != type()) return false;
	return true;
}

//-----------------------------------------------------------------------------
void FEParam::setParent(FEParamContainer* pc) { m_parent = pc; }
FEParamContainer* FEParam::parent() { return m_parent; }

//-----------------------------------------------------------------------------
void FEParam::SetWatchVariable(bool* watchvar)
{
	m_watch = watchvar;
}

//-----------------------------------------------------------------------------
bool* FEParam::GetWatchVariable()
{
	return m_watch;
}

//-----------------------------------------------------------------------------
void FEParam::SetWatchFlag(bool b)
{
	if (m_watch) *m_watch = b;
}

//-----------------------------------------------------------------------------
bool FEParam::IsHidden() const
{
	return (m_flag & FEParamFlag::FE_PARAM_HIDDEN);
}

//-----------------------------------------------------------------------------
bool FEParam::IsVolatile() const
{
	return (m_flag & FEParamFlag::FE_PARAM_VOLATILE);
}

//-----------------------------------------------------------------------------
FEParam* FEParam::MakeVolatile(bool b)
{
	if (b) m_flag = (m_flag | FEParamFlag::FE_PARAM_VOLATILE);
	else m_flag = (m_flag & ~FEParamFlag::FE_PARAM_VOLATILE);
	return this;
}

//-----------------------------------------------------------------------------
bool FEParam::IsTopLevel() const
{
	return (m_flag & FEParamFlag::FE_PARAM_TOPLEVEL);
}

//-----------------------------------------------------------------------------
FEParam* FEParam::MakeTopLevel(bool b)
{
	if (b) m_flag = (m_flag | FEParamFlag::FE_PARAM_TOPLEVEL);
	else m_flag = (m_flag & ~FEParamFlag::FE_PARAM_TOPLEVEL);
	return this;
}

//-----------------------------------------------------------------------------
FEParam* FEParam::SetFlags(unsigned int flags) { m_flag = flags; return this; }
unsigned int FEParam::GetFlags() const { return m_flag; }

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
			if      (paramName == "x") return FEParamValue(param, &v->x, FE_PARAM_DOUBLE);
			else if (paramName == "y") return FEParamValue(param, &v->y, FE_PARAM_DOUBLE);
			else if (paramName == "z") return FEParamValue(param, &v->z, FE_PARAM_DOUBLE);
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
