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
#include "FEParameterList.h"
#include "FECoreKernel.h"
#include "DumpStream.h"
#include "tens3d.h"
#include "FEModelParam.h"
#include <string>
#include <assert.h>
using namespace std;

//-----------------------------------------------------------------------------
FEParameterList::FEParameterList(FEParamContainer* pc) : m_pc(pc) 
{
	m_currentGroup = -1;
}

//-----------------------------------------------------------------------------
FEParameterList::~FEParameterList() 
{
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

		if (s.type() != d.type()) { assert(false); return; }
		if (s.dim() != d.dim()) { assert(false); return; }
		if (s.dim() == 1)
		{
			switch (s.type())
			{
			case FE_PARAM_INT   : d.value<int   >() = s.value<int   >(); break;
			case FE_PARAM_BOOL  : d.value<bool  >() = s.value<bool  >(); break;
			case FE_PARAM_DOUBLE: d.value<double>() = s.value<double>(); break;
			case FE_PARAM_VEC3D : d.value<vec3d >() = s.value<vec3d >(); break;
			case FE_PARAM_MAT3D : d.value<mat3d >() = s.value<mat3d >(); break;
			case FE_PARAM_MAT3DS: d.value<mat3ds>() = s.value<mat3ds>(); break;
			case FE_PARAM_TENS3DRS: d.value<tens3drs>() = s.value<tens3drs>(); break;
			case FE_PARAM_STD_STRING: d.value<std::string>() = s.value<std::string>(); break;
			case FE_PARAM_DOUBLE_MAPPED:
			{
				FEParamDouble& mat3d = d.value<FEParamDouble>();
				FEParamDouble& src = s.value<FEParamDouble>();
				mat3d.setValuator(src.valuator()->copy());
			}
			break;
			case FE_PARAM_MAT3D_MAPPED:
			{
				FEParamMat3d& mat3d = d.value<FEParamMat3d>();
				FEParamMat3d& src = s.value<FEParamMat3d>();
				mat3d.setValuator(src.valuator()->copy());
			}
			break;
			case FE_PARAM_STD_VECTOR_VEC2D:
			{
				d.value< std::vector<vec2d> >() = s.value< std::vector<vec2d> >();
			}
			break;
			default:
				assert(false);
			}
		}
		else
		{
			switch (s.type())
			{
			case FE_PARAM_INT:
				{
					for (int i=0; i<s.dim(); ++i) d.pvalue<int>()[i] = s.pvalue<int>()[i];
				}
				break;
			case FE_PARAM_DOUBLE:
				{
					for (int i=0; i<s.dim(); ++i) d.pvalue<double>()[i] = s.pvalue<double>()[i];
				}
				break;
			case FE_PARAM_DOUBLE_MAPPED:
				{
					for (int i=0; i<s.dim(); ++i)
					{
						FEParamDouble& pd = d.value<FEParamDouble>(i);
						FEParamDouble& ps = s.value<FEParamDouble>(i);
						assert(ps.isConst());
						pd = ps.constValue();
					}
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
FEParam* FEParameterList::AddParameter(void *pv, FEParamType itype, int ndim, const char *sz, bool* watch)
{
	// sanity checks
	assert(pv);
	assert(sz);

	// create a new parameter object
	FEParam p(pv, itype, ndim, sz, watch);
	p.SetParamGroup(m_currentGroup);

	// add the parameter to the list
	m_pl.push_back(p);

	return &(m_pl.back());
}
//-----------------------------------------------------------------------------
// This function adds a parameter to the parameter list
FEParam* FEParameterList::AddParameter(void *pv, FEParamType itype, int ndim, FEParamRange rng, double fmin, double fmax, const char *sz)
{
	assert(pv);
	assert(sz);

	// create a new parameter object
	FEParam p(pv, itype, ndim, sz);
	p.SetParamGroup(m_currentGroup);

	// set the range
	// (range checking is only supported for int and double params)
	if (rng != FE_DONT_CARE)
	{
		if      (itype == FE_PARAM_INT) p.SetValidator(new FEIntValidator(rng, (int) fmin, (int) fmax));
		else if (itype == FE_PARAM_DOUBLE) p.SetValidator(new FEDoubleValidator(rng, fmin, fmax));
		else if (itype == FE_PARAM_DOUBLE_MAPPED) p.SetValidator(new FEParamDoubleValidator(rng, fmin, fmax));
	}

	// add the parameter to the list
	m_pl.push_back(p);

	return &(m_pl.back());
}

//-----------------------------------------------------------------------------
// Find a parameter using its data pointer
FEParam* FEParameterList::FindFromData(void* pv)
{
	FEParam* pp = 0;
	if (m_pl.empty() == false)
	{
		list<FEParam>::iterator it;
		for (it = m_pl.begin(); it != m_pl.end(); ++it)
		{
			if (it->dim() <= 1)
			{
				if (it->data_ptr() == pv)
				{
					pp = &(*it);
					return pp;
				}
			}
			else
			{
				for (int i = 0; i < it->dim(); ++i)
				{
					void* pd = nullptr;
					switch (it->type())
					{
					case FE_PARAM_DOUBLE_MAPPED: pd = &(it->value<FEParamDouble>(i)); break;
					default:
						assert(false);
					}
					if (pv == pd)
					{
						pp = &(*it);
						return pp;
					}
				}
			}
		}
	}
	return pp;
}

//-----------------------------------------------------------------------------
// This function searches the parameters in the list for a parameter
// with the name given by the input argument
// \param sz name of parameter to find
FEParam* FEParameterList::FindFromName(const char* sz)
{
	if (sz == 0) return 0;

	FEParam* pp = 0;
	if (m_pl.size() > 0)
	{
		list<FEParam>::iterator it;
		for (it = m_pl.begin(); it != m_pl.end(); ++it)
		{
			if (strcmp(it->name(), sz) == 0)
			{
				pp = &(*it);
				break;
			}
		}
	}

	return pp;
}

//-----------------------------------------------------------------------------
int FEParameterList::SetActiveGroup(const char* szgroup)
{
	if (szgroup == nullptr) m_currentGroup = -1;
	else
	{
		m_currentGroup = -1;
		for (size_t i = 0; i < m_pg.size(); ++i)
		{
			if (strcmp(m_pg[i], szgroup) == 0)
			{
				m_currentGroup = i;
			}
		}
		if (m_currentGroup == -1)
		{
			m_currentGroup = (int)m_pg.size();
			m_pg.push_back(szgroup);
		}
	}
	return m_currentGroup;
}

//-----------------------------------------------------------------------------
int FEParameterList::GetActiveGroup()
{
	return m_currentGroup;
}

//-----------------------------------------------------------------------------
int FEParameterList::ParameterGroups() const
{
	return (int)m_pg.size();
}

//-----------------------------------------------------------------------------
const char* FEParameterList::GetParameterGroupName(int i)
{
	return m_pg[i];
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
const FEParameterList& FEParamContainer::GetParameterList() const
{
	assert(m_pParam);
	return *m_pParam;
}

//-----------------------------------------------------------------------------
// Find a parameter from its name
FEParam* FEParamContainer::GetParameter(const char* szname)
{
	FEParameterList& pl = GetParameterList();
	return pl.FindFromName(szname);
}

//-----------------------------------------------------------------------------
// Find a parameter from its name
FEParam* FEParamContainer::FindParameter(const ParamString& s)
{
	FEParameterList& pl = GetParameterList();
	return pl.FindFromName(s.c_str());
}

//-----------------------------------------------------------------------------
FEParam* FEParamContainer::FindParameterFromData(void* pv)
{
	FEParameterList& pl = GetParameterList();
	return pl.FindFromData(pv);
}

//-----------------------------------------------------------------------------
//! This function will be overridden by each class that defines a parameter list
void FEParamContainer::BuildParamList() 
{
}

//-----------------------------------------------------------------------------
// Add a parameter to the parameter list
FEParam* FEParamContainer::AddParameter(void* pv, FEParamType itype, int ndim, const char* sz, bool* watch)
{
	assert(m_pParam);
	FEParam* p = m_pParam->AddParameter(pv, itype, ndim, sz, watch);
	p->setParent(this);
	return p;
}

//-----------------------------------------------------------------------------
// Add a parameter to the parameter list
FEParam* FEParamContainer::AddParameter(void* pv, FEParamType itype, int ndim, RANGE rng, const char* sz)
{
	assert(m_pParam);
	FEParam* p = m_pParam->AddParameter(pv, itype, ndim, rng.m_rt, rng.m_fmin, rng.m_fmax, sz);
	p->setParent(this);
	return p;
}

//-----------------------------------------------------------------------------
FEParam* FEParamContainer::AddParameter(int&                      v, const char* sz) { return AddParameter(&v, FE_PARAM_INT, 1, sz); }
FEParam* FEParamContainer::AddParameter(bool&                     v, const char* sz) { return AddParameter(&v, FE_PARAM_BOOL, 1, sz); }
FEParam* FEParamContainer::AddParameter(double&                   v, const char* sz) { return AddParameter(&v, FE_PARAM_DOUBLE, 1, sz); }
FEParam* FEParamContainer::AddParameter(vec2d&                    v, const char* sz) { return AddParameter(&v, FE_PARAM_VEC2D, 1, sz); }
FEParam* FEParamContainer::AddParameter(vec3d&                    v, const char* sz) { return AddParameter(&v, FE_PARAM_VEC3D, 1, sz); }
FEParam* FEParamContainer::AddParameter(mat3d&                    v, const char* sz) { return AddParameter(&v, FE_PARAM_MAT3D, 1, sz); }
FEParam* FEParamContainer::AddParameter(mat3ds&                   v, const char* sz) { return AddParameter(&v, FE_PARAM_MAT3DS, 1, sz); }
FEParam* FEParamContainer::AddParameter(FEParamDouble&            v, const char* sz) { return AddParameter(&v, FE_PARAM_DOUBLE_MAPPED, 1, sz); }
FEParam* FEParamContainer::AddParameter(FEParamVec3&              v, const char* sz) { return AddParameter(&v, FE_PARAM_VEC3D_MAPPED, 1, sz); }
FEParam* FEParamContainer::AddParameter(FEParamMat3d&             v, const char* sz) { return AddParameter(&v, FE_PARAM_MAT3D_MAPPED, 1, sz); }
FEParam* FEParamContainer::AddParameter(FEParamMat3ds&            v, const char* sz) { return AddParameter(&v, FE_PARAM_MAT3DS_MAPPED, 1, sz); }
FEParam* FEParamContainer::AddParameter(FEDataArray&              v, const char* sz) { return AddParameter(&v, FE_PARAM_DATA_ARRAY, 1, sz); }
FEParam* FEParamContainer::AddParameter(tens3drs& 		          v, const char* sz) { return AddParameter(&v, FE_PARAM_TENS3DRS, 1, sz); }
FEParam* FEParamContainer::AddParameter(std::string&              v, const char* sz) { return AddParameter(&v, FE_PARAM_STD_STRING, 1, sz); }
FEParam* FEParamContainer::AddParameter(std::vector<int>&         v, const char* sz) { return AddParameter(&v, FE_PARAM_STD_VECTOR_INT, 1, sz); }
FEParam* FEParamContainer::AddParameter(std::vector<double>&      v, const char* sz) { return AddParameter(&v, FE_PARAM_STD_VECTOR_DOUBLE, 1, sz); }
FEParam* FEParamContainer::AddParameter(std::vector<vec2d>&       v, const char* sz) { return AddParameter(&v, FE_PARAM_STD_VECTOR_VEC2D, 1, sz); }
FEParam* FEParamContainer::AddParameter(std::vector<std::string>& v, const char* sz) { return AddParameter(&v, FE_PARAM_STD_VECTOR_STRING, 1, sz); }
FEParam* FEParamContainer::AddParameter(FEMaterialPointProperty&  v, const char* sz) { return AddParameter(&v, FE_PARAM_MATERIALPOINT, 1, sz); }
//FEParam* FEParamContainer::AddParameter(Image& v                   , const char* sz) { return AddParameter(&v, FE_PARAM_IMAGE_3D, 1, sz); }

FEParam* FEParamContainer::AddParameter(int&           v, RANGE rng, const char* sz) { return AddParameter(&v, FE_PARAM_INT, 1, rng, sz); }
FEParam* FEParamContainer::AddParameter(double&        v, RANGE rng, const char* sz) { return AddParameter(&v, FE_PARAM_DOUBLE, 1, rng, sz); }
FEParam* FEParamContainer::AddParameter(FEParamDouble& v, RANGE rng, const char* sz) { return AddParameter(&v, FE_PARAM_DOUBLE_MAPPED, 1, rng, sz); }

FEParam* FEParamContainer::AddParameter(double&        v, const char* sz, bool& watch) { return AddParameter(&v, FE_PARAM_DOUBLE, 1, sz, &watch); }

FEParam* FEParamContainer::AddParameter(int*           v, int ndim, const char* sz) { return AddParameter(v, FE_PARAM_INT, ndim, sz); }
FEParam* FEParamContainer::AddParameter(double*        v, int ndim, const char* sz) { return AddParameter(v, FE_PARAM_DOUBLE, ndim, sz); }
FEParam* FEParamContainer::AddParameter(FEParamDouble* v, int ndim, const char* sz) { return AddParameter(v, FE_PARAM_DOUBLE_MAPPED, ndim, sz); }

FEParam* FEParamContainer::AddParameter(int*           v, int ndim, RANGE rng, const char* sz) { return AddParameter(v, FE_PARAM_INT, ndim, rng, sz); }
FEParam* FEParamContainer::AddParameter(double*        v, int ndim, RANGE rng, const char* sz) { return AddParameter(v, FE_PARAM_DOUBLE, ndim, rng, sz); }
FEParam* FEParamContainer::AddParameter(FEParamDouble* v, int ndim, RANGE rng, const char* sz) { return AddParameter(v, FE_PARAM_DOUBLE_MAPPED, ndim, rng, sz); }

//-----------------------------------------------------------------------------
FEParam* FEParamContainer::AddParameter(int& v, const char* sz, unsigned int flags, const char* szenum)
{
	FEParam* p = AddParameter(&v, FE_PARAM_INT, 1, sz);
	p->setParent(this);
	p->SetFlags(flags);
	p->setEnums(szenum);
	return p;
}

//-----------------------------------------------------------------------------
FEParam* FEParamContainer::AddParameter(std::vector<int>& v, const char* sz, unsigned int flags, const char* szenum)
{
	FEParam* p = AddParameter(&v, FE_PARAM_STD_VECTOR_INT, 1, sz);
	p->setParent(this);
	p->SetFlags(flags);
	p->setEnums(szenum);
	return p;
}

//-----------------------------------------------------------------------------
FEParam* FEParamContainer::AddParameter(std::string& s, const char* sz, unsigned int flags, const char* szenum)
{
	FEParam* p = AddParameter(&s, FE_PARAM_STD_STRING, 1, sz);
	p->setParent(this);
	p->SetFlags(flags);
	p->setEnums(szenum);
	return p;
}

//-----------------------------------------------------------------------------
// Serialize parameters to archive
void FEParamContainer::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

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
			ar << p;
		}
	}
	else
	{
		int NP = 0;
		ar >> NP;
		if (NP)
		{
			FEParameterList& pl = GetParameterList();
			if (NP != pl.Parameters()) throw DumpStream::ReadError();
			list<FEParam>::iterator it = pl.first();
			for (int i=0; i<NP; ++i)
			{
				FEParam& p = *it++;
				ar >> p;
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function validates all parameters. 
//! It fails on the first parameter that is outside its allowed range.
//! Use fecore_get_error_string() to find out which parameter failed validation.
bool FEParamContainer::Validate()
{
	FEParameterList& pl = GetParameterList();
	int N = pl.Parameters();
	list<FEParam>::iterator pi = pl.first();
	for (int i=0; i<N; ++i, pi++)
	{
		FEParam& p = *pi;
		if (p.is_valid() == false)
		{
            return false;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEParamContainer::CopyParameterListState(const FEParameterList& pl)
{
	FEParameterList& pl_this = GetParameterList();
	assert(pl_this.Parameters() == pl.Parameters());
	int NP = pl.Parameters();
	FEParamIteratorConst it_s = pl.first();
	FEParamIterator it_d = pl_this.first();
	for (int i=0; i<NP; ++i, ++it_s, ++it_d)
	{
		const FEParam& ps = *it_s;
		FEParam& pd = *it_d;
		if (pd.CopyState(ps) == false) { assert(false); }
	}
}

//-----------------------------------------------------------------------------
void FEParamContainer::BeginParameterGroup(const char* szname)
{
	FEParameterList& pl_this = GetParameterList();
	pl_this.SetActiveGroup(szname);
}

//-----------------------------------------------------------------------------
void FEParamContainer::EndParameterGroup()
{
	FEParameterList& pl_this = GetParameterList();
	pl_this.SetActiveGroup(nullptr);
}
