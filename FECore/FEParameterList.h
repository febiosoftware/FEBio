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



#pragma once
#include "vec2d.h"
#include "vec3d.h"
#include "mat3d.h"
#include <assert.h>
#include <list>
#include <memory>
#include "FEParam.h"
#include "FEParamValidator.h"
#include "ParamString.h"
#include "MathObject.h"
#include "fecore_api.h"
#include <stdio.h>
using namespace std;

//-----------------------------------------------------------------------------
class DumpStream;
class FEParamContainer;
class FEParamDouble;
class FEParamVec3;
class FEParamMat3d;
class FEParamMat3ds;
class FEDataArray;
class tens3drs;
class FEMaterialPointProperty;
class Image;

//-----------------------------------------------------------------------------
typedef std::list<FEParam>::iterator FEParamIterator;
typedef std::list<FEParam>::const_iterator FEParamIteratorConst;

//-----------------------------------------------------------------------------
//! A list of material parameters
class FECORE_API FEParameterList
{
public:
	FEParameterList(FEParamContainer* pc);
	virtual ~FEParameterList();

	//! assignment operator
	void operator = (FEParameterList& l);

	//! Add a parameter to the list
	FEParam* AddParameter(void* pv, FEParamType itype, int ndim, const char* sz, bool* watch = nullptr);

	//! Add a parameter to the list
	FEParam* AddParameter(void* pv, FEParamType type, int ndim, FEParamRange rng, double fmin, double fmax, const char* sz);

	//! find a parameter using the data pointer
	FEParam* FindFromData(void* pv);

	//! find a parameter using its name (the safe way)
	FEParam* FindFromName(const char* sz);

	//! get a parameter (the dangerous way)
	FEParam& operator [] (const char* sz) { return *FindFromName(sz); }

	//! returs the first parameter
	FEParamIterator first() { return m_pl.begin(); }

	//! returs the first parameter
	FEParamIteratorConst first() const { return m_pl.begin(); }

	//! number of parameters
	int Parameters() const { return (int) m_pl.size(); }

	//! return the parameter container
	FEParamContainer* GetContainer() { return m_pc; }

public:
	int SetActiveGroup(const char* szgroup);
	int GetActiveGroup();
	int ParameterGroups() const;
	const char* GetParameterGroupName(int i);

protected:
	FEParamContainer*	m_pc;	//!< parent container
	std::list<FEParam>		m_pl;	//!< the actual parameter list
	std::vector<const char*>	m_pg;	//!< parameter groups
	int	m_currentGroup;	//!< active parameter group (new parameters are assigned to the current group; can be -1)
};

//-----------------------------------------------------------------------------
//! Base class for classes that wish to support parameter lists
class FECORE_API FEParamContainer
{
public:
	//! constructor
	FEParamContainer();

	//! destructor
	virtual ~FEParamContainer();

	//! return the material's parameter list
	FEParameterList& GetParameterList();
	const FEParameterList& GetParameterList() const;

	//! find a parameter using it's name
	virtual FEParam* GetParameter(const char* szname);
	virtual FEParam* FindParameter(const ParamString& s);

	//! find a parameter using a pointer to the variable
	virtual FEParam* FindParameterFromData(void* pv);

	//! serialize parameter data
	virtual void Serialize(DumpStream& ar);

	//! validate material parameters.
	//! This function returns false on the first parameter encountered
	//! that is not valid (i.e. is outside its defined range). 
	//! Overload this function to do additional validation. Make sure to always call the base class.
	//! Use fecore_get_error_string() to find out which parameter failed validation.
	virtual bool Validate();

public:
	//! This copies the state of a parameter list (i.e. assigned load curve IDs)
	//! This function assumes that there is a one-to-one correspondence between
	//! source and target parameter lists.
	void CopyParameterListState(const FEParameterList& pl);

public:
	void BeginParameterGroup(const char* szname);
	void EndParameterGroup();

public:
	//! This function will be overridden by each class that defines a parameter list
	virtual void BuildParamList();

	//! Add a parameter to the list
	FEParam* AddParameter(void* pv, FEParamType itype, int ndim, const char* sz, bool* watch = nullptr);

	//! Add a parameter to the list
	FEParam* AddParameter(void* pv, FEParamType type, int ndim, RANGE rng, const char* sz);

public:
	FEParam* AddParameter(int&                 v, const char* sz);
	FEParam* AddParameter(bool&                v, const char* sz);
	FEParam* AddParameter(double&              v, const char* sz);
	FEParam* AddParameter(vec2d&               v, const char* sz);
	FEParam* AddParameter(vec3d&               v, const char* sz);
	FEParam* AddParameter(mat3d&               v, const char* sz);
	FEParam* AddParameter(mat3ds&              v, const char* sz);
	FEParam* AddParameter(FEParamDouble&       v, const char* sz);
	FEParam* AddParameter(FEParamVec3&         v, const char* sz);
	FEParam* AddParameter(FEParamMat3d&        v, const char* sz);
	FEParam* AddParameter(FEParamMat3ds&       v, const char* sz);
	FEParam* AddParameter(FEDataArray&         v, const char* sz);
	FEParam* AddParameter(tens3drs& 		   v, const char* sz);
	FEParam* AddParameter(std::string&         v, const char* sz);
	FEParam* AddParameter(std::vector<int>& v   , const char* sz);
	FEParam* AddParameter(std::vector<double>& v, const char* sz);
	FEParam* AddParameter(std::vector<vec2d>&  v, const char* sz);
	FEParam* AddParameter(std::vector<std::string>& v, const char* sz);
	FEParam* AddParameter(FEMaterialPointProperty& v, const char* sz);
	FEParam* AddParameter(MSimpleExpression& m, const char* sz);
	FEParam* AddParameter(Image& im           , const char* sz);

	FEParam* AddParameter(int&           v, RANGE rng, const char* sz);
	FEParam* AddParameter(double&        v, RANGE rng, const char* sz);
	FEParam* AddParameter(FEParamDouble& v, RANGE rng, const char* sz);

	FEParam* AddParameter(double& v, const char* sz, bool& watch);

	FEParam* AddParameter(int*           v, int ndim, const char* sz);
	FEParam* AddParameter(double*        v, int ndim, const char* sz);
	FEParam* AddParameter(FEParamDouble* v, int ndim, const char* sz);

	FEParam* AddParameter(int*           v, int ndim, RANGE rng, const char* sz);
	FEParam* AddParameter(double*        v, int ndim, RANGE rng, const char* sz);
	FEParam* AddParameter(FEParamDouble* v, int ndim, RANGE rng, const char* sz);

	FEParam* AddParameter(int& v, const char* sz, unsigned int flags, const char* szenum);
	FEParam* AddParameter(std::vector<int>& v, const char* sz, unsigned int flags, const char* szenum);
	FEParam* AddParameter(std::string& s, const char* sz, unsigned int flags, const char* szenum = nullptr);

	template <typename T> void SetParameter(const char* sz, T v);

private:
	FEParameterList*	m_pParam;	//!< parameter list
};

//-----------------------------------------------------------------------------
template <typename T> void FEParamContainer::SetParameter(const char* sz, T v)
{
	FEParam* p = m_pParam->FindFromName(sz); p->value<T>() = v;
}

//-----------------------------------------------------------------------------
// To add parameter list to a class, simply do the following two steps
// 1) add the DECLARE_FECORE_CLASS macro in the material class declaration
// 2) use the BEGIN_FECORE_CLASS, ADD_PARAM and END_FECORE_CLASS to
//    define a parameter list

// the following macro declares the parameter list for a material
#define DECLARE_FECORE_CLASS() \
	public: void BuildParamList() override;

#define FECORE_BASE_CLASS(theClass) \
	public: static const char* BaseClassName() { return #theClass; } \

// the BEGIN_FECORE_CLASS defines the beginning of a parameter list
#define BEGIN_FECORE_CLASS(theClass, baseClass) \
	void theClass::BuildParamList() { \
			baseClass::BuildParamList(); \

// the ADD_PARAMETER macro adds a parameter to the parameter list
#define ADD_PARAMETER(...) \
	AddParameter(__VA_ARGS__)

// the END_FECORE_CLASS defines the end of a parameter list
#define END_FECORE_CLASS() \
	}

// macro for starting a parameter group
#define BEGIN_PARAM_GROUP(a) BeginParameterGroup(a)

// macro for ending a parameter group
#define END_PARAM_GROUP()    EndParameterGroup()
