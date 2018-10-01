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
#include "fecore_api.h"
#include <stdio.h>
using namespace std;

//-----------------------------------------------------------------------------
class DumpStream;
class FEParamContainer;
class FEParamDouble;
class FEParamVec3;
class FEDataArray;
class FEFunction1D;
class tens3drs;

//-----------------------------------------------------------------------------
typedef list<FEParam>::iterator FEParamIterator;
typedef list<FEParam>::const_iterator FEParamIteratorConst;

//-----------------------------------------------------------------------------
//! A list of material parameters
class FECORE_API FEParameterList
{
public:
	FEParameterList(FEParamContainer* pc) : m_pc(pc) {}
	virtual ~FEParameterList(){}

	//! assignment operator
	void operator = (FEParameterList& l);

	//! Add a parameter to the list
	FEParam* AddParameter(void* pv, FEParamType itype, int ndim, const char* sz);

	//! Add a parameter to the list
	FEParam* AddParameter(void* pv, FEParamType type, int ndim, FEParamRange rng, double fmin, double fmax, const char* sz);

	//! find a parameter using the data pointer
	FEParam* FindFromData(void* pv);

	//! find a parameter using it's name (the safe way)
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

protected:
	FEParamContainer*	m_pc;	//!< parent container
	list<FEParam>		m_pl;	//!< the actual parameter list
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
	//! This function is called after the parameter was read in from the input file.
	//! It can be used to do additional processing when a parameter is read in.
	virtual void SetParameter(FEParam& p) {}

	//! If a parameter has attributes, this function will be called
	virtual bool SetParameterAttribute(FEParam& p, const char* szatt, const char* szval) { return false; }

	//! This copies the state of a parameter list (i.e. assigned load curve IDs)
	//! This function assumes that there is a one-to-one correspondence between
	//! source and target parameter lists.
	void CopyParameterListState(const FEParameterList& pl);

public:
	//! This function will be overridden by each class that defines a parameter list
	virtual void BuildParamList() {}

	//! Add a parameter to the list
	void AddParameter(void* pv, FEParamType itype, int ndim, const char* sz);

	//! Add a parameter to the list
	void AddParameter(void* pv, FEParamType type, int ndim, RANGE rng, const char* sz);

public:
	void AddParameter(int&                 v, const char* sz) { AddParameter(&v, FE_PARAM_INT              , 1, sz); }
	void AddParameter(bool&                v, const char* sz) { AddParameter(&v, FE_PARAM_BOOL             , 1, sz); }
	void AddParameter(double&              v, const char* sz) { AddParameter(&v, FE_PARAM_DOUBLE           , 1, sz); }
	void AddParameter(vec2d&               v, const char* sz) { AddParameter(&v, FE_PARAM_VEC2D            , 1, sz); }
	void AddParameter(vec3d&               v, const char* sz) { AddParameter(&v, FE_PARAM_VEC3D            , 1, sz); }
	void AddParameter(mat3d&               v, const char* sz) { AddParameter(&v, FE_PARAM_MAT3D            , 1, sz); }
	void AddParameter(mat3ds&              v, const char* sz) { AddParameter(&v, FE_PARAM_MAT3DS           , 1, sz); }
	void AddParameter(FEParamDouble&       v, const char* sz) { AddParameter(&v, FE_PARAM_DOUBLE_MAPPED    , 1, sz); }
	void AddParameter(FEParamVec3&         v, const char* sz) { AddParameter(&v, FE_PARAM_VEC3D_MAPPED     , 1, sz); }
	void AddParameter(FEDataArray&         v, const char* sz) { AddParameter(&v, FE_PARAM_DATA_ARRAY       , 1, sz); }
	void AddParameter(FEFunction1D&        v, const char* sz) { AddParameter(&v, FE_PARAM_FUNC1D           , 1, sz); }
	void AddParameter(tens3drs& 		   v, const char* sz) { AddParameter(&v, FE_PARAM_TENS3DRS         , 1, sz); }
	void AddParameter(std::string&         v, const char* sz) { AddParameter(&v, FE_PARAM_STRING           , 1, sz); }
	void AddParameter(std::vector<double>& v, const char* sz) { AddParameter(&v, FE_PARAM_STD_VECTOR_DOUBLE, 1, sz);  }

	void AddParameter(int&           v, RANGE rng, const char* sz) { AddParameter(&v, FE_PARAM_INT          , 1, rng, sz); }
	void AddParameter(double&        v, RANGE rng, const char* sz) { AddParameter(&v, FE_PARAM_DOUBLE       , 1, rng, sz); }
	void AddParameter(FEParamDouble& v, RANGE rng, const char* sz) { AddParameter(&v, FE_PARAM_DOUBLE_MAPPED, 1, rng, sz); }

	void AddParameter(int*           v, int ndim, const char* sz) { AddParameter(v, FE_PARAM_INT          , ndim, sz); }
	void AddParameter(bool*          v, int ndim, const char* sz) { AddParameter(v, FE_PARAM_BOOL         , ndim, sz); }
	void AddParameter(double*        v, int ndim, const char* sz) { AddParameter(v, FE_PARAM_DOUBLE       , ndim, sz); }
	void AddParameter(vec2d*         v, int ndim, const char* sz) { AddParameter(v, FE_PARAM_VEC2D        , ndim, sz); }
	void AddParameter(vec3d*         v, int ndim, const char* sz) { AddParameter(v, FE_PARAM_VEC3D        , ndim, sz); }
	void AddParameter(mat3d*         v, int ndim, const char* sz) { AddParameter(v, FE_PARAM_MAT3D        , ndim, sz); }
	void AddParameter(mat3ds*        v, int ndim, const char* sz) { AddParameter(v, FE_PARAM_MAT3DS       , ndim, sz); }
	void AddParameter(std::string*   v, int ndim, const char* sz) { AddParameter(v, FE_PARAM_STRING       , ndim, sz); }

	void AddParameter(int*           v, int ndim, RANGE rng, const char* sz) { AddParameter(v, FE_PARAM_INT          , ndim, rng, sz); }
	void AddParameter(double*        v, int ndim, RANGE rng, const char* sz) { AddParameter(v, FE_PARAM_DOUBLE       , ndim, rng, sz); }

private:
	FEParameterList*	m_pParam;	//!< parameter list
};

//-----------------------------------------------------------------------------
// To add parameter list to a material, simply do the following two steps
// 1) add the DECLARE_PARAMETER_LIST macro in the material class declaration
// 2) use the BEGIN_PARAMETER_LIST, ADD_PARAM and END_PARAMETER_LIST to
//    define a parameter list

// the following macro declares the parameter list for a material
#define DECLARE_PARAMETER_LIST() \
protected: \
	void BuildParamList() override;

// the BEGIN_PARAMETER_LIST defines the beginning of a parameter list
#define BEGIN_PARAMETER_LIST(theClass, baseClass) \
	void theClass::BuildParamList() { \
			baseClass::BuildParamList(); \

// the ADD_PARAMETER macro adds a parameter to the parameter list
#define ADD_PARAMETER(...) \
	AddParameter(__VA_ARGS__);

// the END_PARAMETER_LIST defines the end of a parameter list
#define END_PARAMETER_LIST() \
	}
