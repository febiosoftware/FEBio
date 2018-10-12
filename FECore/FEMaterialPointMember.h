#pragma once
#include "FECoreBase.h"
#include "FEMaterialPoint.h"
#include "fecore_type.h"

class FEMaterialPointMember
{
public:
	FEMaterialPointMember(FEDataType dataType) : m_type(dataType) {}
	virtual ~FEMaterialPointMember() {}

	FEDataType dataType() const { return m_type; }

	virtual void set(FEMaterialPoint& mp, const double& v) {}
	virtual void set(FEMaterialPoint& mp, const vec3d&  v) {}
	virtual void set(FEMaterialPoint& mp, const mat3d&  v) {}

	virtual void get(FEMaterialPoint& mp, double& v) {}
	virtual void get(FEMaterialPoint& mp, vec3d&  v) {}
	virtual void get(FEMaterialPoint& mp, mat3d&  v) {}

private:
	FEDataType	m_type;
};


template <class T, class A> class FEMaterialPointMember_T : public FEMaterialPointMember
{
public:
	FEMaterialPointMember_T() : FEMaterialPointMember(fecoreType<A>::type()) {}
	void set(FEMaterialPoint& mp, const A& Q)
	{
		T& pt = *mp.ExtractData<T>();
		data(pt) = Q;
	}

	void get(FEMaterialPoint& mp, A& Q)
	{
		T& pt = *mp.ExtractData<T>();
		Q = data(pt);
	}

	virtual A& data(T& pt) = 0;
};
