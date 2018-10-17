#pragma once
#include "FECoreBase.h"
#include "FEMaterialPoint.h"
#include "fecore_type.h"

//--------------------------------------------------------------------
// This class defines a material point property. Material point properties
// are used to move data in and out of material points.
class FEMaterialPointProperty
{
public:
	FEMaterialPointProperty(FEDataType dataType) : m_type(dataType) {}
	virtual ~FEMaterialPointProperty() {}

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

//--------------------------------------------------------------------
// Template class for defining material point property classes. 
// Derived classes need to override the data member. 
template <class T, class A> class FEMaterialPointProperty_T : public FEMaterialPointProperty
{
public:
	FEMaterialPointProperty_T() : FEMaterialPointProperty(fecoreType<A>::type()) {}
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
