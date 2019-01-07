#pragma once
#include "FEScalarValuator.h"
#include "FEVec3dValuator.h"
#include "FEMat3dValuator.h"
#include "FEItemList.h"

//---------------------------------------------------------------------------------------
// Base for model parameters.
class FECORE_API FEModelParam
{
public:
	FEModelParam();
	virtual ~FEModelParam() {}

	// set the domain
	void SetItemList(FEItemList* itemList) { m_dom = itemList; }

	// get the domain list
	FEItemList* GetItemList() { return m_dom;  }

	// set the scale factor
	void SetScaleFactor(double s) { m_scl = s; }

protected:
	double			m_scl;	//!< scale factor. Used to store load curve value
	FEItemList*		m_dom;
};

//---------------------------------------------------------------------------------------
class FECORE_API FEParamDouble : public FEModelParam
{
public:
	FEParamDouble();

	FEParamDouble(const FEParamDouble& p);

	// set the value
	void operator = (double v);

	// set the valuator
	void setValuator(FEScalarValuator* val);

	// evaluate the parameter at a material point
	double operator () (const FEMaterialPoint& pt) { return m_scl*(*m_val)(pt); }

	// is this a const value
	bool isConst() const;

	// get the const value (returns 0 if param is not const)
	double& constValue();
	double constValue() const;

private:
	FEScalarValuator*	m_val;
};

//=======================================================================================

//---------------------------------------------------------------------------------------
class FECORE_API FEParamVec3 : public FEModelParam
{
public:
	FEParamVec3();

	FEParamVec3(const FEParamVec3& p);

	// set the value
	void operator = (const vec3d& v);

	// set the valuator
	void setValuator(FEVec3dValuator* val);

	// evaluate the parameter at a material point
	vec3d operator () (const FEMaterialPoint& pt) { return (*m_val)(pt)*m_scl; }

	// return a unit vector
	vec3d unitVector(const FEMaterialPoint& pt) { return (*this)(pt).normalized(); }

	// is this a const
	bool isConst() const { return m_val->isConst(); }

	vec3d& constValue() { return *m_val->constValue(); };

private:
	FEVec3dValuator*	m_val;
};

//=======================================================================================

//---------------------------------------------------------------------------------------
class FECORE_API FEParamMat3d : public FEModelParam
{
public:
	FEParamMat3d();

	FEParamMat3d(const FEParamMat3d& p);

	// set the value
	void operator = (const mat3d& v);

	// set the valuator
	void setValuator(FEMat3dValuator* val);

	// evaluate the parameter at a material point
	mat3d operator () (const FEMaterialPoint& pt) { return (*m_val)(pt)*m_scl; }

	// is this a const
	bool isConst() const { return m_val->isConst(); }

	mat3d& constValue() { return *m_val->constValue(); };

private:
	FEMat3dValuator*	m_val;
};
