#pragma once
#include "FEMaterial.h"
#include "MathObject.h"
#include "FESolidDomain.h"
#include "FEDataMap.h"
#include "FEDomainList.h"
#include "FENodeDataMap.h"

//---------------------------------------------------------------------------------------
// Base class for evaluating model parameters
template <class T> class FEValuator
{
public:
	FEValuator() {}
	virtual ~FEValuator() {}

	virtual T operator()(const FEMaterialPoint& pt) = 0;

	virtual FEValuator<T>* copy() { return 0; }
};

//---------------------------------------------------------------------------------------
class FEConstValue : public FEValuator<double>
{
public:
	FEConstValue(double v = 0.0) : m_val(v) {};
	double operator()(const FEMaterialPoint& pt) override { return m_val; }

	double& value() { return m_val; }
	double value() const { return m_val;  }

	FEValuator<double>* copy() override { return new FEConstValue(m_val); }

private:
	double	m_val;
};

//---------------------------------------------------------------------------------------
class FEMathExpression : public FEValuator<double>
{
public:
	FEMathExpression() {}
	FEMathExpression(const std::string& s, FECoreBase* pc = 0);
	~FEMathExpression();
	double operator()(const FEMaterialPoint& pt) override;

	FEValuator<double>* copy() override;

private:
	std::string			m_expr;
	MSimpleExpression	m_math;
	std::vector<FEParam*>	m_vars;
};

//---------------------------------------------------------------------------------------
class FEMappedValue : public FEValuator<double>
{
public:
	FEMappedValue(FEDataMap* val, double scl = 1.0);

	double operator()(const FEMaterialPoint& pt) override;

	FEValuator<double>* copy() override;

	double& GetScale() { return m_scale; }

private:
	double		m_scale;
	FEDataMap*	m_val;
};

//---------------------------------------------------------------------------------------
class FENodeMappedValue : public FEValuator<double>
{
public:
	FENodeMappedValue(FENodeDataMap* val);

	double operator()(const FEMaterialPoint& pt) override;

	FEValuator<double>* copy() override {
		return new FENodeMappedValue(m_val);
	}

private:
	FENodeDataMap*		m_val;
};

//---------------------------------------------------------------------------------------
// Base for model parameters.
class FEModelParam
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
class FEParamDouble : public FEModelParam
{
public:
	FEParamDouble();

	FEParamDouble(const FEParamDouble& p);

	// set the value
	void operator = (double v);

	// set the valuator
	void setValuator(FEValuator<double>* val);

	// evaluate the parameter at a material point
	double operator () (const FEMaterialPoint& pt) { return m_scl*(*m_val)(pt); }

	// is this a const value
	bool isConst() const;

	// get the const value (returns 0 if param is not const)
	double& constValue();
	double constValue() const;

private:
	FEValuator<double>*	m_val;
};

//=======================================================================================

//---------------------------------------------------------------------------------------
class FEConstValueVec3 : public FEValuator<vec3d>
{
public:
	FEConstValueVec3(const vec3d& r) : m_val(r) {};
	vec3d operator()(const FEMaterialPoint& pt) override { return m_val; }

	FEValuator<vec3d>* copy() override { return new FEConstValueVec3(m_val); }

	vec3d& value() { return m_val; }

private:
	vec3d	m_val;
};

//---------------------------------------------------------------------------------------
class FEMathExpressionVec3 : public FEValuator<vec3d>
{
public:
	FEMathExpressionVec3() {}
	FEMathExpressionVec3(const std::string& sx, const std::string& sy, const std::string& sz);
	vec3d operator()(const FEMaterialPoint& pt) override;

	FEValuator<vec3d>* copy() override;

private:
	MSimpleExpression	m_math[3];
};

//---------------------------------------------------------------------------------------
class FEMappedValueVec3 : public FEValuator<vec3d>
{
public:
	FEMappedValueVec3(FEDataMap* val, vec3d scl = vec3d(1,1,1));

	vec3d operator()(const FEMaterialPoint& pt) override;

	FEValuator<vec3d>* copy() override;

	vec3d& GetScale() { return m_scale; }

private:
	vec3d			m_scale;
	FEDataMap*		m_val;
};

//---------------------------------------------------------------------------------------
class FEParamVec3 : public FEModelParam
{
public:
	FEParamVec3();

	FEParamVec3(const FEParamVec3& p);

	// set the value
	void operator = (const vec3d& v);

	// set the valuator
	void setValuator(FEValuator<vec3d>* val);

	// evaluate the parameter at a material point
	vec3d operator () (const FEMaterialPoint& pt) { return (*m_val)(pt)*m_scl; }

	// is this a const
	bool isConst() const;

	vec3d& constValue();

private:
	FEValuator<vec3d>*	m_val;
};
