#pragma once
#include "FEMaterial.h"
#include "MathObject.h"
#include "FESolidDomain.h"
#include "FEDataMap.h"
#include "FEDomainList.h"

//---------------------------------------------------------------------------------------
// Base class for evaluating model parameters
template <class T> 
class FEValuator
{
public:
	FEValuator() {}
	virtual ~FEValuator() {}

	virtual T eval(const FEMaterialPoint& pt) = 0;
};

//---------------------------------------------------------------------------------------
class FEConstValue : public FEValuator<double>
{
public:
	FEConstValue(double v = 0.0) : m_val(v) {};
	double eval(const FEMaterialPoint& pt) override { return m_val; }

private:
	double	m_val;
};

//---------------------------------------------------------------------------------------
class FEMathExpression : public FEValuator<double>
{
public:
	FEMathExpression(const std::string& s);
	~FEMathExpression();
	double eval(const FEMaterialPoint& pt) override;

private:
	std::string			m_expr;
	MSimpleExpression	m_math;
};

//---------------------------------------------------------------------------------------
class FEMappedValue : public FEValuator<double>
{
public:
	FEMappedValue(FEDataMap* val);

	double eval(const FEMaterialPoint& pt) override;

private:
	FEDataMap*		m_val;
};

//---------------------------------------------------------------------------------------
// Base for model parameters.
class FEModelParam
{
public:
	FEModelParam();

	// set the domain
	void addDomain(FEDomain* dom) { m_dom.AddDomain(dom); }

	// get the domain list
	FEDomainList& getDomainList() { return m_dom;  }

	// set the scale factor
	void setScaleFactor(double s) { m_scl = s; }

protected:
	double			m_scl;	//!< scale factor. Used to store load curve value
	FEDomainList	m_dom;
};

//---------------------------------------------------------------------------------------
class FEParamDouble : public FEModelParam
{
public:
	FEParamDouble();

	// set the value
	void operator = (double v);

	// set the valuator
	void setValuator(FEValuator<double>* val);

	// evaluate the parameter at a material point
	double operator () (const FEMaterialPoint& pt) { return m_scl*m_val->eval(pt); }

private:
	FEParamDouble(const FEParamDouble& p) {}

private:
	FEValuator<double>*	m_val;
};

//=======================================================================================

//---------------------------------------------------------------------------------------
class FEConstValueVec3 : public FEValuator<vec3d>
{
public:
	FEConstValueVec3(const vec3d& r) : m_val(r) {};
	vec3d eval(const FEMaterialPoint& pt) override { return m_val; }

private:
	vec3d	m_val;
};

//---------------------------------------------------------------------------------------
class FEMathExpressionVec3 : public FEValuator<vec3d>
{
public:
	FEMathExpressionVec3(const std::string& sx, const std::string& sy, const std::string& sz);
	vec3d eval(const FEMaterialPoint& pt) override;

private:
	MSimpleExpression	m_math[3];
};

//---------------------------------------------------------------------------------------
class FEMappedValueVec3 : public FEValuator<vec3d>
{
public:
	FEMappedValueVec3(FEDataMap* val);

	vec3d eval(const FEMaterialPoint& pt) override;

private:
	FEDataMap*		m_val;
};

//---------------------------------------------------------------------------------------
class FEParamVec3 : public FEModelParam
{
public:
	FEParamVec3();

	// set the value
	void operator = (const vec3d& v);

	// set the valuator
	void setValuator(FEValuator<vec3d>* val);

	// evaluate the parameter at a material point
	vec3d operator () (const FEMaterialPoint& pt) { return m_val->eval(pt)*m_scl; }

private:
	FEParamVec3(const FEParamVec3& p) {}

private:
	FEValuator<vec3d>*	m_val;
};
