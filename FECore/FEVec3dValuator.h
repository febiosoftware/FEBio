#pragma once
#include "FEValuator.h"
#include "FEDataMap.h"
#include "MathObject.h"

//---------------------------------------------------------------------------------------
// Base class for evaluating vec3d parameters
class FECORE_API FEVec3dValuator : public FEValuator
{
	FECORE_SUPER_CLASS

public:
	// evaluate value at a material point
	virtual vec3d operator()(const FEMaterialPoint& pt) = 0;

	// create a copy of the valuator
	virtual FEVec3dValuator* copy() = 0;

	// is the valuator constant
	virtual bool isConst() { return false; }

	// return the const value
	virtual vec3d* constValue() { return nullptr; }

public:
	FEVec3dValuator(FEModel* fem) : FEValuator(fem) {};
};

//---------------------------------------------------------------------------------------
// A constant valuator
class FECORE_API FEConstValueVec3 : public FEVec3dValuator
{
public:
	FEConstValueVec3(FEModel* fem);

	FEVec3dValuator* copy() override;

	vec3d operator()(const FEMaterialPoint& pt) override { return m_val; }

	// is this a const value
	bool isConst() override { return true; }

	// get the const value (returns 0 if param is not const)
	vec3d* constValue() override { return &m_val; }

	vec3d& value() { return m_val; }

private:
	vec3d	m_val;

	DECLARE_FECORE_CLASS();
};

//---------------------------------------------------------------------------------------
// The value is calculated using a mathematical expression
class FECORE_API FEMathValueVec3 : public FEVec3dValuator
{
public:
	FEMathValueVec3(FEModel* fem);
	vec3d operator()(const FEMaterialPoint& pt) override;

	bool Init() override;

	bool create(const std::string& sx, const std::string& sy, const std::string& sz);

	FEVec3dValuator* copy() override;

private:
	std::string			m_expr;
	MSimpleExpression	m_math[3];

	DECLARE_FECORE_CLASS();
};

//---------------------------------------------------------------------------------------
// The value is determined by a data map
class FECORE_API FEMappedValueVec3 : public FEVec3dValuator
{
public:
	FEMappedValueVec3(FEModel* fem);

	void setDataMap(FEDataMap* val, vec3d scl = vec3d(1, 1, 1));

	vec3d operator()(const FEMaterialPoint& pt) override;

	FEVec3dValuator* copy() override;

	vec3d& GetScale() { return m_scale; }

	// get the const value (returns 0 if param is not const)
	vec3d* constValue() override { return &m_scale; }

	void Serialize(DumpStream& ar);

private:
	vec3d			m_scale;
	FEDataMap*		m_val;
};


//-----------------------------------------------------------------------------
// This class calculates a vector based on local element node numbering
class FECORE_API FELocalVectorGenerator : public FEVec3dValuator
{
public:
	FELocalVectorGenerator(FEModel* fem);

	bool Init() override;

	vec3d operator () (const FEMaterialPoint& mp) override;

	FEVec3dValuator* copy() override;

protected:
	int	m_n[2];

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FECORE_API FECylindricalVectorGenerator : public FEVec3dValuator
{
public:
	FECylindricalVectorGenerator(FEModel* fem);

	bool Init() override;

	vec3d operator () (const FEMaterialPoint& mp) override;

	FEVec3dValuator* copy() override;

protected:
	vec3d	m_center;		// center of map
	vec3d	m_axis;			// cylinder axis
	vec3d	m_vector;		// reference direction

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FECORE_API FESphericalVectorGenerator : public FEVec3dValuator
{
public:
	FESphericalVectorGenerator(FEModel* fem);

	bool Init() override;

	vec3d operator () (const FEMaterialPoint& mp) override;

	FEVec3dValuator* copy() override;

protected:
	vec3d	m_center;
	vec3d	m_vector;

	DECLARE_FECORE_CLASS();
};
