#pragma once
#include <vector>
#include <string>
#include "fecore_api.h"
#include "vec3d.h"
#include "FECoreBase.h"

class FEMaterialPoint;
class FEModel;

// NOTE: Don't change the order of these values!
// Also make sure this list matches the union in FEValue!
enum class FEValueType
{
	Invalid,
	Bool,
	Int,
	Double,
	Vec2d,
	Vec3d,
	Mat2d,
	Mat3d
};

struct FEValue
{
	FEValueType type = FEValueType::Invalid;
	union {
		bool b;
		int i;
		double d;
		vec2d v2;
		vec3d v3;
		mat2d m2;
		mat3d m3;
	};

	FEValue() {}

	FEValue(const FEValue& other)
	{
		type = other.type;
		switch (type)
		{
		case FEValueType::Bool  : b = other.b; break;
		case FEValueType::Int   : i = other.i; break;
		case FEValueType::Double: d = other.d; break;
		case FEValueType::Vec2d : v2 = other.v2; break;
		case FEValueType::Vec3d : v3 = other.v3; break;
		case FEValueType::Mat2d : m2 = other.m2; break;
		case FEValueType::Mat3d : m3 = other.m3; break;
		}
	}

	FEValue& operator = (const FEValue& other)
	{
		type = other.type;
		switch (type)
		{
		case FEValueType::Bool: b = other.b; break;
		case FEValueType::Int: i = other.i; break;
		case FEValueType::Double: d = other.d; break;
		case FEValueType::Vec2d: v2 = other.v2; break;
		case FEValueType::Vec3d: v3 = other.v3; break;
		case FEValueType::Mat2d: m2 = other.m2; break;
		case FEValueType::Mat3d: m3 = other.m3; break;
		}
		return *this;
	}

	void operator = (bool v) { type = FEValueType::Bool; b = v; }
	void operator = (int v) { type = FEValueType::Int; i = v; }
	void operator = (double v) { type = FEValueType::Double; d = v; }
	void operator = (const vec2d& v) { type = FEValueType::Vec2d; v2 = v; }
	void operator = (const vec3d& v) { type = FEValueType::Vec3d; v3 = v; }
	void operator = (const mat2d& v) { type = FEValueType::Mat2d; m2 = v; }
	void operator = (const mat3d& v) { type = FEValueType::Mat3d; m3 = v; }
};

struct FECORE_API ScriptContext
{
	struct Variable
	{
		std::string name;
		FEValueType type;
		bool differentiable;
	};

	FEValueType returnType;
	std::vector<Variable> variables;

	void addVariable(const std::string& name, FEValueType type, bool differentiable)
	{
		variables.push_back({ name, type, differentiable });
	}

	bool operator == (const ScriptContext& other) const
	{
		if (returnType != other.returnType) return false;
		if (variables.size() != other.variables.size()) return false;
		for (size_t i = 0; i < variables.size(); ++i)
		{
			if (variables[i].name != other.variables[i].name) return false;
			if (variables[i].type != other.variables[i].type) return false;
			if (variables[i].differentiable != other.variables[i].differentiable) return false;
		}
		return true;
	}

	bool operator != (const ScriptContext& other) const
	{
		return !(*this == other);
	}
};

class FECORE_API FEScriptedBehavior : public FECoreBase
{
	class Imp; // PIMPL for hiding implementation details

	FECORE_SUPER_CLASS(FESCRIPT_ID)
	FECORE_BASE_CLASS(FEScriptedBehavior)

public:
	FEScriptedBehavior(FEModel* fem);
	virtual ~FEScriptedBehavior() {}

	void SetScriptName(const std::string& scriptName);

	void SetScriptContext(const ScriptContext& context);

	ScriptContext GetScriptContext() const;

	bool Init();

	bool HasDerivative(int id) const;

public:
	FEValue Value(const FEMaterialPoint& mp, const std::vector<FEValue>& vars);
	double Value(const FEMaterialPoint& mp, const std::vector<double>& vars);

	FEValue DerivValue(const FEMaterialPoint& mp, const std::vector<FEValue>& vars, int varIndex);
	double DerivValue(const FEMaterialPoint& mp, const std::vector<double>& vars, int varIndex);

private:
	std::string m_scriptName;

private:
	Imp& m;
};

// helper function to see if a script compiles. This is used in FEBio studio.
FECORE_API bool ValidateScript(const std::string& script, const ScriptContext& context, std::string& err);

struct FECORE_API ScriptInputVariable
{
	std::string name;
	FEValueType type;
};

FECORE_API std::vector<ScriptInputVariable> GetScriptInputVariables(const std::string& script, const ScriptContext& context, bool& ok);

// use this as a base class for scripted model components.
template <class Base> 
class FEScripted : public Base
{
public:
	FEScripted(FEModel* fem) : Base(fem), m_script(fem) {}

	void BuildParamList() override {
		Base::BuildParamList();

		ADD_PROPERTY(m_script, "script");
	}

	void SetScriptContext(const ScriptContext& context)
	{
		m_script.SetScriptContext(context);
	}

	FEValue Value(const FEMaterialPoint& mp, const std::vector<FEValue>& vars)
	{
		return m_script.Value(mp, vars);
	}

	FEValue DerivValue(const FEMaterialPoint& mp, const std::vector<FEValue>& vars, int varIndex)
	{
		return m_script.DerivValue(mp, vars, varIndex);
	}

	bool HasDerivative(int varIndex) const
	{
		return m_script.HasDerivative(varIndex);
	}

private:
	FEScriptedBehavior m_script;
};
