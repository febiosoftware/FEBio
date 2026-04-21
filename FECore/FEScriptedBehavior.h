#pragma once
#include <vector>
#include <string>
#include "fecore_api.h"
#include "vec3d.h"
#include "FECoreBase.h"

class FEMaterialPoint;
class FEModel;

enum class FEValueType
{
	Invalid,
	Bool,
	Int,
	Double,
	Vec3d,
	Mat3d
};

struct FEValue
{
	FEValueType type = FEValueType::Invalid;
	union {
		bool b;
		int i;
		double d;
		vec3d v3;
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
		case FEValueType::Vec3d : v3 = other.v3; break;
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
		case FEValueType::Vec3d: v3 = other.v3; break;
		case FEValueType::Mat3d: m3 = other.m3; break;
		}
		return *this;
	}

	void operator = (bool v) { type = FEValueType::Bool; b = v; }
	void operator = (int v) { type = FEValueType::Int; i = v; }
	void operator = (double v) { type = FEValueType::Double; d = v; }
	void operator = (const vec3d& v) { type = FEValueType::Vec3d; v3 = v; }
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
};

class FECORE_API FEScriptedBehavior
{
	class Imp; // PIMPL for hiding implementation details

public:
	FEScriptedBehavior(FEModel* fem);
	virtual ~FEScriptedBehavior() {}

	void SetSibling(FECoreBase* pc);

	void SetScriptName(const std::string& scriptName);

	virtual ScriptContext GetScriptContext() const = 0;

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

FECORE_API std::vector<ScriptInputVariable> GetScriptInputVariables(const std::string& script, const ScriptContext& context);
