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

class FECORE_API FEScriptedBehavior
{
	class Imp; // PIMPL for hiding implementation details

public:
	FEScriptedBehavior(FEModel* fem);
	virtual ~FEScriptedBehavior() {}

	void SetSibling(FECoreBase* pc);

	void SetScriptName(const std::string& scriptName);

	void SetProgramReturnType(FEValueType type);

	void AddVariable(const std::string& varName, FEValueType type, bool differentiable = true);

	bool Init();

	bool HasDerivative(int id) const;

public:
	FEValue Value(const FEMaterialPoint& mp, const std::vector<FEValue>& vars);
	double Value(const FEMaterialPoint& mp, const std::vector<double>& vars);

	FEValue DerivValue(const FEMaterialPoint& mp, const std::vector<FEValue>& vars, int varIndex);
	double DerivValue(const FEMaterialPoint& mp, const std::vector<double>& vars, int varIndex);

protected:
	std::string m_scriptName;

private:
	Imp& m;
};

struct FECORE_API ScriptContext
{
	FEValueType returnType;
	std::vector<std::pair<std::string, FEValueType>> variables;

	void addVariable(const std::string& name, FEValueType type)
	{
		variables.push_back({ name, type });
	}
};

// helper function to see if a script compiles. This is used in FEBio studio.
FECORE_API bool ValidateScript(const std::string& script, const ScriptContext& context, std::string& err);
