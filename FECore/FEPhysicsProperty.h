#pragma once
#include <vector>
#include <string>
#include "fecore_api.h"
#include "vec3d.h"
#include "FECoreBase.h"

class FEModel;

enum class FEValueType
{
	Invalid,
	Bool,
	Int,
	Double,
	Vec3d
};

struct FEValue
{
	FEValueType type = FEValueType::Invalid;
	union {
		bool b;
		int i;
		double d;
		vec3d v3;
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
		}
		return *this;
	}
};

class FECORE_API FEPhysicsProperty
{
	class Imp; // PIMPL for hiding implementation details

public:
	FEPhysicsProperty(FEModel* fem);
	virtual ~FEPhysicsProperty() {}

	void SetSibling(FECoreBase* pc);

	void SetScriptName(const std::string& scriptName);

	void AddVariable(const std::string& varName, FEValueType type);

	bool Init();

public:
	double Value(const std::vector<FEValue>& vars);
	double Value(const std::vector<double>& vars);

	double DerivValue(const std::vector<FEValue>& vars, int varIndex);
	double DerivValue(const std::vector<double>& vars, int varIndex);

protected:
	std::string m_scriptName;

private:
	Imp& m;
};
