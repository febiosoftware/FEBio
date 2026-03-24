#pragma once
#include <vector>
#include <string>

class FEModel;

class FEPhysicsProperty
{
	class Imp; // PIMPL for hiding implementation details

public:
	FEPhysicsProperty(FEModel* fem);

	void SetScriptName(const std::string& scriptName);

	void AddVariable(const std::string& varName);

	bool Init();

public:
	double Value(const std::vector<double>& vars);

	double DerivValue(const std::vector<double>& vars, int varIndex);

private:
	Imp& m;
};
