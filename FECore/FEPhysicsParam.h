#pragma once
#include <FECore/FEModelParam.h>
#include "fecore_api.h"

class FEModel;

class FECORE_API FEPhysicsParam : public FEParamDouble
{
	class Imp; // PIMPL for hiding implementation details

public:
	FEPhysicsParam(FEModel* fem);

	bool Init() override;

	void AddVariable(const std::string& varName);

	void operator = (double v) { FEParamDouble::operator=(v); }

public:
	double Value(const FEMaterialPoint& pt, const std::vector<double>& vars);

	double DerivValue(const FEMaterialPoint& pt, const std::vector<double>& vars, int varIndex);

private:
	Imp& m;
};
