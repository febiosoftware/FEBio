#include "FEPhysicsParam.h"
#include <FECore/FECodeValuator.h>
#include <FECore/log.h>

class FEPhysicsParam::Imp
{
public:
	struct Deriv {
		std::unique_ptr<FECodeValuator> derivVal;
		std::vector<int> slots;
	};

public:
	FEModel* fem = nullptr;

	FECodeValuator* code = nullptr;
	std::vector<std::string> varNames;
	std::vector<int> slots;

	std::vector<Deriv> valDeriv; //!< valuators for derivatives
};

FEPhysicsParam::FEPhysicsParam(FEModel* fem) : m(*new Imp())
{
	m.fem = fem;
}

void FEPhysicsParam::AddVariable(const std::string& varName)
{
	m.varNames.push_back(varName);
}

bool FEPhysicsParam::Init()
{
	m.code = dynamic_cast<FECodeValuator*>(valuator());
	if (m.code)
	{
		std::string scriptName = m.code->GetScriptName();

		// get the number of variables
		int nvars = (int)m.varNames.size();

		// prep the variable names (prefix with underscore)
		std::vector<std::string> varNamesPrefixed(nvars);
		for (int i = 0; i < nvars; ++i)
		{
			varNamesPrefixed[i] = "_" + m.varNames[i];
		}

		// add the variables to the valuator's global list
		m.slots.resize(nvars, -1);
		for (int i=0; i<nvars; ++i)
		{
			m.slots[i] = m.code->AddGlobalDouble(varNamesPrefixed[i]);
		}

		// construct the derivatives
		for (int i = 0; i < nvars; ++i)
		{
			m.valDeriv.emplace_back();
			m.valDeriv[i].derivVal.reset(fecore_new_class<FECodeValuator>("FECodeValuator", m.fem));
			m.valDeriv[i].slots = std::vector<int>(nvars, -1);
			m.valDeriv[i].derivVal->SetScriptName(scriptName);
			m.valDeriv[i].derivVal->CompileDerivative(varNamesPrefixed[i]);

			// add the variables to the derivative valuator's global list
			for (int j = 0; j < nvars; ++j)
			{
				m.valDeriv[i].slots[j] = m.valDeriv[i].derivVal->AddGlobalDouble(varNamesPrefixed[j]);
			}
		}

		// initialize the derivates
		for (int i = 0; i < m.valDeriv.size(); ++i)
		{
			if (m.valDeriv[i].derivVal->Init() == false)
			{
				feLogErrorEx(m.fem, "Failed to initialize derivative valuator for variable %d", i);
				return false;
			}
		}
	}
	return FEParamDouble::Init();
}

double FEPhysicsParam::Value(const FEMaterialPoint& pt, const std::vector<double>& vars)
{
	if (m.code)
	{
		std::vector<std::pair<int, double>> globals(vars.size());
		for (int i = 0; i < m.slots.size(); ++i)
		{
			globals[i].first = m.slots[i];
			globals[i].second = vars[i];
		}

		double r = m.code->run(pt, globals);
		return m_scl*r;
	}
	else
		return operator()(pt);
}

double FEPhysicsParam::DerivValue(const FEMaterialPoint& pt, const std::vector<double>& vars, int varIndex)
{
	if ((varIndex >= 0) && (varIndex < m.valDeriv.size()))
	{
		Imp::Deriv& deriv_i = m.valDeriv[varIndex];

		std::vector<std::pair<int, double>> globals(deriv_i.slots.size());
		for (int i = 0; i < deriv_i.slots.size(); ++i)
		{
			globals[i].first = deriv_i.slots[i];
			globals[i].second = vars[i];
		}

		double dr = deriv_i.derivVal->run(pt, globals);
		return m_scl*dr;
	}
	return 0.0;
}
