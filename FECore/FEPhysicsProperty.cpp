#include "FEPhysicsProperty.h"
#include <febcode/vm.h>
#include <febcode/parser.h>
#include <febcode/module_vec3.h>
#include <febcode/module_vec2.h>
#include <febcode/module_math.h>
#include <febcode/compiler.h>
#include <febcode/differentiator.h>
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <FECore/FECoreClass.h>
#include <sstream>

class FEPhysicsProperty::Imp
{
public:
	struct Global {
		int slot = -1;
		std::string name;
	};

	struct Script {
		std::vector<std::pair<std::string, FEValueType>> vars;
		std::vector<Global> globals;
		std::string script;
		febcode::Program program;
		FECoreBase* pc = nullptr;

		int AddGlobalDouble(const std::string& name)
		{
			Global g;
			g.name = name;
			g.slot = program.injectGlobal(name, program.types.getDoubleType());
			globals.push_back(g);
			return g.slot;
		}

		int AddGlobalVec3(const std::string& name)
		{
			Global g;
			g.name = name;
			g.slot = program.injectGlobal(name, program.types.getVec3Type());
			globals.push_back(g);
			return g.slot;
		}

		bool Init()
		{
			febcode::ParseSource(program, script);

			febcode::Compiler compiler(program);
			compiler.compile();

			// see if the globals were actually used
			for (auto& g : globals)
			{
				if (program.globals[g.slot].refcount == 0)
				{
					g.slot = -1; // mark as unused
				}
			}
			return true;
		}

		double Run(const std::vector<double>& vars)
		{
			thread_local febcode::VM vm;
			vm.setProgram(program);

			assert(globals.size() == vars.size());
			for (int i = 0; i < (int)vars.size(); ++i)
			{
				Imp::Global& g = globals[i];
				if (g.slot >= 0)
				{
					vm.setGlobal(g.slot, vars[i]);
				}
			}

			if (pc)
			{
				FEParameterList& pl = pc->GetParameterList();
				auto it = pl.first();
				for (int i=0; i < pl.Parameters(); ++i, it++)
				{
					FEParam& pi = *it;
					if (pi.GetFlags() & FEParamFlag::FE_PARAM_USER)
					{
						vm.setInput(pi.name(), pi.value<double>());
					}
				}
			}

			febcode::Value v = vm.run();
			return febcode::getDouble(v);
		}

		double Run(const std::vector<FEValue>& vars)
		{
			thread_local febcode::VM vm;
			vm.setProgram(program);

			assert(globals.size() == vars.size());
			for (int i = 0; i < (int)vars.size(); ++i)
			{
				Imp::Global& g = globals[i];
				if (g.slot >= 0)
				{
					switch (vars[i].type)
					{
					case FEValueType::Bool  : vm.setGlobal(g.slot, (double)vars[i].b); break;
					case FEValueType::Int   : vm.setGlobal(g.slot, (double)vars[i].i); break;
					case FEValueType::Double: vm.setGlobal(g.slot, vars[i].d); break;
					case FEValueType::Vec3d:
						{
							const vec3d& v3 = vars[i].v3;
							febcode::vec3 v3_febcode(v3.x, v3.y, v3.z);
							vm.setGlobal(g.slot, v3_febcode);
							break;
						}
					}
				}
			}

			if (pc)
			{
				FEParameterList& pl = pc->GetParameterList();
				auto it = pl.first();
				for (int i = 0; i < pl.Parameters(); ++i, it++)
				{
					FEParam& pi = *it;
					if (pi.GetFlags() & FEParamFlag::FE_PARAM_USER)
					{
						vm.setInput(pi.name(), pi.value<double>());
					}
				}
			}

			febcode::Value v = vm.run();
			return febcode::getDouble(v);
		}
	};

	struct Derive {
		Script code;
		std::string varName; // variable with respect to which we are taking the derivative
		bool isNullProgram = false; // flag to indicate if the derivative program is null (i.e. the original program does not depend on the variable we're differentiating with respect to)

		bool Init()
		{
			febcode::ParseSource(code.program, code.script);

			febcode::Differentiator diff(code.program);
			auto diffAST = diff.differentiate(*code.program.ast, varName);

			if (!diff.DependencyFound())
			{
				// no dependency was found on the variable we are differentiating with respect to
				isNullProgram = true;
				code.program.ast.reset();
			}
			else
			{
				code.program.ast = std::move(diffAST);

				febcode::Compiler compiler(code.program);
				compiler.compile();

				// see if the globals were actually used
				for (auto& g : code.globals)
				{
					if (code.program.globals[g.slot].refcount == 0)
					{
						g.slot = -1; // mark as unused
					}
				}
			}
			return true;
		}
	};

public:
	FEModel* fem = nullptr;
	FECoreBase* pc = nullptr;

	Script code;
	std::vector<Derive> valDeriv; //!< programs for derivatives
};

FEPhysicsProperty::FEPhysicsProperty(FEModel* fem) : m(*new Imp())
{
	m.fem = fem;
}

void FEPhysicsProperty::SetSibling(FECoreBase* pc)
{
	m.pc = pc;
}

void FEPhysicsProperty::SetScriptName(const std::string& scriptName)
{
	m_scriptName = scriptName;
}

void FEPhysicsProperty::AddVariable(const std::string& varName, FEValueType type)
{
	m.code.vars.push_back({ varName, type });
}

bool FEPhysicsProperty::HasDerivative(int id) const
{
	if ((id < 0) || (id >= m.valDeriv.size())) return false;
	return !m.valDeriv[id].isNullProgram;
}

bool FEPhysicsProperty::Init()
{
#ifndef NDEBUG
	feLogEx(m.fem, "compiling script \"%s\":\n", m_scriptName.c_str());
#endif

	m.code.script = m.fem->GetScript(m_scriptName);
	if (m.code.script.empty())
	{
		feLogErrorEx(m.fem, "Script \"%s\" is empty", m_scriptName.c_str());
		return false;
	}

	m.code.pc = m.pc;

	// get the number of variables
	int nvars = (int)m.code.vars.size();

	// prep the variable names (prefix with underscore)
	std::vector<std::string> varNamesPrefixed(nvars);
	for (int i = 0; i < nvars; ++i)
	{
		varNamesPrefixed[i] = "_" + m.code.vars[i].first;
	}

	// add the variables to the globals list
	for (int i = 0; i < nvars; ++i)
	{
		switch (m.code.vars[i].second)
		{
		case FEValueType::Double: m.code.AddGlobalDouble(varNamesPrefixed[i]); break;
		case FEValueType::Vec3d : m.code.AddGlobalVec3  (varNamesPrefixed[i]); break;
		default:
			feLogErrorEx(m.fem, "Unsupported variable type for variable \"%s\"", m.code.vars[i].first.c_str());
			return false;
		}
	}

	// construct the derivatives
	m.valDeriv.resize(nvars);
	for (int i = 0; i < nvars; ++i)
	{
		m.valDeriv[i].code.script = m.code.script;
		m.valDeriv[i].varName = varNamesPrefixed[i];
		m.valDeriv[i].code.pc = m.pc;

		// add the variables to the derivative code's global list
		for (int j = 0; j < nvars; ++j)
		{
			switch (m.code.vars[j].second)
			{
			case FEValueType::Double: m.valDeriv[i].code.AddGlobalDouble(varNamesPrefixed[j]); break;
			case FEValueType::Vec3d: m.valDeriv[i].code.AddGlobalVec3(varNamesPrefixed[j]); break;
			}
		}

		if (m.code.vars[i].second != FEValueType::Double)
		{
			m.valDeriv[i].isNullProgram = true; // we only support derivatives with respect to double variables for now
		}
	}

	// compile the main program
	if (m.code.Init() == false) return false;

	// initialize the derivates
	for (int i = 0; i < m.valDeriv.size(); ++i)
	{
		Imp::Derive& deriv_i = m.valDeriv[i];
		if (!deriv_i.isNullProgram && deriv_i.Init() == false)
		{
			feLogErrorEx(m.fem, "Failed to initialize derivative valuator for variable %d", i);
			return false;
		}

#ifndef NDEBUG
		if (deriv_i.isNullProgram)
		{
			feLogEx(m.fem, "Derivative AST w.r.t %s is null (no dependency found)\n\n", deriv_i.varName.c_str());
		}
		else
		{
			feLogEx(m.fem, "Derivative AST w.r.t %s :\n>>>\n", deriv_i.varName.c_str());
			std::stringstream ss;
			febcode::prettyPrintAST(ss, *deriv_i.code.program.ast);
			feLogEx(m.fem, "%s\n<<<\n\n", ss.str().c_str());
		}
#endif
	}
	return true;
}

double FEPhysicsProperty::Value(const std::vector<FEValue>& vars)
{
	return m.code.Run(vars);
}

double FEPhysicsProperty::Value(const std::vector<double>& vars)
{
	return m.code.Run(vars);
}

double FEPhysicsProperty::DerivValue(const std::vector<FEValue>& vars, int varIndex)
{
	if ((varIndex >= 0) && (varIndex < m.valDeriv.size()))
	{
		Imp::Derive& deriv_i = m.valDeriv[varIndex];
		if (deriv_i.isNullProgram)
		{
			// this derivative was optimized out, so we return zero
			return 0.0;
		}

		double dr = deriv_i.code.Run(vars);
		return dr;
	}
	return 0.0;
}

double FEPhysicsProperty::DerivValue(const std::vector<double>& vars, int varIndex)
{
	if ((varIndex >= 0) && (varIndex < m.valDeriv.size()))
	{
		Imp::Derive& deriv_i = m.valDeriv[varIndex];
		if (deriv_i.isNullProgram)
		{
			// this derivative was optimized out, so we return zero
			return 0.0;
		}

		double dr = deriv_i.code.Run(vars);
		return dr;
	}
	return 0.0;
}
