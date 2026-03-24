#include "FEPhysicsProperty.h"
#include <febcode/vm.h>
#include <febcode/parser.h>
#include <febcode/module_vec3.h>
#include <febcode/module_vec2.h>
#include <febcode/module_math.h>
#include <febcode/compiler.h>
#include <febcode/differentiator.h>
#include <FECore/FEMaterialPoint.h>
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <sstream>

class FEPhysicsProperty::Imp
{
public:
	struct Global {
		int slot = -1;
		std::string name;
	};

	struct Script {
		std::vector<std::string> varNames;
		std::vector<Global> globals;
		std::string script;
		febcode::Program program;

		int AddGlobalDouble(const std::string& name)
		{
			Global g;
			g.name = name;
			g.slot = program.injectGlobal(name, program.types.getDoubleType());
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

		double Run(const std::vector<double>& args)
		{
			thread_local febcode::VM vm;
			vm.setProgram(program);
			assert(args.size() == globals.size());

			for (int i=0; i<(int)args.size(); ++i)
			{
				Global& g = globals[i];
				if (g.slot >= 0)
					vm.setGlobal(g.slot, args[i]);
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
	std::string scriptName;

	Script code;
	std::vector<Derive> valDeriv; //!< programs for derivatives
};

FEPhysicsProperty::FEPhysicsProperty(FEModel* fem) : m(*new Imp())
{
	m.fem = fem;
}

void FEPhysicsProperty::SetScriptName(const std::string& scriptName)
{
	m.scriptName = scriptName;
}

void FEPhysicsProperty::AddVariable(const std::string& varName)
{
	m.code.varNames.push_back(varName);
}

bool FEPhysicsProperty::Init()
{
#ifndef NDEBUG
	feLogEx(m.fem, "compiling script \"%s\":\n", m.scriptName.c_str());
#endif

	m.code.script = m.fem->GetScript(m.scriptName);
	if (m.code.script.empty())
	{
		feLogErrorEx(m.fem, "Script \"%s\" is empty", m.scriptName.c_str());
		return false;
	}

	// get the number of variables
	int nvars = (int)m.code.varNames.size();

	// prep the variable names (prefix with underscore)
	std::vector<std::string> varNamesPrefixed(nvars);
	for (int i = 0; i < nvars; ++i)
	{
		varNamesPrefixed[i] = "_" + m.code.varNames[i];
	}

	// add the variables to the globals list
	for (int i = 0; i < nvars; ++i)
	{
		m.code.AddGlobalDouble(varNamesPrefixed[i]);
	}

	// construct the derivatives
	for (int i = 0; i < nvars; ++i)
	{
		m.valDeriv.emplace_back();
		m.valDeriv[i].code.script = m.code.script;
		m.valDeriv[i].varName = varNamesPrefixed[i];

		// add the variables to the derivative code's global list
		for (int j = 0; j < nvars; ++j)
		{
			m.valDeriv[i].code.AddGlobalDouble(varNamesPrefixed[j]);
		}
	}

	// compile the main program
	if (m.code.Init() == false) return false;

	// initialize the derivates
	for (int i = 0; i < m.valDeriv.size(); ++i)
	{
		Imp::Derive& deriv_i = m.valDeriv[i];
		if (deriv_i.Init() == false)
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

double FEPhysicsProperty::Value(const std::vector<double>& vars)
{
	double r = m.code.Run(vars);
	return r;
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
