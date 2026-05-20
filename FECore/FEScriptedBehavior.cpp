#include "FEScriptedBehavior.h"
#include <febcode/vm.h>
#include <febcode/parser.h>
#include <febcode/compiler.h>
#include <febcode/differentiator.h>
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <FECore/FEModelParam.h>
#include <sstream>

class FEScriptedBehavior::Imp
{
public:
	struct Global {
		int slot = -1;
		std::string name;
	};

	struct Script {
		std::vector<Global> globals;
		std::string script;
		febcode::Program program;
		FECoreBase* pc = nullptr;

		int AddGlobalDouble(const std::string& name)
		{
			Global g;
			g.name = name;
			g.slot = program.injectGlobal(name, program.types.Double());
			globals.push_back(g);
			return g.slot;
		}

		int AddGlobalVec3(const std::string& name)
		{
			Global g;
			g.name = name;
			g.slot = program.injectGlobal(name, program.types.Vec3());
			globals.push_back(g);
			return g.slot;
		}

		void SetReturnType(FEValueType type)
		{
			switch (type)
			{
			case FEValueType::Bool  : program.returnType = program.types.Bool(); break;
			case FEValueType::Int   : program.returnType = program.types.Int(); break;
			case FEValueType::Double: program.returnType = program.types.Double(); break;
			case FEValueType::Vec3d : program.returnType = program.types.Vec3(); break;
			case FEValueType::Mat3d : program.returnType = program.types.Mat3(); break;
			default:
				assert(false);
			}
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

		double Run(const FEMaterialPoint& mp, const std::vector<double>& vars)
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
					switch (pi.type())
					{
					case FEParamType::FE_PARAM_BOOL: vm.setInput(pi.name(), pi.value<bool  >()); break;
					case FEParamType::FE_PARAM_INT: vm.setInput(pi.name(), pi.value<int   >()); break;
					case FEParamType::FE_PARAM_DOUBLE: vm.setInput(pi.name(), pi.value<double>()); break;
					case FEParamType::FE_PARAM_DOUBLE_MAPPED:
					{
						FEParamDouble& p = pi.value<FEParamDouble>();
						vm.setInput(pi.name(), p(mp));
						break;
					}
					default:
						assert(false);
					}
				}
			}

			febcode::Value v = vm.run();
			return febcode::getDouble(v);
		}

		FEValue Run(const FEMaterialPoint& mp, const std::vector<FEValue>& vars)
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
					case FEValueType::Vec3d : vm.setGlobal(g.slot, vars[i].v3); break;
					case FEValueType::Mat3d : vm.setGlobal(g.slot, vars[i].m3); break;
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
						switch (pi.type())
						{
						case FEParamType::FE_PARAM_BOOL  : vm.setInput(pi.name(), pi.value<bool  >()); break;
						case FEParamType::FE_PARAM_INT   : vm.setInput(pi.name(), pi.value<int   >()); break;
						case FEParamType::FE_PARAM_DOUBLE: vm.setInput(pi.name(), pi.value<double>()); break;
						case FEParamType::FE_PARAM_VEC3D : vm.setInput(pi.name(), pi.value<vec3d >()); break;
						case FEParamType::FE_PARAM_DOUBLE_MAPPED:
						{
							FEParamDouble& p = pi.value<FEParamDouble>();
							vm.setInput(pi.name(), p(mp));
							break;
						}
						case FEParamType::FE_PARAM_VEC3D_MAPPED:
						{
							FEParamVec3& p = pi.value<FEParamVec3>();
							vm.setInput(pi.name(), p(mp));
							break;
						}
						default:
							assert(false);
						}
					}
				}
			}

			febcode::Value v = vm.run();

			FEValue result;
			switch (v.index)
			{
			case febcode::ValueIndex::BOOL:
				result.type = FEValueType::Bool;
				result.b = v.b;
				break;
			case febcode::ValueIndex::INT:
				result.type = FEValueType::Int;
				result.i = v.i;
				break;
			case febcode::ValueIndex::DOUBLE:
				result.type = FEValueType::Double;
				result.d = v.d;
				break;
			case febcode::ValueIndex::VEC3:
				result.type = FEValueType::Vec3d;
				result.v3 = v.vec3Value;
				break;
			case febcode::ValueIndex::MAT3:
				result.type = FEValueType::Mat3d;
				febcode::mat3& m = v.mat3Value;
				result.m3 = m;
				break;
			}
			return result;
		}
	};

	struct Derive {
		Script code;
		std::string varName; // variable with respect to which we are taking the derivative
		bool isNullProgram = false; // flag to indicate if the derivative program is null (i.e. the original program does not depend on the variable we're differentiating with respect to)

		bool Init()
		{
			febcode::ParseSource(code.program, code.script);

			// get the type of the variable
			if (code.program.globalIndices.count(varName) == 0)
			{
				feLogErrorEx(code.pc->GetFEModel(), "Internal error: variable \"%s\" not found in global indices after differentiation", varName.c_str());
				return false;
			}
			febcode::Type varType = code.program.globalType(varName);
			if (varType == nullptr) return false;

			febcode::Differentiator diff(code.program);
			diff.differentiate(varName);

			if (!diff.DependencyFound())
			{
				// no dependency was found on the variable we are differentiating with respect to
				isNullProgram = true;
				code.program.ast.reset();
			}
			else
			{
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

	ScriptContext context;

	Script code;
	std::vector<Derive> valDeriv; //!< programs for derivatives
};

FEScriptedBehavior::FEScriptedBehavior(FEModel* fem) : FECoreBase(fem), m(*new Imp())
{
	m.fem = fem;
}

void FEScriptedBehavior::SetScriptContext(const ScriptContext& context)
{
	m.context = context;
}

void FEScriptedBehavior::SetDefaultContext(FEValueType returnType, FECoreBase* parent)
{
	ScriptContext ctx;
	ctx.returnType = returnType;
	ctx.addVariable("pos0", FEValueType::Vec3d, false);
	ctx.addVariable("time", FEValueType::Double, false);
	SetScriptContext(ctx);
}

ScriptContext FEScriptedBehavior::GetScriptContext() const
{
	return m.context;
}

void FEScriptedBehavior::SetScriptName(const std::string& scriptName)
{
	m_scriptName = scriptName;
}

bool FEScriptedBehavior::HasDerivative(int id) const
{
	if ((id < 0) || (id >= m.valDeriv.size())) return false;
	return !m.valDeriv[id].isNullProgram;
}

bool FEScriptedBehavior::Init()
{
	if (m_scriptName.empty())
		m_scriptName = GetName();

#ifndef NDEBUG
	feLogEx(m.fem, "compiling script \"%s\":\n", m_scriptName.c_str());
#endif

	m.code.script = m.fem->GetScript(m_scriptName).script;
	if (m.code.script.empty())
	{
		feLogErrorEx(m.fem, "Script \"%s\" is empty", m_scriptName.c_str());
		return false;
	}

	ScriptContext ctx = m.context;

	m.code.SetReturnType(ctx.returnType);
	m.code.pc = this;

	// get the number of variables
	int nvars = (int)ctx.variables.size();

	// prep the variable names (prefix with underscore)
	std::vector<std::string> varNamesPrefixed(nvars);
	for (int i = 0; i < nvars; ++i)
	{
		varNamesPrefixed[i] = "_" + ctx.variables[i].name;
	}

	// add the variables to the globals list
	for (int i = 0; i < nvars; ++i)
	{
		switch (ctx.variables[i].type)
		{
		case FEValueType::Double: m.code.AddGlobalDouble(varNamesPrefixed[i]); break;
		case FEValueType::Vec3d : m.code.AddGlobalVec3  (varNamesPrefixed[i]); break;
		default:
			feLogErrorEx(m.fem, "Unsupported variable type for variable \"%s\"", ctx.variables[i].name.c_str());
			return false;
		}
	}

	// construct the derivatives
	m.valDeriv.resize(nvars);
	for (int i = 0; i < nvars; ++i)
	{
		m.valDeriv[i].code.script = m.code.script;
		m.valDeriv[i].varName = varNamesPrefixed[i];
		m.valDeriv[i].code.pc = m.code.pc;

		// add the variables to the derivative code's global list
		for (int j = 0; j < nvars; ++j)
		{
			switch (ctx.variables[j].type)
			{
			case FEValueType::Double: m.valDeriv[i].code.AddGlobalDouble(varNamesPrefixed[j]); break;
			case FEValueType::Vec3d: m.valDeriv[i].code.AddGlobalVec3(varNamesPrefixed[j]); break;
			}
		}
	}

	// compile the main program
	m.code.SetReturnType(ctx.returnType);
	if (m.code.Init() == false) return false;

	// initialize the derivates
	for (int i = 0; i < m.valDeriv.size(); ++i)
	{
		Imp::Derive& deriv_i = m.valDeriv[i];

		// Set the return type of the original code. 
		// The derivative code's return type will be set during differentiation based on the type of the variable 
		// we're differentiating with respect to.
		deriv_i.code.SetReturnType(ctx.returnType); 

		// if the corresponding variable is not differentiable, we mark the derivative program as null and skip initialization
		if (!ctx.variables[i].differentiable)
		{
			deriv_i.isNullProgram = true;
			continue;
		}

		bool success = true;
		try {

			if (!deriv_i.isNullProgram && deriv_i.Init() == false)
			{
				feLogErrorEx(m.fem, "Failed to initialize derivative valuator for variable %d", i);
				success = false;
			}
		}
		catch (std::runtime_error e)
		{
			feLogErrorEx(m.fem, "Runtime error while initializing derivative valuator for variable %d: %s", i, e.what());
			success = false;
		}
		catch (...)
		{
			feLogErrorEx(m.fem, "Unknown error while initializing derivative valuator for variable %d", i);
			success = false;
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
		if (!success) return false;
	}
	return true;
}

double FEScriptedBehavior::Value(double v)
{
	thread_local std::vector<double> vars(1);
	vars[0] = v;

	static FEMaterialPoint dummy;

	return Value(dummy, vars);
}

FEValue FEScriptedBehavior::Value(const FEMaterialPoint& mp)
{
	thread_local std::vector<FEValue> vars(2);
	vars[0] = mp.m_r0; // pos0
	vars[1] = GetFEModel()->GetCurrentTime(); // time
	return Value(mp, vars);
}

FEValue FEScriptedBehavior::Value(const FEMaterialPoint& mp, const std::vector<FEValue>& vars)
{
	return m.code.Run(mp, vars);
}

double FEScriptedBehavior::Value(const FEMaterialPoint& mp, const std::vector<double>& vars)
{
	return m.code.Run(mp, vars);
}

FEValue FEScriptedBehavior::DerivValue(const FEMaterialPoint& mp, const std::vector<FEValue>& vars, int varIndex)
{
	if ((varIndex >= 0) && (varIndex < m.valDeriv.size()))
	{
		Imp::Derive& deriv_i = m.valDeriv[varIndex];
		if (deriv_i.isNullProgram)
		{
			// this derivative was optimized out, so we return zero
			switch (deriv_i.code.program.returnType->kind)
			{
			case febcode::TypeKind::Double:
			{
				FEValue result;
				result.type = FEValueType::Double;
				result.d = 0.0;
				return result;
			}
			case febcode::TypeKind::Vec3:
			{
				FEValue result;
				result.type = FEValueType::Vec3d;
				result.v3 = vec3d(0.0, 0.0, 0.0);
				return result;
			}
			case febcode::TypeKind::Mat3:
			{
				FEValue result;
				result.type = FEValueType::Mat3d;
				result.m3 = mat3dd(0.0);
				return result;
			}
			}
		}

		FEValue dr = deriv_i.code.Run(mp, vars);
		return dr;
	}
	return FEValue();
}

double FEScriptedBehavior::DerivValue(const FEMaterialPoint& mp, const std::vector<double>& vars, int varIndex)
{
	if ((varIndex >= 0) && (varIndex < m.valDeriv.size()))
	{
		Imp::Derive& deriv_i = m.valDeriv[varIndex];
		if (deriv_i.isNullProgram)
		{
			// this derivative was optimized out, so we return zero
			return 0.0;
		}

		double dr = deriv_i.code.Run(mp, vars);
		return dr;
	}
	return 0.0;
}

bool ValidateScript(febcode::Program& program, const std::string& script, const ScriptContext& context, std::string& err)
{
	err.clear();

	switch (context.returnType)
	{
	case FEValueType::Bool  : program.returnType = program.types.Bool(); break;
	case FEValueType::Int   : program.returnType = program.types.Int(); break;
	case FEValueType::Double: program.returnType = program.types.Double(); break;
	case FEValueType::Vec3d : program.returnType = program.types.Vec3(); break;
	case FEValueType::Mat3d : program.returnType = program.types.Mat3(); break;
	default:
		err = "Invalid return type specified in script context";
		return false;
	}

	try {

		std::string varPrefix = "_";
		for (const auto& var : context.variables)
		{
			switch (var.type)
			{
			case FEValueType::Bool  : program.injectGlobal(varPrefix + var.name, program.types.Bool()); break;
			case FEValueType::Int   : program.injectGlobal(varPrefix + var.name, program.types.Int()); break;
			case FEValueType::Double: program.injectGlobal(varPrefix + var.name, program.types.Double()); break;
			case FEValueType::Vec3d : program.injectGlobal(varPrefix + var.name, program.types.Vec3()); break;
			case FEValueType::Mat3d : program.injectGlobal(varPrefix + var.name, program.types.Mat3()); break;
			default:
				err = "Unsupported variable type for variable \"" + var.name + "\" in script context";
				return false;
			}
		}

		febcode::ParseSource(program, script);

		febcode::Compiler compiler(program);

		compiler.compile();
	}
	catch (const std::exception& e)
	{
		err = e.what();
		return false;
	}
	catch (...)
	{
		err = "Unknown error compiling code";
		return false;
	}
	return true;
}

bool ValidateScript(const std::string& script, const ScriptContext& context, std::string& err)
{
	febcode::Program program;
	return ValidateScript(program, script, context, err);
}

FECORE_API std::vector<ScriptInputVariable> GetScriptInputVariables(const std::string& script, const ScriptContext& context, bool& ok)
{
	std::vector<ScriptInputVariable> inputs;

	febcode::Program program;
	std::string err;
	bool success = ValidateScript(program, script, context, /*out*/ err);
	ok = success;
	if (success)
	{
		for (auto& input : program.inputs)
		{
			FEValueType valueType;
			switch (input.type->kind)
			{
			case febcode::TypeKind::Bool  : valueType = FEValueType::Bool  ; break;
			case febcode::TypeKind::Int   : valueType = FEValueType::Int   ; break;
			case febcode::TypeKind::Double: valueType = FEValueType::Double; break;
			case febcode::TypeKind::Vec3  : valueType = FEValueType::Vec3d ; break;
			case febcode::TypeKind::Mat3  : valueType = FEValueType::Mat3d ; break;
			default:
				assert(false);
				continue;
			}
			inputs.push_back({input.name, valueType});
		}
	}
	return inputs;
}
