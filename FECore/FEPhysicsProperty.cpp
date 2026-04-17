#include "FEPhysicsProperty.h"
#include <febcode/vm.h>
#include <febcode/parser.h>
#include <febcode/compiler.h>
#include <febcode/differentiator.h>
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

	struct Variable {
		std::string name;
		FEValueType type;
		bool differentiable = true;
	};

	struct Script {
		std::vector<Variable> vars;
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

		FEValue Run(const std::vector<FEValue>& vars)
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
					case FEValueType::Mat3d:
						{
							const mat3d& m3 = vars[i].m3;
							febcode::mat3 m3_febcode(m3(0,0), m3(0, 1), m3(0, 2), m3(1, 0), m3(1, 1), m3(1, 2), m3(2, 0), m3(2, 1), m3(2, 2));
							vm.setGlobal(g.slot, m3_febcode);
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
				result.v3 = vec3d(v.vec3Value.x, v.vec3Value.y, v.vec3Value.z);
				break;
			case febcode::ValueIndex::MAT3:
				result.type = FEValueType::Mat3d;
				febcode::mat3& m = v.mat3Value;
				result.m3 = mat3d(m.m[0][0], m.m[0][1], m.m[0][2], m.m[1][0], m.m[1][1], m.m[1][2], m.m[2][0], m.m[2][1], m.m[2][2]);
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
			auto diffAST = diff.differentiate(*code.program.ast, varName);

			// figure out the return type of the derivative program based on the type of the variable
			// we're differentiating with respect to and the return type of the original program.
			if (code.program.returnType == nullptr) return false;
			febcode::TypeKind returnTypeKind = code.program.returnType->kind;
			if (returnTypeKind == febcode::TypeKind::Double)
			{
				switch (varType->kind)
				{
				case febcode::TypeKind::Double:
					code.program.returnType = code.program.types.Double();
					break;
				case febcode::TypeKind::Vec3:
					code.program.returnType = code.program.types.Vec3();
					break;
				default:
					feLogErrorEx(code.pc->GetFEModel(), "Unsupported variable type for differentiation: only double variables are supported for differentiation when the return type is double");
					return false;
				}
			}
			else if (returnTypeKind == febcode::TypeKind::Vec3)
			{
				switch (varType->kind)
				{
				case febcode::TypeKind::Double:
					code.program.returnType = code.program.types.Vec3();
					break;
				case febcode::TypeKind::Vec3:
					code.program.returnType = code.program.types.Mat3();
					break;
				default:
					feLogErrorEx(code.pc->GetFEModel(), "Unsupported variable type for differentiation: only double variables are supported for differentiation when the return type is double");
					return false;
				}
			}
			else
			{
				feLogErrorEx(code.pc->GetFEModel(), "Unsupported return type for differentiation: only double and vec3 return types are supported");
				return false;
			}

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

	FEValueType returnType = FEValueType::Double;

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

void FEPhysicsProperty::SetProgramReturnType(FEValueType type)
{
	m.returnType = type;
}

void FEPhysicsProperty::SetScriptName(const std::string& scriptName)
{
	m_scriptName = scriptName;
}

void FEPhysicsProperty::AddVariable(const std::string& varName, FEValueType type, bool differentiable)
{
	m.code.vars.push_back({ varName, type, differentiable });
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
		varNamesPrefixed[i] = "_" + m.code.vars[i].name;
	}

	// add the variables to the globals list
	for (int i = 0; i < nvars; ++i)
	{
		switch (m.code.vars[i].type)
		{
		case FEValueType::Double: m.code.AddGlobalDouble(varNamesPrefixed[i]); break;
		case FEValueType::Vec3d : m.code.AddGlobalVec3  (varNamesPrefixed[i]); break;
		default:
			feLogErrorEx(m.fem, "Unsupported variable type for variable \"%s\"", m.code.vars[i].name.c_str());
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
			switch (m.code.vars[j].type)
			{
			case FEValueType::Double: m.valDeriv[i].code.AddGlobalDouble(varNamesPrefixed[j]); break;
			case FEValueType::Vec3d: m.valDeriv[i].code.AddGlobalVec3(varNamesPrefixed[j]); break;
			}
		}
	}

	// compile the main program
	m.code.SetReturnType(m.returnType);
	if (m.code.Init() == false) return false;

	// initialize the derivates
	for (int i = 0; i < m.valDeriv.size(); ++i)
	{
		Imp::Derive& deriv_i = m.valDeriv[i];

		// Set the return type of the original code. 
		// The derivative code's return type will be set during differentiation based on the type of the variable 
		// we're differentiating with respect to.
		deriv_i.code.SetReturnType(m.returnType); 

		// if the corresponding variable is not differentiable, we mark the derivative program as null and skip initialization
		if (!m.code.vars[i].differentiable)
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

FEValue FEPhysicsProperty::Value(const std::vector<FEValue>& vars)
{
	return m.code.Run(vars);
}

double FEPhysicsProperty::Value(const std::vector<double>& vars)
{
	return m.code.Run(vars);
}

FEValue FEPhysicsProperty::DerivValue(const std::vector<FEValue>& vars, int varIndex)
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

		FEValue dr = deriv_i.code.Run(vars);
		return dr;
	}
	return FEValue();
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
