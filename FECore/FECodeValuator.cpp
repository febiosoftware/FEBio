/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2026 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/
#include "FECodeValuator.h"
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
#include <FECore/FESurfaceLoad.h>
#include <febcode/types.h>
#include <iostream>
#include <sstream>

class FECodeValuator::Imp {
public:
	struct Global {
		int slot = -1;
		std::string name;
	};

	std::string scriptCode;

	febcode::Program program;
	std::vector<Global> globals;

	bool compileDeriv = false;
	std::string derivVarName;

	bool isNullProgram = false;

	FESurface* surf = nullptr; // for surface valuators

	int AddGlobal(const std::string& name, febcode::TypeKind typeKind)
	{
		// make sure the name is not defined yet
		for (auto& g : globals)
		{
			if (g.name == name)
			{
				return -1;
			}
		}

		febcode::Type type = program.types.getTypeFromKind(typeKind);
		if (type == nullptr)
		{
			return -1;
		}

		// add the new global
		Global g({ -1, name });
		g.slot = program.injectGlobal(g.name, type);
		globals.push_back(g);
		return g.slot;
	}
};

BEGIN_FECORE_CLASS(FECodeValuator, FEScalarValuator)
	ADD_PARAMETER(m_scriptName, "code");
END_FECORE_CLASS()

FECodeValuator::FECodeValuator(FEModel* fem) : FEScalarValuator(fem), m(*new Imp())
{
	febcode::Vec2Module vec2Module;
	vec2Module.Register(m.program);
	febcode::Vec3Module vec3Module;
	vec3Module.Register(m.program);
	febcode::MathModule mathModule;
	mathModule.Register(m.program);

	m.AddGlobal("_pos0", febcode::TypeKind::Vec3);
	m.AddGlobal("_time", febcode::TypeKind::Double);
	m.AddGlobal("_norm0", febcode::TypeKind::Vec3);
}

FECodeValuator::~FECodeValuator()
{

}

bool FECodeValuator::IsNullProgram() const
{
	return m.isNullProgram;
}

void FECodeValuator::CompileDerivative(const std::string& varName)
{
	m.compileDeriv = true;
	m.derivVarName = varName;
}

bool FECodeValuator::CompileScript()
{
	if (m.scriptCode.empty()) return false;

	FECoreBase* parent = GetParent();
	try {
		febcode::Type vec3_t = m.program.types.getVec3Type();
		if (!vec3_t) {
			feLogError("Error compiling code: 'vec3' type not defined");
			return false;
		}

		febcode::ParseSource(m.program, m.scriptCode);

		if (m.compileDeriv)
		{
			febcode::Differentiator diff(m.program);
			auto diffAST = diff.differentiate(*m.program.ast, m.derivVarName);

			if (!diff.DependencyFound())
			{
				// no dependency was found on the variable we are differentiating with respect to
				m.isNullProgram = true;
			}

			m.program.ast = std::move(diffAST);

#ifndef NDEBUG
			feLog("Derivative AST w.r.t %s :\n>>>\n", m.derivVarName.c_str());
			std::stringstream ss;
			febcode::prettyPrintAST(ss, *m.program.ast);
			feLog("%s\n<<<\n\n", ss.str().c_str());
#endif
		}

		febcode::Compiler compiler(m.program);

		compiler.compile();

		// see if the globals were actually used
		for (auto& g : m.globals)
		{
			if (m.program.globals[g.slot].refcount == 0)
			{
				g.slot = -1; // mark as unused
			}
		}
	}
	catch (const std::exception& e)
	{
		feLogError("Error compiling code: %s", e.what());
		return false;
	}
	catch (...)
	{
		feLogError("Unknown error compiling code");
		return false;
	}

	return true;
}

int FECodeValuator::AddGlobalDouble(const std::string& name)
{
	return m.AddGlobal(name, febcode::TypeKind::Double);
}

bool FECodeValuator::Init()
{
	if (FEScalarValuator::Init() == false) return false;

	FECoreBase* parent = GetParent();
	if (auto s = dynamic_cast<FESurfaceLoad*>(parent))
	{
		m.surf = &s->GetSurface();
	}

	// get the script
	m.scriptCode = GetFEModel()->GetScript(m_scriptName);
	if (m.scriptCode.empty())
	{
		feLogError("Script \"%s\" not found or is empty.", m_scriptName.c_str());
		return false;
	}

	// compile the script
	if (CompileScript() == false) return false;

	return true;
}

double FECodeValuator::operator()(const FEMaterialPoint& pt)
{
	thread_local febcode::VM vm;
	vm.setProgram(m.program);
	if (m.globals[0].slot >= 0) vm.setGlobal(m.globals[0].slot, febcode::vec3(pt.m_r0.x, pt.m_r0.y, pt.m_r0.z));
	if (m.globals[1].slot >= 0) vm.setGlobal(m.globals[1].slot, GetFEModel()->GetTime().currentTime);
	if (m.surf && (m.globals[2].slot >= 0))
	{
		FESurfaceElement* el = dynamic_cast<FESurfaceElement*>(pt.m_elem);
		if (el)
		{
			vec3d n = m.surf->SurfaceNormal(*el, pt.m_index);
			vm.setGlobal(m.globals[2].slot, febcode::vec3(n.x, n.y, n.z));
		}
	}
	febcode::Value v = vm.run();
	return febcode::getDouble(v);
}

double FECodeValuator::run(const FEMaterialPoint& pt, const std::vector<int>& slots, const std::vector<double>& values) const
{
	thread_local febcode::VM vm;
	vm.setProgram(m.program);
	if (m.globals[0].slot >= 0) vm.setGlobal(m.globals[0].slot, febcode::vec3(pt.m_r0.x, pt.m_r0.y, pt.m_r0.z));
	if (m.globals[1].slot >= 0) vm.setGlobal(m.globals[1].slot, GetFEModel()->GetTime().currentTime);
	if (m.surf && (m.globals[2].slot >= 0))
	{
		FESurfaceElement* el = dynamic_cast<FESurfaceElement*>(pt.m_elem);
		if (el)
		{
			vec3d n = m.surf->SurfaceNormal(*el, pt.m_index);
			vm.setGlobal(m.globals[2].slot, febcode::vec3(n.x, n.y, n.z));
		}
	}
	for (int i = 0; i < slots.size(); ++i)
	{
		if (slots[i] >= 0)
			vm.setGlobal(slots[i], values[i]);
	}
	febcode::Value v = vm.run();
	return febcode::getDouble(v);
}

FEScalarValuator* FECodeValuator::copy()
{
	FECodeValuator* p = new FECodeValuator(GetFEModel());
	p->m_scriptName = m_scriptName;
	p->m.scriptCode = m.scriptCode;
	return p;
}

void FECodeValuator::Serialize(DumpStream& ar)
{
	FEScalarValuator::Serialize(ar);
	if (!ar.IsShallow())
	{
		ar& m.scriptCode;

		if (!ar.IsSaving())
			CompileScript();
	}
}

bool ValidateScript(const std::string& script, std::string& err)
{
	err.clear();
	febcode::Program program;
	try {
		febcode::Vec2Module vec2Module;
		vec2Module.Register(program);
		febcode::Vec3Module vec3Module;
		vec3Module.Register(program);
		febcode::MathModule mathModule;
		mathModule.Register(program);

		febcode::Type vec3 = program.types.getVec3Type();
		if (!vec3) {
			err = "Error compiling code: 'vec3' type not defined";
			return false;
		}

		febcode::ParseSource(program, script);

		febcode::Compiler compiler(program);

		int globals[3] = { -1, -1, -1 };
		globals[0] = program.injectGlobal("_pos0" , vec3);
		globals[1] = program.injectGlobal("_time" , program.types.getDoubleType());
		globals[2] = program.injectGlobal("_norm0", vec3);

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
