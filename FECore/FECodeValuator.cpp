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
#include <FECore/FEMaterialPoint.h>
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <FECore/FESurfaceLoad.h>
#include <febcode/types.h>

class FECodeValuator::Imp {
public:
	std::string scriptCode;

	febcode::Program program;
	int globals[3]; // for _pos0, _time, _norm0

	FESurface* surf = nullptr; // for surface valuators
};

BEGIN_FECORE_CLASS(FECodeValuator, FEScalarValuator)
	ADD_PARAMETER(m_scriptName, "code");
END_FECORE_CLASS()

FECodeValuator::FECodeValuator(FEModel* fem) : FEScalarValuator(fem), m(*new Imp())
{
}

FECodeValuator::~FECodeValuator()
{

}

bool FECodeValuator::CompileScript()
{
	if (m.scriptCode.empty()) return false;

	FECoreBase* parent = GetParent();
	try {

		febcode::Vec2Module vec2Module;
		vec2Module.Register(m.program);
		febcode::Vec3Module vec3Module;
		vec3Module.Register(m.program);
		febcode::MathModule mathModule;
		mathModule.Register(m.program);

		febcode::Type vec3_t = m.program.types.getVec3Type();
		if (!vec3_t) {
			feLogError("Error compiling code: 'vec3' type not defined");
			return false;
		}

		febcode::ParseSource(m.program, m.scriptCode);

		febcode::Compiler compiler(m.program);


		m.globals[0] = m.program.addGlobal("_pos0", vec3_t, febcode::vec3(0., 0., 0.), true);
		m.globals[1] = m.program.addGlobal("_time", 0.0);
		m.globals[2] = -1;

		if (auto s = dynamic_cast<FESurfaceLoad*>(parent))
		{
			m.globals[2] = m.program.addGlobal("_norm0", vec3_t, febcode::vec3(0., 0., 0.), true);
			m.surf = &s->GetSurface();
		}

		compiler.compile();

		// see if the globals were actually used
		if (m.program.globals["_pos0"].refcount == 0) m.globals[0] = -1;
		if (m.program.globals["_time"].refcount == 0) m.globals[1] = -1;
		if ((m.globals[2] >= 0) && (m.program.globals["_norm0"].refcount == 0)) m.globals[2] = -1;
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


bool FECodeValuator::Init()
{
	if (FEScalarValuator::Init() == false) return false;

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
	febcode::VM vm(m.program);
	if (m.globals[0] >= 0) vm.setGlobal(m.globals[0], febcode::vec3(pt.m_r0.x, pt.m_r0.y, pt.m_r0.z));
	if (m.globals[1] >= 0) vm.setGlobal(m.globals[1], GetFEModel()->GetTime().currentTime);
	if (m.surf && (m.globals[2] >= 0))
	{
		FESurfaceElement* el = dynamic_cast<FESurfaceElement*>(pt.m_elem);
		if (el)
		{
			vec3d n = m.surf->SurfaceNormal(*el, pt.m_index);
			vm.setGlobal(m.globals[2], febcode::vec3(n.x, n.y, n.z));
		}
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
		globals[0] = program.addGlobal("_pos0", vec3, febcode::vec3( 0., 0., 0. ));
		globals[1] = program.addGlobal("_time", 0.0);
		globals[2] = program.addGlobal("_norm0", vec3, febcode::vec3(0., 0., 0.));

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
