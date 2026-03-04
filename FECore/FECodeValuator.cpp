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

BEGIN_FECORE_CLASS(FECodeValuator, FEScalarValuator)
	ADD_PARAMETER(m_script, "code");
END_FECORE_CLASS()

FECodeValuator::FECodeValuator(FEModel* fem) : FEScalarValuator(fem)
{
}

FECodeValuator::~FECodeValuator()
{

}

bool FECodeValuator::CompileScript()
{
	FECoreBase* parent = GetParent();
	try {

		febcode::Vec2Module vec2Module;
		vec2Module.Register(m_program);
		febcode::Vec3Module vec3Module;
		vec3Module.Register(m_program);
		febcode::MathModule mathModule;
		mathModule.Register(m_program);

		febcode::Type vec3 = m_program.types.getStructType("vec3");
		if (!vec3) {
			feLogError("Error compiling code: 'vec3' type not defined");
			return false;
		}

		febcode::ParseSource(m_program, m_script);

		febcode::Compiler compiler(m_program);


		globals[0] = m_program.addGlobal("_pos0", vec3, { 0., 0., 0. });
		globals[1] = m_program.addGlobal("_time", 0.0);
		globals[2] = -1;

		if (auto s = dynamic_cast<FESurfaceLoad*>(parent))
		{
			globals[2] = m_program.addGlobal("_norm0", vec3, { 0., 0., 0. });
			m_surf = &s->GetSurface();
		}

		compiler.compile();

		// see if the globals were actually used
		if (m_program.globals["_pos0"].refcount == 0) globals[0] = -1;
		if (m_program.globals["_time"].refcount == 0) globals[1] = -1;
		if ((globals[2] >= 0) && (m_program.globals["_norm0"].refcount == 0)) globals[2] = -1;
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

	if (CompileScript() == false) return false;

	return true;
}

double FECodeValuator::operator()(const FEMaterialPoint& pt)
{
	febcode::VM vm(m_program);
	if (globals[0] >= 0) vm.setGlobal(globals[0], { pt.m_r0.x, pt.m_r0.y, pt.m_r0.z });
	if (globals[1] >= 0) vm.setGlobal(globals[1], GetFEModel()->GetTime().currentTime);
	if (m_surf && (globals[2] >= 0))
	{
		FESurfaceElement* el = dynamic_cast<FESurfaceElement*>(pt.m_elem);
		if (el)
		{
			vec3d n = m_surf->SurfaceNormal(*el, pt.m_index);
			vm.setGlobal(globals[2], { n.x, n.y, n.z });
		}
	}
	febcode::Value v = vm.run();
	return febcode::getDouble(v);
}

FEScalarValuator* FECodeValuator::copy()
{
	FECodeValuator* p = new FECodeValuator(GetFEModel());
	p->m_script = m_script;
	return p;
}

void FECodeValuator::Serialize(DumpStream& ar)
{
	if (!ar.IsShallow() && !ar.IsSaving())
	{
		FEScalarValuator::Serialize(ar);
		CompileScript();
	}
}
