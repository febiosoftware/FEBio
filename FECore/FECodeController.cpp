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
#include "FECodeController.h"
#include <febcode/parser.h>
#include <febcode/vm.h>
#include <FECore/log.h>

BEGIN_FECORE_CLASS(FECodeController, FELoadController)
	ADD_PARAMETER(script, "script");
END_FECORE_CLASS();

FECodeController::FECodeController(FEModel* fem) : FELoadController(fem)
{
}

bool FECodeController::Init()
{
	if (FELoadController::Init() == false) return false;
	return CompileScript();
}

bool FECodeController::CompileScript()
{
	if (script.empty()) {
		feLogError("Error compiling code: script is empty");
		return false;
	}

	try {
		febcode::ParseSource(program, script);

		febcode::Compiler compiler(program);

		globals[0] = program.injectGlobal("_time", program.types.Double());

		compiler.compile();

		// see if the globals were actually used
		if (program.globals[ globals[0] ].refcount == 0) globals[1] = -1;
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

double FECodeController::GetValue(double time)
{
	febcode::VM vm(program);
	if (globals[0] >= 0) vm.setGlobal(globals[0], time);
	febcode::Value v = vm.run();
	return febcode::getDouble(v);
}

void FECodeController::Serialize(DumpStream& ar)
{
	if (!ar.IsShallow())
	{
		ar& script;

		if (ar.IsLoading())
		{
			bool b = CompileScript(); assert(b);
		}
	}
}

void FECodeController::Reset()
{
}
