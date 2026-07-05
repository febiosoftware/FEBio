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
#include "FEScriptedLoadController.h"
#include <FECore/log.h>

FEScriptedLoadController::FEScriptedLoadController(FEModel* fem) : FEScripted<FELoadController>(fem)
{
	ScriptContext context;
	context.returnType = FEValueType::Double;
	context.allowMappedInputs = false;	// can't assign maps to input parameters
	context.allowVolatileInputs = false; // can't assign load controllers to input parameters
	context.addVariable("time", FEValueType::Double, true);
	SetScriptContext(context);
}

double FEScriptedLoadController::GetValue(double time)
{
	return Value(time);
}

void FEScriptedLoadController::Serialize(DumpStream& ar)
{
	// TODO: implement this
}

void FEScriptedLoadController::Reset()
{
	// TODO: implement this
}
