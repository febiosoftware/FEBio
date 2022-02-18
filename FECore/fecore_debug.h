/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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



#pragma once

// Load debugger template constructions.
// Don't use anything in this file directly.
// Instead use the functions and macros defined below
#include "fecore_debug_t.h"

//-----------------------------------------------------------------------------
// Set a break point. This will bring up the debugger prompt.
// (type help for a list of available debugger commands)
#define fecore_break() { \
	static FECoreBreakPoint _br_##__LINE__; \
	if (_br_##__LINE__.IsActive()) { \
		cout << "breakpoint " << _br_##__LINE__.GetID() << " : " << __FILE__ <<  "(line " << __LINE__ << ")\n"; \
		_br_##__LINE__.Break(); }}

//-----------------------------------------------------------------------------
// Add a variable to the watch list.
// A variable on the watch list can be inspected from the debugger 
// (type print var on the debugger prompt, where var is the variable name)
#define fecore_watch(a) FECoreWatchVariable _##a(create_watch_variable(&a, #a));

//-----------------------------------------------------------------------------
// print the contents of a variable to the screen
#define fecore_print(a) { \
	cout << #a << endl << typeid(a).name() << endl; \
	fecore_print_T(&a); \
	cout << endl; }
