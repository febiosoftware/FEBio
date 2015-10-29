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
