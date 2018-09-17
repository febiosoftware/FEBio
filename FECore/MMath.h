#pragma once
#include "MItem.h"

// Take the derivative with respect of x
MITEM MDerive(MITEM& e, MVariable& x);

// Take the nth-derivative with respect of x
MITEM MDerive(MITEM& e, MVariable& x, int n);

// Take the partial derivative with respect to x, y, z, ...
MITEM MDerive(MITEM& e, MSequence& x);

// Replace an expression with an expression
MITEM MReplace(MITEM& e, MVariable& x, MITEM& s);

// Replace an expression with an expression (assuming x is an equality)
MITEM MReplace(MITEM& e, MITEM& x);

// Replace an expression with an expression
MITEM MReplace(MITEM& e, MITEM& x, MITEM& s);

// replace multiple expressions with other expressions
MITEM MReplace(MITEM& e, MSequence& x, MSequence& s);

// Calculate the taylor expansion of an expression
MITEM MTaylor(MITEM& e, MVariable& x, double z, int n);

// Calculate the indefinite integral (without constant)
MITEM MIntegral(MITEM& e, MVariable& x);

// Calculate the definite integral
MITEM MIntegral(MITEM& e, MVariable& x, MITEM& a, MITEM& b);

// Expand an expression
MITEM MExpand(MITEM& e);

// Expand an expression but don't expand the second argument
MITEM MExpand(MITEM& e, MITEM& s);

// solve an expression
MITEM MSolve(MITEM& e, MITEM& v);

// collect terms
MITEM MCollect(MITEM& e, MITEM& x);

// simplify an expression
MITEM MSimplify(MITEM& e);
