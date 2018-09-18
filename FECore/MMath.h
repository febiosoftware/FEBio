#pragma once
#include "MItem.h"

// Take the derivative with respect of x
MITEM MDerive(const MITEM& e, const MVariable& x);

// Take the nth-derivative with respect of x
MITEM MDerive(const MITEM& e, const MVariable& x, int n);

// Take the partial derivative with respect to x, y, z, ...
MITEM MDerive(const MITEM& e, const MSequence& x);

// Replace an expression with an expression
MITEM MReplace(const MITEM& e, const MVariable& x, const MITEM& s);

// Replace an expression with an expression (assuming x is an equality)
MITEM MReplace(const MITEM& e, const MITEM& x);

// Replace an expression with an expression
MITEM MReplace(const MITEM& e, const MITEM& x, const MITEM& s);

// replace multiple expressions with other expressions
MITEM MReplace(const MITEM& e, const MSequence& x, const MSequence& s);

// Calculate the taylor expansion of an expression
MITEM MTaylor(const MITEM& e, const MVariable& x, double z, int n);

// Calculate the indefinite integral (without constant)
MITEM MIntegral(const MITEM& e, const MVariable& x);

// Calculate the definite integral
MITEM MIntegral(const MITEM& e, const MVariable& x, const MITEM& a, const MITEM& b);

// Expand an expression
MITEM MExpand(const MITEM& e);

// Expand an expression but don't expand the second argument
MITEM MExpand(const MITEM& e, const MITEM& s);

// solve an expression
MITEM MSolve(const MITEM& e, const MITEM& v);

// collect terms
MITEM MCollect(const MITEM& e, const MITEM& x);

// simplify an expression
MITEM MSimplify(const MITEM& e);
