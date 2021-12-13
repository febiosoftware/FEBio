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
#include "MathObject.h"
#include <string>

// this class converts a MathObject to a string
class FECORE_API MObj2String
{
public:
	std::string Convert(const MathObject& o);

protected:
	std::string Convert(const MItem* pi);

	std::string Constant (const MConstant* pc);
	std::string Fraction (const MFraction* pc);
	std::string NamedCt  (const MNamedCt*  pc);
	std::string Variable (const MVarRef*   pv);
	std::string OpNeg    (const MNeg*      po);
	std::string OpAdd    (const MAdd*      po);
	std::string OpSub    (const MSub*      po);
	std::string OpMul    (const MMul*      po);
	std::string OpDiv    (const MDiv*      po);
	std::string OpPow    (const MPow*      po);
	std::string OpEqual  (const MEquation* po);
	std::string OpFnc1D  (const MFunc1D*   po);
	std::string OpFnc2D  (const MFunc2D*   po);
	std::string OpFncND  (const MFuncND*   po);
	std::string OpSFnc   (const MSFuncND*  po);
};
