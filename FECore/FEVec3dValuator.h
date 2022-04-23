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
#include "FEValuator.h"
#include "FEDataMap.h"
#include "MathObject.h"

//---------------------------------------------------------------------------------------
// Base class for evaluating vec3d parameters
class FECORE_API FEVec3dValuator : public FEValuator
{
	FECORE_SUPER_CLASS(FEVEC3DVALUATOR_ID)
	FECORE_BASE_CLASS(FEVec3dValuator)

public:
	FEVec3dValuator(FEModel* fem);

public:
	// evaluate value at a material point
	virtual vec3d operator()(const FEMaterialPoint& pt) = 0;

	// create a copy of the valuator
	virtual FEVec3dValuator* copy() = 0;

	// is the valuator constant
	virtual bool isConst() { return false; }

	// return the const value
	virtual vec3d* constValue() { return nullptr; }

	// return a unit vector
	vec3d unitVector(const FEMaterialPoint& pt)
	{
		vec3d v = operator () (pt);
		return v.Normalized();
	}
};
