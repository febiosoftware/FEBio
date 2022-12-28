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

class FEDataMap;

//---------------------------------------------------------------------------------------
// Base class for evaluating vec3d parameters
class FECORE_API FEMat3dsValuator : public FEValuator
{
	FECORE_SUPER_CLASS(FEMAT3DSVALUATOR_ID)
	FECORE_BASE_CLASS(FEMat3dsValuator)

public:
	FEMat3dsValuator(FEModel* fem) : FEValuator(fem) {};

	virtual mat3ds operator()(const FEMaterialPoint& pt) = 0;

	virtual FEMat3dsValuator* copy() = 0;

	virtual bool isConst() { return false; }

	virtual mat3ds* constValue() { return nullptr; }
};

//-----------------------------------------------------------------------------
// A constant valuator
class FECORE_API FEConstValueMat3ds : public FEMat3dsValuator
{
public:
	FEConstValueMat3ds(FEModel* fem);

	FEMat3dsValuator* copy() override;

	mat3ds operator()(const FEMaterialPoint& pt) override { return m_val; }

	// is this a const value
	bool isConst() override { return true; }

	// get the const value (returns 0 if param is not const)
	mat3ds* constValue() override { return &m_val; }

	mat3ds& value() { return m_val; }

private:
	mat3ds	m_val;

	DECLARE_FECORE_CLASS();
};

//---------------------------------------------------------------------------------------
class FECORE_API FEMappedValueMat3ds : public FEMat3dsValuator
{
public:
	FEMappedValueMat3ds(FEModel* fem);

	void setDataMap(FEDataMap* val);

	mat3ds operator()(const FEMaterialPoint& pt) override;

	FEMat3dsValuator* copy() override;

	void Serialize(DumpStream& ar) override;

private:
	std::string	m_mapName;

private:
	FEDataMap*	m_val;

	DECLARE_FECORE_CLASS();
};
