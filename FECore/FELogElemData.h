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
#include "FELogData.h"
#include "FEFunction1D.h"

class FEElement;

class FECORE_API FELogElemSource : public FELogData
{
public:
	FELogElemSource(FEModel* fem) : FELogData(fem) {}
	virtual double value(FEElement& el) = 0;
};

//! Base class for element log data
class FECORE_API FELogElemData : public FELogElemSource
{
	FECORE_SUPER_CLASS(FELOGELEMDATA_ID)
	FECORE_BASE_CLASS(FELogElemData)

public:
	FELogElemData(FEModel* fem);
	virtual ~FELogElemData();
};

//! Base class for element log data definitions 
class FECORE_API FELogElemDefinition : public FELogElemSource
{
	FECORE_SUPER_CLASS(FELOGELEMDEFINITION_ID)
	FECORE_BASE_CLASS(FELogElemData)

public:
	FELogElemDefinition(FEModel* fem) : FELogElemSource(fem) {}
};

class FECORE_API FELogElemAlias : public FELogElemDefinition
{
public:
	FELogElemAlias(FEModel* fem);

	bool Init() override;

	double value(FEElement& el) override;

private:
	std::string m_var;

	FELogElemData* m_pdata = nullptr;

	DECLARE_FECORE_CLASS();
};

class FECORE_API FELogElemFunction : public FELogElemDefinition
{
public:
	FELogElemFunction(FEModel* fem);
	bool Init() override;
	double value(FEElement& el) override;

private:
	std::string m_var;
	FEFunction1D* m_func = nullptr;

	FELogElemSource* m_pdata = nullptr;

	DECLARE_FECORE_CLASS();
};
