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
#include "FELoadController.h"
#include "PointCurve.h"

//-----------------------------------------------------------------------------
// Base class for load curves.
// Load curves are used to manipulate the time dependency of model parameters.
class FECORE_API FELoadCurve : public FELoadController
{
public:
	// constructor
	FELoadCurve(FEModel* fem);
	FELoadCurve(const FELoadCurve& lc);

	void operator = (const FELoadCurve& lc);

	// destructor
	virtual ~FELoadCurve();

	void Serialize(DumpStream& ar) override;

	bool CopyFrom(FELoadCurve* lc);

	void Add(double time, double value);

	void Clear();
    
    bool Init() override;

	PointCurve& GetFunction() { return m_fnc; }

	int GetInterpolation() const { return m_int; }
	void SetInterpolation(PointCurve::INTFUNC f);

	int GetExtendMode() const { return m_ext; }
	void SetExtendMode(PointCurve::EXTMODE f);

	std::vector<vec2d> GetPoints() const { return m_points; }

	double GetValue(double time) override;

private:
	int		m_int;
	int		m_ext;
	std::vector<vec2d>	m_points;

private:
	PointCurve	m_fnc;	//!< function to evaluate

	DECLARE_FECORE_CLASS();
};
