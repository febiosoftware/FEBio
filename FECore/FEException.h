/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "fecore_api.h"
#include <vector>
#include <string>
using namespace std;

class FEElement;

class FECORE_API FEException
{
public:
	FEException(const char* msg = nullptr);
	virtual ~FEException();

	const char* what();

	void what(const char* msg, ...);

private:
	std::string	m_what;
};

class FECORE_API NegativeJacobian : public FEException
{
public:
	NegativeJacobian(int iel, int ng, double vol, FEElement* pe = 0);

	int		m_iel;	// element where the jacobian was negative
	int		m_ng;	// integration point
	double	m_vol;	// volume
	FEElement*	m_pel;	// pointer to element

	static bool DoOutput();

public:
	static bool m_boutput;	//!< set to false to suppress output of negative jacobians
};

class FECORE_API ZeroDiagonal : public FEException
{
private:
	struct EQUATION
	{
		int	node;	// node
		int	dof;	// degree of node
	};

public:
	ZeroDiagonal(int node, int dof);

	char m_szerr[256];	// the error message
};

class FECORE_API EnergyDiverging : public FEException {};

class FECORE_API MaxStiffnessReformations : public FEException {};

class FECORE_API ZeroLinestepSize : public FEException {};

class FECORE_API ForceConversion {};

class FECORE_API IterationFailure {};

class FECORE_API MaxResidualError {};

class FECORE_API NANDetected {};

class FECORE_API FatalError {};

class FECORE_API FEMultiScaleException
{
public:
	FEMultiScaleException(int eid, int gpt) : elemId(eid), gptIndex(gpt) {}

public:
	int elemId;
	int gptIndex;
};

class FECORE_API LinearSolverFailed {};

class FECORE_API DoRunningRestart{};
