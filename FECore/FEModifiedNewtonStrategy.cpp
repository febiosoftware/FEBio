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
#include "stdafx.h"
#include "FEModifiedNewtonStrategy.h"
#include "FESolver.h"
#include "FEException.h"
#include "FENewtonSolver.h"

//-----------------------------------------------------------------------------
FEModifiedNewtonStrategy::FEModifiedNewtonStrategy(FEModel* fem) : FENewtonStrategy(fem)
{
	m_plinsolve = nullptr;
	m_maxups = 0;
}

//-----------------------------------------------------------------------------
bool FEModifiedNewtonStrategy::Init()
{
	if (m_pns == nullptr) return false;
	m_plinsolve = m_pns->GetLinearSolver();
	return true;
}

//-----------------------------------------------------------------------------
bool FEModifiedNewtonStrategy::Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1)
{
	// nothing to do here.
	return true;
}

//-----------------------------------------------------------------------------
void FEModifiedNewtonStrategy::SolveEquations(vector<double>& x, vector<double>& b)
{
	// perform a backsubstitution
	if (m_plinsolve->BackSolve(x, b) == false)
	{
		throw LinearSolverFailed();
	}
}
