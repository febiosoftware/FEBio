/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in
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
#include <FECore/FESolver.h>
#include "FEFluidSolver.h"
#include "FESolutesSolver.h"
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! The FEFluidSolutesSolver class solves fluid mechanics problems
//! with mass transport (solutes)
//! It can deal with quasi-static and dynamic problems
//!
class FEBIOFLUID_API FEFluidSolutesSolver2 : public FESolver
{
public:
	//! constructor
	FEFluidSolutesSolver2(FEModel* pfem);

	//! destructor
	~FEFluidSolutesSolver2();

	//! Initializes data structures
	bool Init() override;

	//! initialize the step
	bool InitStep(double time) override;

	//! Initialize linear equation system
	bool InitEquations() override;

	bool SolveStep() override;

	//! Serialization
	void Serialize(DumpStream& ar) override;

protected:
	void MapVelocitySolution();

private:
	FEFluidSolver	m_fldSolver;
	FESolutesSolver	m_sltSolver;

// these parameters must be transferred to fluid and solute solvers
private:
	// shared parameters
	int		m_maxRef;
	int		m_maxUps;
	bool	m_divergeReform;
	bool	m_reformTimeStep;
	double	m_etol;
	double	m_rtol;
	double	m_lstol;
	double	m_minRes;
	double	m_rhoi;
	int		m_qnmethod;

	// fluid solver parameters
	double	m_vtol;
	double	m_ftol;

	// solute solver parameters
	double	m_ctol;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
