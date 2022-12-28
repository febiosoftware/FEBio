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
#include "FECoreBase.h"
class FEAnalysis;
class DumpStream;

//-------------------------------------------------------------------
// Class to control the time step
class FECORE_API FETimeStepController : public FECoreBase
{
	FECORE_SUPER_CLASS(FETIMECONTROLLER_ID)
	FECORE_BASE_CLASS(FETimeStepController)

public:
	FETimeStepController(FEModel* fem);

	void SetAnalysis(FEAnalysis* step);

	// initialization
	bool Init() override;

	//! reset
	void Reset();

	//! serialize
	void Serialize(DumpStream& ar) override;

	//! copy from
	void CopyFrom(FETimeStepController* tc);

public:
	//! Do a running restart
	void Retry();

	//! Update Time step
	void AutoTimeStep(int niter);

	//! Adjust for must points
	double CheckMustPoints(double t, double dt);

private:
	FEAnalysis*	m_step;

public:
	int		m_nretries;		//!< nr of retries tried so far
	int		m_maxretries;	//!< max nr of retries allowed per time step
	int		m_naggr;		//!< aggressivness parameter
	int		m_nmplc;		//!< must point load curve number
	int		m_nmust;		//!< current must-point
	int		m_next_must;	//!< next must-point to visit
	int		m_iteopt;		//!< optimum nr of iterations
	double	m_dtmin;		//!< min time step size
	double	m_dtmax;		//!< max time step size
	double	m_cutback;		//!< cut back factor used in aggressive time stepping

	std::vector<double>	m_must_points;	//!< the list of must-points
	bool				m_mp_repeat;	//!< repeat must-points
	double				m_mp_toff;		//!< offset for repeat must-points

private:
	double	m_ddt;			//!< used by auto-time stepper
	double	m_dtp;			//!< previous time step size

	bool	m_dtforce;		//!< force max time step

	DECLARE_FECORE_CLASS();
};
