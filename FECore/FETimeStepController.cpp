/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include "FETimeStepController.h"
#include "FELoadCurve.h"
#include "FEAnalysis.h"
#include "FEModel.h"
#include "FEPointFunction.h"
#include "DumpStream.h"
#include "log.h"

#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)>(b) ? (a) : (b))

REGISTER_SUPER_CLASS(FETimeStepController, FETIMECONTROLLER_ID);

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FETimeStepController, FEParamContainer)
	ADD_PARAMETER(m_maxretries, "max_retries");
	ADD_PARAMETER(m_iteopt    , "opt_iter");
	ADD_PARAMETER(m_dtmin     , "dtmin");
	ADD_PARAMETER(m_dtmax     , "dtmax");
	ADD_PARAMETER(m_naggr     , "aggressiveness");
	ADD_PARAMETER(m_dtforce   , "dtforce");
	ADD_PARAMETER(m_must_points, "must_points");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FETimeStepController::FETimeStepController(FEModel* fem) : FECoreBase(fem)
{
	m_step = nullptr; // must be set with SetAnalysis
	m_nretries = 0;
	m_maxretries = 5;
	m_naggr = 0;
	m_nmust = -1;
	m_next_must = -1;
	m_nmplc = -1;
	m_iteopt = 11;
	m_dtmax = m_dtmin = 0;

	m_ddt = 0;
	m_dtp = 0;

	m_dtforce = false;
}

//-----------------------------------------------------------------------------
void FETimeStepController::SetAnalysis(FEAnalysis* step)
{
	m_step = step;
}

//-----------------------------------------------------------------------------
//! copy from
void FETimeStepController::CopyFrom(FETimeStepController* tc)
{
	assert(m_step);

	m_naggr = tc->m_naggr;
	m_nmplc = tc->m_nmplc;
	m_iteopt = tc->m_iteopt;
	m_dtmin = tc->m_dtmin;
	m_dtmax = tc->m_dtmax;

	m_ddt = tc->m_ddt;
	m_dtp = tc->m_dtp;

	m_must_points = tc->m_must_points;
}

//-----------------------------------------------------------------------------
// initialization
bool FETimeStepController::Init()
{
	// make sure we have a step assigned
	if (m_step == 0) return false;

	// steal the load curve param from the dtmax parameter
	FEParam* p = FindParameterFromData((void*) &m_dtmax); assert(p);
	FEModel* fem = GetFEModel();
	FELoadController* plc = fem->GetLoadController(p);
	if (plc)
	{
		m_nmplc = plc->GetID();
		fem->DetachLoadController(p);

		// if a must-point curve is defined and the must-points are empty,
		// we copy the load curve points to the must-points
		if (m_must_points.empty())
		{
			FELoadCurve* lc = dynamic_cast<FELoadCurve*>(plc);
			if (lc)
			{
				FEPointFunction& f = lc->GetFunction();
				for (int i = 0; i < f.Points(); ++i)
				{
					double ti = f.LoadPoint(i).time;
					m_must_points.push_back(ti);
				}
			}
		}
	}

	// initialize "previous" time step
	m_dtp = m_step->m_dt0;

	return true;
}

//-----------------------------------------------------------------------------
//! reset
void FETimeStepController::Reset()
{
	m_dtp = m_step->m_dt0;
}

//-----------------------------------------------------------------------------
//! Restores data for a running restart

void FETimeStepController::Retry()
{
	FEModel* fem = m_step->GetFEModel();
	feLogEx(fem, "Retrying time step. Retry attempt %d of max %d\n\n", m_nretries + 1, m_maxretries);

	// adjust time step
	double dt = m_step->m_dt;
	if (m_nretries == 0) m_ddt = (dt) / (m_maxretries + 1);

	double dtn;
	if (m_naggr == 0) dtn = dt - m_ddt;
	else dtn = dt*0.5;

	feLogEx(fem, "\nAUTO STEPPER: retry step, dt = %lg\n\n", dtn);

	// increase retry counter
	m_nretries++;

	// the new time step cannot be a must-point
	if (m_nmust != -1)
	{
		// if we were at a must-point, make sure we can hit this must-point again
		m_next_must--;
		m_nmust = -1;
	}

	m_dtp = dtn;
	
	m_step->m_dt = dtn;
}

//-----------------------------------------------------------------------------
//! Adjusts the time step size based on the convergence information.
//!	If the previous time step was able to converge in less than
//! m_fem.m_iteopt iterations the step size is increased, else it
//! is decreased.

void FETimeStepController::AutoTimeStep(int niter)
{
	FEModel* fem = m_step->GetFEModel();
	double dt = m_step->m_dt;

	double dtn = m_dtp;
	double told = fem->GetCurrentTime();

	// make sure the timestep size is at least the minimum
	if (dtn < m_dtmin) dtn = m_dtmin;

	// get the max time step
	double dtmax = m_dtmax;

	// If we have a must-point load curve
	// we take the max step size from the lc
	if (m_nmplc >= 0)
	{
		FELoadCurve& mpc = *(dynamic_cast<FELoadCurve*>(fem->GetLoadController(m_nmplc)));
		FEPointFunction& lc = mpc.GetFunction();
		dtmax = lc.value(told);
	}

	// adjust time step size
	if (m_dtforce)
	{
		// if the force flag is set, we just set the time step to the max value
		dtn = dtmax;
	}
	else if (niter > 0)
	{
		double scale = sqrt((double)m_iteopt / (double)niter);

		// Adjust time step size
		if (scale >= 1)
		{
			dtn = dtn + (dtmax - dtn)*MIN(.20, scale - 1);
			dtn = MIN(dtn, 5.0*m_dtp);
			dtn = MIN(dtn, dtmax);
		}
		else
		{
			dtn = dtn - (dtn - m_dtmin)*(1 - scale);
			dtn = MAX(dtn, m_dtmin);
			dtn = MIN(dtn, dtmax);
		}

		// Report new time step size
		if (dtn > dt)
			feLogEx(fem, "\nAUTO STEPPER: increasing time step, dt = %lg\n\n", dtn);
		else if (dtn < dt)
			feLogEx(fem, "\nAUTO STEPPER: decreasing time step, dt = %lg\n\n", dtn);
	}

	// Store this time step value. This is the value that will be used to evaluate
	// the next time step increment. This will not include adjustments due to the must-point
	// controller since this could create really small time steps that may be difficult to
	// recover from. 
	m_dtp = dtn;

	// check for mustpoints
	if (m_must_points.empty() == false) dtn = CheckMustPoints(told, dtn);

	// make sure we are not exceeding the final time
	if (told + dtn > m_step->m_tend)
	{
		dtn = m_step->m_tend - told;
		feLogEx(fem, "MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtn);
	}

	// store time step size
	m_step->m_dt = dtn;
}

//-----------------------------------------------------------------------------
//! This function makes sure that no must points are passed. It returns an
//! updated value (less than dt) if t + dt would pass a must point. Otherwise
//! it returns dt.
//! \param t current time
//! \param dt current time step
//! \return updated time step.
double FETimeStepController::CheckMustPoints(double t, double dt)
{
	FEModel* fem = m_step->GetFEModel();

	double tnew = t + dt;
	double dtnew = dt;
	const double eps = m_step->m_tend*1e-12;
	double tmust = tnew + eps;
	m_nmust = -1;
	const int points = (int)m_must_points.size();
	if (m_next_must < points)
	{
		double lp;
		if (m_next_must < 0)
		{
			// find the first must-point that is on or past this time
			m_next_must = 0;
			bool bfound = false;
			do
			{
				lp = m_must_points[m_next_must];
				if ((tmust > lp) && (fabs(tnew - lp) > 1e-12)) ++m_next_must;
				else bfound = true;
			} while ((bfound == false) && (m_next_must < points));

			// make sure we did not pass all must points
			if (m_next_must >= points) return dt;
		}
		else lp = m_must_points[m_next_must];

		// TODO: what happens when dtnew < dtmin and the next time step fails??
		if (tmust > lp)
		{
			dtnew = lp - t;
			feLogEx(fem, "MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtnew);
			m_nmust = m_next_must++;
		}
		else if (fabs(tnew - lp) < 1e-12)
		{
			m_nmust = m_next_must++;
			tnew = lp;
			dtnew = tnew - t;
			feLogEx(fem, "MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtnew);
		}
		else if (tnew > m_step->m_tend)
		{
			dtnew = m_step->m_tend - t;
			feLogEx(fem, "MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtnew);
			m_nmust = m_next_must++;
		}
	}
	return dtnew;
}

//-----------------------------------------------------------------------------
//! serialize
void FETimeStepController::Serialize(DumpStream& ar)
{
	FECoreBase::Serialize(ar);
	ar & m_nretries;
	ar & m_nmplc;
	ar & m_nmust;
	ar & m_next_must;
	ar & m_ddt & m_dtp;
	ar & m_step;
	ar & m_must_points;
}
