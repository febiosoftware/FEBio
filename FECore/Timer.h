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
#include "fecore_api.h"

class FEModel;

//! This class implements a simple timer. 

//! The start function starts the timer, the stop
//! function stops it and the GetTime function returns the time elapsed between
//! the call to start and stop

class FECORE_API Timer
{
	struct Imp;

public:
	//! constructor
	Timer();

	//! constructor
	~Timer();

	//! Start the timer
	void start();

	//! Stop the timer
	void stop();

	//! pause the timer
	void pause();

	//! continue
	void unpause();

	void reset();

	//! Get the elapsed time
	void GetTime(int& nhour, int& nmin, int& nsec);

	//! Get the time in seconds
	double GetTime();

	//! Get the exclusive time (i.e. time when not paused)
	double GetExclusiveTime();

	//! return the time as a text string
	void time_str(char* sz);

	//! check the elapsed time
	double peek();

	//! see if the timer is running
	bool isRunning() const;

public:
	static void time_str(double fsec, char* sz);
	static void GetTime(double fsec, int& nhour, int& nmin, int& nsec);
	static Timer* activeTimer();

private:
	Imp*	m;	//!< local timing data (using PIMPL ididom)
};

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
// Timer IDs
enum TimerID {
	Timer_Init,
	Timer_Update,
	Timer_LinSol_Factor,
	Timer_LinSol_Backsolve,
	Timer_Reform,
	Timer_Residual,
	Timer_Stiffness,
	Timer_QNUpdate,
	Timer_Serialize,
	Timer_ModelSolve,
	Timer_Callback,
	Timer_USER1,
	Timer_USER2,
	Timer_USER3,
	Timer_USER4,
	TIMER_COUNT // leave this at the end so that it equals the nr. of timers we need
};

//-----------------------------------------------------------------------------
// This is helper class that can be used to ensure that a timer is stopped when
// the function that is being timed exits. That way, the Timer::stop member does not 
// have to be called at every exit point of a function.
// In addition, it will also check if the timer is already running (e.g. from a function
// higher in the call stack) in which case it will track the timer. 
class FECORE_API TimerTracker
{
public:
	TimerTracker(FEModel* fem, int timerId);
	TimerTracker(Timer* timer)
	{
		if (timer && !timer->isRunning()) { m_timer = timer; timer->start(); }
		else m_timer = nullptr;
	}
	~TimerTracker() { if (m_timer) m_timer->stop(); }

private:
	Timer*	m_timer;
};

#define TRACK_TIME(timerId) TimerTracker _trackTimer(GetFEModel(), timerId);
