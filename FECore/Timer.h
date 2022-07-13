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
#include "FECoreKernel.h"
#include <vector>
#include <string>
#include <chrono>

using dseconds = std::chrono::duration<double>;
using dminutes = std::chrono::duration<double, std::ratio<60>>;
using dhours = std::chrono::duration<double, std::ratio<3600>>;

//-----------------------------------------------------------------------------
//! This class implements a simple timer. 

//! The start function starts the timer, the stop
//! function stops it and the GetTime function returns the time elapsed between
//! the call to start and stop

class FECORE_API Timer
{
public:
	//! constructor
	Timer();

	//! constructor
	~Timer();

	//! Start the timer
	void start();

	//! Stop the timer
	void stop();

	//! Reset the timer
	void reset();

	//! Get the elapsed time
	void GetTime(int& nhour, int& nmin, int& nsec);

	//! Get the time in seconds
	double GetTime();

	//! return the time as a text string
	void time_str(char* sz);

	//! check the elapsed time
	double peek();

	//! see if the timer is running
	bool isRunning() const { return m_brunning; }

public:
	static void time_str(double fsec, char* sz);
	static void GetTime(double fsec, int& nhour, int& nmin, int& nsec);

private:
	void*	m_pimpl;	//!< local timing data (using PIMPL ididom to hide OS specifics)

	bool	m_brunning;	//!< flag indicating whether start was called
    dseconds m_total; //!< accumulated time so far in seconds
};

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
// Timer IDs
enum TimerID {
	Timer_Update,
	Timer_LinSolve,
	Timer_Reform,
	Timer_Residual,
	Timer_Stiffness,
	Timer_QNUpdate,
	Timer_ModelSolve
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
	TimerTracker(Timer* timer);
	~TimerTracker();

private:
	Timer*	m_timer;
};

#define TRACK_TIME(timerId) TimerTracker _trackTimer(GetFEModel(), timerId);
