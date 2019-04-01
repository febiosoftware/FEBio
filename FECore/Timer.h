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
#include "FECoreKernel.h"
#include <vector>
#include <string>

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

public:
	static void time_str(double fsec, char* sz);
	static void GetTime(double fsec, int& nhour, int& nmin, int& nsec);

private:
	void*	m_pimpl;	//!< local timing data (using PIMPL ididom to hide OS specifics)

	bool	m_brunning;	//!< flag indicating whether start was called
	double	m_sec;		//!< accumulated time so far in seconds
};

//-----------------------------------------------------------------------------
// This is helper class that can be used to ensure that a timer is stopped when
// the function that is being timed exits. That way, the Timer::stop member does not 
// have to be called at every exit point of a function.
class FECORE_API TimerTracker
{
public:
	TimerTracker(Timer* timer) : m_timer(timer) { timer->start(); };
	~TimerTracker() { m_timer->stop(); }

private:
	Timer*	m_timer;
};

#define TRACK_TIME(timerId) TimerTracker _trackTimer(GetFEModel()->GetTimer(timerId));
