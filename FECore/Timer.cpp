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
#include "Timer.h"
#include <stdio.h>
#include <string>
#include <chrono>
#include "FEModel.h"
using namespace std::chrono;

using dseconds = duration<double>;

// data storing timing info
struct	Timer::Imp {
	time_point<steady_clock>	start;	//!< time at start
	time_point<steady_clock>	stop;	//!< time at last stop
	time_point<steady_clock>	pause;	//!< time when paused

	bool	isRunning = false;	//!< flag indicating whether start was called
	dseconds total; //!< total accumulated time (between start and stop)

	bool isPaused = false;
	dseconds paused; //!< accumulated time while paused

	Timer* parent = nullptr; // the timer that was active when this timer starts

	static Timer* activeTimer;
};

Timer* Timer::Imp::activeTimer = nullptr;

Timer::Timer()
{
	m = new Imp;
	reset();
}

Timer::~Timer()
{
	delete m;
}

Timer* Timer::activeTimer()
{
	return Imp::activeTimer;
}

void Timer::start()
{
	m->parent = m->activeTimer;
	if (m->parent) m->parent->pause();
	m->activeTimer = this;
	m->start = steady_clock::now();
	assert(m->isRunning == false);
	m->isRunning = true;
}

void Timer::stop()
{
	m->stop = steady_clock::now();
	assert(m->isRunning == true);
	m->isRunning = false;
	m->total += m->stop - m->start;

	assert(m->activeTimer == this);
	m->activeTimer = m->parent;
	if (m->parent) m->parent->unpause();
}

void Timer::pause()
{
	assert(!m->isPaused);
	m->pause = steady_clock::now();
	m->isPaused = true;
}

void Timer::unpause()
{
	assert(m->isPaused);
	auto tmp = steady_clock::now();
	m->paused += tmp - m->pause;
	m->isPaused = false;
}

void Timer::reset()
{
	m->total = dseconds(0);
	m->paused = dseconds(0);
	m->isRunning = false;
	m->isPaused = false;
	m->parent = nullptr;
}

bool Timer::isRunning() const { return m->isRunning; }

double Timer::peek()
{
	if (m->isRunning)
	{
		time_point<steady_clock> tmp = steady_clock::now();
		return duration_cast<dseconds>(m->total + (tmp - m->start)).count();
	}
	else 
	{
		return m->total.count();
	}
}

void Timer::GetTime(int& nhour, int& nmin, int& nsec)
{
	double sec = peek();
	GetTime(sec, nhour, nmin, nsec);
}

double Timer::GetExclusiveTime()
{
	assert(!m->isPaused);
	double sec = peek() - m->paused.count();
	return sec;
}

void Timer::GetTime(double fsec, int& nhour, int& nmin, int& nsec)
{
	nhour = (int) (fsec / 3600.0); fsec -= nhour*3600;
	nmin  = (int) (fsec /   60.0); fsec -= nmin*60;
	nsec  = (int) (fsec + 0.5);
}

double Timer::GetTime()
{
	return (m->isRunning? peek() : m->total.count());
}

void Timer::time_str(char* sz)
{
	int nhour, nmin, nsec;
	GetTime(nhour, nmin, nsec);
	snprintf(sz, 64, "%d:%02d:%02d", nhour, nmin, nsec);
}

void Timer::time_str(double fsec, char* sz)
{
	int nhour, nmin, nsec;
	GetTime(fsec, nhour, nmin, nsec);
	snprintf(sz, 64, "%d:%02d:%02d", nhour, nmin, nsec);
}

//============================================================================
TimerTracker::TimerTracker(FEModel* fem, int timerId) : TimerTracker(fem->GetTimer(timerId)) {}
