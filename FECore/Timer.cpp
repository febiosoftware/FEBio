// Timer.cpp: implementation of the Timer class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Timer.h"
#include <stdio.h>
#include <time.h>

//-----------------------------------------------------------------------------
// Define the data types used for measuring times.
// This depends on the system
#ifdef WIN32
#define TIMER_TYPE clock_t
#else
#define TIMER_TYPE	time_t
#endif

//-----------------------------------------------------------------------------
// forward declaration of the functions to retrieve timing info
void sys_get_time(TIMER_TYPE& t);
double sys_diff_time(TIMER_TYPE& t1, TIMER_TYPE& t0);

//-----------------------------------------------------------------------------
// OS-dependent timing functions
#ifdef WIN32
void sys_get_time(TIMER_TYPE& t) { t = clock(); }
double sys_diff_time(TIMER_TYPE& t1, TIMER_TYPE& t0) { return (double) (t1 - t0) / CLOCKS_PER_SEC; }
#else
void sys_get_time(TIMER_TYPE& t) { time(&t); }
double sys_diff_time(TIMER_TYPE& t1, TIMER_TYPE& t0) { return difftime(t1, t0); }
#endif

//-----------------------------------------------------------------------------
// data storing timing info
struct	timer_data {
	TIMER_TYPE	m_start;	//!< time at start
	TIMER_TYPE	m_stop;		//!< time at last stop
};

//-----------------------------------------------------------------------------
Timer::Timer()
{
	m_pimpl = new timer_data;
	reset(); 
}

//-----------------------------------------------------------------------------
Timer::~Timer()
{
	delete (timer_data*)m_pimpl;
}

//-----------------------------------------------------------------------------
void Timer::start()
{
	timer_data& t = *(static_cast<timer_data*>(m_pimpl));
	sys_get_time(t.m_start);
	m_brunning = true;
}

//-----------------------------------------------------------------------------
void Timer::stop()
{
	timer_data& t = *(static_cast<timer_data*>(m_pimpl));
	sys_get_time(t.m_stop);
	m_brunning = false;

	m_sec += sys_diff_time(t.m_stop, t.m_start);
}

//-----------------------------------------------------------------------------
void Timer::reset()
{
	m_sec = 0;
	m_brunning = false;
}

//-----------------------------------------------------------------------------
double Timer::peek()
{
	TIMER_TYPE pause;
	sys_get_time(pause);
	timer_data& t = *(static_cast<timer_data*>(m_pimpl));
	
	return m_sec + sys_diff_time(pause, t.m_start);
}

//-----------------------------------------------------------------------------
void Timer::GetTime(int& nhour, int& nmin, int& nsec)
{
	double sec = (m_brunning? peek() : m_sec);
	GetTime(sec, nhour, nmin, nsec);
}

//-----------------------------------------------------------------------------
void Timer::GetTime(double fsec, int& nhour, int& nmin, int& nsec)
{
	nhour = (int) (fsec / 3600.0); fsec -= nhour*3600;
	nmin  = (int) (fsec /   60.0); fsec -= nmin*60;
	nsec  = (int) (fsec + 0.5);
}

//-----------------------------------------------------------------------------
double Timer::GetTime()
{
	return (m_brunning? peek() : m_sec);
}

//-----------------------------------------------------------------------------
void Timer::time_str(char* sz)
{
	int nhour, nmin, nsec;
	GetTime(nhour, nmin, nsec);
	sprintf(sz, "%d:%02d:%02d", nhour, nmin, nsec);
}

//-----------------------------------------------------------------------------
void Timer::time_str(double fsec, char* sz)
{
	int nhour, nmin, nsec;
	GetTime(fsec, nhour, nmin, nsec);
	sprintf(sz, "%d:%02d:%02d", nhour, nmin, nsec);
}

//-----------------------------------------------------------------------------
const std::string& Timer::name() const
{
	return m_name;
}

//-----------------------------------------------------------------------------
void Timer::setName(const std::string& name)
{
	m_name = name;
}


//=================================================================================================

std::vector<Timer*> TimerManager::m_timers;


Timer* TimerManager::findTimer(const std::string& name)
{
	// see if the timer already exists
	for (size_t i=0; i<m_timers.size(); ++i)
	{
		if (m_timers[i]->name() == name) return m_timers[i];
	}

	// create new timer
	Timer* newTimer = new Timer;
	newTimer->setName(name);

	// add it to the list
	m_timers.push_back(newTimer);

	// return it
	return newTimer;
}

int TimerManager::Timers()
{
	return (int) m_timers.size();
}

Timer* TimerManager::getTimer(int i)
{
	return m_timers[i];
}
