// Timer.h: interface for the Timer class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TIMER_H__5C4CDF72_0B19_4C5B_9B21_DB7B85FCEC4D__INCLUDED_)
#define AFX_TIMER_H__5C4CDF72_0B19_4C5B_9B21_DB7B85FCEC4D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

//-----------------------------------------------------------------------------
//! This class implements a simple timer. 

//! The start function starts the timer, the stop
//! function stops it and the GetTime function returns the time elapsed between
//! the call to start and stop

class Timer  
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
class TimerTracker
{
public:
 	TimerTracker(Timer& timer) : m_timer(timer) { timer.start(); }
	~TimerTracker() { m_timer.stop(); }
private:
	Timer&	m_timer;
};

#endif // !defined(AFX_TIMER_H__5C4CDF72_0B19_4C5B_9B21_DB7B85FCEC4D__INCLUDED_)
