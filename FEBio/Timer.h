// Timer.h: interface for the Timer class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TIMER_H__5C4CDF72_0B19_4C5B_9B21_DB7B85FCEC4D__INCLUDED_)
#define AFX_TIMER_H__5C4CDF72_0B19_4C5B_9B21_DB7B85FCEC4D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <time.h>

//-----------------------------------------------------------------------------
//! This implements a simple timer. 

//! The start function starts the timer, the stop
//! function stops it and the GetTime function returns the time elapsed between
//! the call to start and stop

class Timer  
{
public:
	//! constructor
	Timer() { reset(); }

	//! Start the timer
	void start();

	//! Stop the timer
	void stop();

	//! Reset the timer
	void reset();

	//! Get the elapsed time
	void GetTime(int& nhour, int& nmin, int& nsec);

	//! return the time as a text string
	void time_str(char* sz);

private:
	time_t	m_start;	//!< time at start
	time_t	m_stop;		//!< time at last stop

	double	m_sec;		//!< accumulated time so far in seconds
};

#endif // !defined(AFX_TIMER_H__5C4CDF72_0B19_4C5B_9B21_DB7B85FCEC4D__INCLUDED_)
