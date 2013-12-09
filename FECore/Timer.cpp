// Timer.cpp: implementation of the Timer class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Timer.h"
#include <stdio.h>

//-----------------------------------------------------------------------------
Timer::Timer()
{ 
	reset(); 
}

//-----------------------------------------------------------------------------
void Timer::start()
{
	time(&m_start);
	m_brunning = true;
}

//-----------------------------------------------------------------------------
void Timer::stop()
{
	time(&m_stop);
	m_brunning = false;

	m_sec += (double) difftime(m_stop, m_start);
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
	time_t pause;
	time(&pause);
	return (m_sec + (double) difftime(pause, m_start));
}

//-----------------------------------------------------------------------------
void Timer::GetTime(int& nhour, int& nmin, int& nsec)
{
	double sec = (m_brunning? peek() : m_sec);
	nhour = (int) (sec / 3600.0); sec -= nhour*3600;
	nmin  = (int) (sec /   60.0); sec -= nmin*60;
	nsec  = (int) (sec);
}

//-----------------------------------------------------------------------------
void Timer::time_str(char* sz)
{
	int nhour, nmin, nsec;

	GetTime(nhour, nmin, nsec);
	sprintf(sz, "%d:%02d:%02d", nhour, nmin, nsec);
}
