#pragma once

// break point class
struct BREAK_POINT
{
	int		nwhen;		// break on callback
	double	time;		// break at time (nwhen == -1)
	int		flag;		// flag to see if break point was hit
};

// clear all break points
void clear_break_points();

// clear a particular break point
void clear_break_points(int n);

// add a break point
void add_break_point(const char* szcond);

// add a break point that breaks at (sim) time t
void add_time_break_point(double t);

// add a break point that breats at an event
void add_event_break_point(const char* szwhen);

// break points cb
void print_break_points();

// see if a break point is reached (return -1 if break point was not reached)
int check_break(int nwhen, double t);
