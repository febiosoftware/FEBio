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
#include "breakpoint.h"
#include <FECore/Callback.h>
#include <vector>
#include <iostream>
using namespace std;

//-----------------------------------------------------------------------------
// list of current break points
std::vector<BREAK_POINT> break_points;

//-----------------------------------------------------------------------------
void clear_break_points()
{
	break_points.clear();
}

//-----------------------------------------------------------------------------
void add_time_break_point(double t)
{
	BREAK_POINT bp;
	bp.nwhen = -1;
	bp.time = t;
	bp.flag = 1;
	break_points.push_back(bp);
}

//-----------------------------------------------------------------------------
void add_event_break_point(int nwhen)
{
	BREAK_POINT bp;
	bp.nwhen = nwhen;
	bp.time = 0;
	bp.flag = 1;
	break_points.push_back(bp);
}

//-----------------------------------------------------------------------------
#ifdef WIN32
#define szcmp	_stricmp
#else
#define szcmp	strcasecmp
#endif

//-----------------------------------------------------------------------------
void add_break_point(const char* szcond)
{
	if      (szcmp(szcond, "ALWAYS"       ) == 0) add_event_break_point(CB_ALWAYS);
	else if (szcmp(szcond, "INIT"         ) == 0) add_event_break_point(CB_INIT);
	else if (szcmp(szcond, "STEP_ACTIVE"  ) == 0) add_event_break_point(CB_STEP_ACTIVE);
	else if (szcmp(szcond, "MAJOR_ITERS"  ) == 0) add_event_break_point(CB_MAJOR_ITERS);
	else if (szcmp(szcond, "MINOR_ITERS"  ) == 0) add_event_break_point(CB_MINOR_ITERS);
	else if (szcmp(szcond, "SOLVED"       ) == 0) add_event_break_point(CB_SOLVED);
	else if (szcmp(szcond, "UPDATE_TIME"  ) == 0) add_event_break_point(CB_UPDATE_TIME);
	else if (szcmp(szcond, "AUGMENT"      ) == 0) add_event_break_point(CB_AUGMENT);
	else if (szcmp(szcond, "STEP_SOLVED"  ) == 0) add_event_break_point(CB_STEP_SOLVED);
	else if (szcmp(szcond, "REFORM"       ) == 0) add_event_break_point(CB_MATRIX_REFORM);
	else if (szcmp(szcond, "MATRIX_SOLVE" ) == 0) add_event_break_point(CB_PRE_MATRIX_SOLVE);
	else
	{
		double f = atof(szcond);
		add_time_break_point(f);
	}
}

//-----------------------------------------------------------------------------
void clear_break_points(int n)
{
	if (n == -1)
	{
		break_points.clear();
	}
	else
	{
		if ((n >= 0) && (n < break_points.size())) break_points.erase(break_points.begin() + n);
	}
}

//-----------------------------------------------------------------------------
// break points cb
void print_break_points()
{
	int nbp = (int)break_points.size();
	if (nbp == 0)
	{
		cout << "No breakpoints defined.\n";
		return;
	}

	for (int i = 0; i < nbp; ++i)
	{
		BREAK_POINT& bpi = break_points[i];
		cout << i + 1 << ": ";
		if (bpi.nwhen > 0)
		{
			cout << "event = ";
			switch (bpi.nwhen)
			{
			case CB_ALWAYS: cout << "ALWAYS"; break;
			case CB_INIT: cout << "INIT"; break;
			case CB_STEP_ACTIVE: cout << "STEP ACTIVE"; break;
			case CB_MAJOR_ITERS: cout << "MAJOR_ITERS"; break;
			case CB_MINOR_ITERS: cout << "MINOR_ITERS"; break;
			case CB_SOLVED: cout << "SOLVED"; break;
			case CB_UPDATE_TIME: cout << "UPDATE_TIME"; break;
			case CB_AUGMENT: cout << "AUGMENT"; break;
			case CB_STEP_SOLVED: cout << "STEP_SOLVED"; break;
			case CB_MATRIX_REFORM: cout << "MATRIX_REFORM"; break;
			default:
				cout << "(unknown)"; break;
			}
			cout << endl;
		}
		else
		{
			cout << "time = " << bpi.time << endl;
		}
	}
}

// see if a break point is reached
int check_break(int nwhen, double t)
{
	const double eps = 1e-12;

	int nbp = (int)break_points.size();
	for (int i = 0; i<nbp; ++i)
	{
		BREAK_POINT& bpi = break_points[i];
		if (bpi.nwhen > 0)
		{
			if ((int)nwhen & bpi.nwhen)
			{
				return i;
			}
		}
		else if (nwhen == CB_MAJOR_ITERS)
		{
			if (bpi.flag && (bpi.time <= t + eps))
			{
				bpi.flag = 0;
				return i;
			}
		}
	}

	return -1;
}
