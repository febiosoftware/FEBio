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

struct ModelStats {
	int		ntimeSteps = 0;		//!< total nr of time steps
	int		ntotalIters = 0;	//!< total nr of equilibrium iterations
	int		ntotalRHS = 0;		//!< total nr of right hand side evaluations
	int		ntotalReforms = 0;	//!< total nr of stiffness reformations
};

struct TimingInfo {
	double total_time;
	double input_time;
	double init_time;
	double solve_time;
	double io_time;
	double total_ls_factor;
	double total_ls_backsolve;
	double total_reform;
	double total_stiff;
	double total_rhs;
	double total_update;
	double total_qn;
	double total_serialize;
	double total_callback;
	double total_other;
};

