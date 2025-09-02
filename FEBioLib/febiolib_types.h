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
	double total_time = 0;
	double input_time = 0;
	double init_time = 0;
	double solve_time = 0;
	double io_time = 0;
	double total_ls_factor = 0;
	double total_ls_backsolve = 0;
	double total_reform = 0;
	double total_stiff = 0;
	double total_rhs = 0;
	double total_update = 0;
	double total_qn = 0;
	double total_serialize = 0;
	double total_callback = 0;
	double total_other = 0;
};

struct TimeStepStats {
	int iters = 0;
	int nrhs = 0;
	int refs = 0;
	int status = 0; // 0 = failed, 1 = converged
};
