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
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
// The FECoreTask class is the base class for all tasks.
// A task is simply the highest level module which defines what the code will do.
class FECORE_API FECoreTask : public FECoreBase
{
	FECORE_SUPER_CLASS(FETASK_ID)
	FECORE_BASE_CLASS(FECoreTask)

public:
	FECoreTask(FEModel* fem);
	virtual ~FECoreTask(void);

	//! initialize the task
	//! make sure to call FEModel::Init at some point
	virtual bool Init(const char* szfile) = 0;

	//! Run the task.
	virtual bool Run() = 0;
};
