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
#include "FEOptimize.h"
#include "FEOptimizeData.h"
#include "FECore/FECoreKernel.h"
#include "FECore/log.h"

//-----------------------------------------------------------------------------
#define VERSION 2
#define SUB_VERSION	0

//-----------------------------------------------------------------------------
FEOptimize::FEOptimize(FEModel* pfem) : FECoreTask(pfem), m_opt(pfem)
{

}

//-----------------------------------------------------------------------------
bool FEOptimize::Init(const char* szfile)
{
	char szversion[32] = { 0 };
	sprintf(szversion, "version %d.%d", VERSION, SUB_VERSION);
	feLog("P A R A M E T E R   O P T I M I Z A T I O N   M O D U L E\n%s\n\n", szversion);

	// read the data from the xml input file
	if (m_opt.Input(szfile) == false) return false;

	// do initialization
	bool ret = m_opt.Init();

	if (ret == false)
	{
		feLogErrorEx(m_opt.GetFEModel(), "Failed to initialize the optimization data.");
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEOptimize::Run()
{
	// solve the problem
	bool bret = m_opt.Solve();

	if (bret)
		feLog("\n\n N O R M A L   T E R M I N A T I O N\n\n");
	else 
		feLog("\n\n E R R O R   T E R M I N A T I O N\n\n");

	return bret;
}
