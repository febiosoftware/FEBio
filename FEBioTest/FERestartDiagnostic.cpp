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
#include "FERestartDiagnostics.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/DumpFile.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
FERestartDiagnostic::FERestartDiagnostic(FEModel*pfem) : FECoreTask(pfem), m_dmp(*pfem)
{
	m_bok = false;
	m_bfile = false;
	strcpy(m_szdmp, "out.dmp");
}

//-----------------------------------------------------------------------------
bool restart_test_cb(FEModel* pfem, unsigned int nwen, void* pd)
{
	FERestartDiagnostic* ptask = (FERestartDiagnostic*) pd;

	// try to create archive
	if (ptask->m_bfile)
	{
		DumpFile ar(*pfem);
		if (ar.Create(ptask->m_szdmp) == false)
		{
			feLogErrorEx(pfem, "FAILED CREATING RESTART DUMP FILE.\n");
			return false;
		}
		
		// serialize the data
		pfem->Serialize(ar);

		// close the dump file.
		ar.Close();
	}
	else
	{
		ptask->m_dmp.Open(true, false);
		pfem->Serialize(ptask->m_dmp);
	}

	// set the ok flag to tell the diagnostic that the restart 
	// file was created successfully.
	ptask->m_bok = true;

	// suppress the error output
	pfem->BlockLog();

	return false;
}

//-----------------------------------------------------------------------------
// initialize the diagnostic
bool FERestartDiagnostic::Init(const char* sz)
{
	FEModel& fem = *GetFEModel();

	// copy the file name (if any)
	if (sz && (sz[0] != 0)) strcpy(m_szdmp, sz);

	// Make sure that restart flag is off.
	// This is because we are hijacking restart and we don't
	// want regular restart to interfere.
	// TODO: is this really necessary? This creates a circular link between FEBioTest and FEBioLib
//	fem.SetDumpLevel(FE_DUMP_NEVER);

	// Add the restart callback
	fem.AddCallback(restart_test_cb, CB_MAJOR_ITERS, this);

	// do the FE initialization
	return fem.Init();
}

//-----------------------------------------------------------------------------
// run the diagnostic
bool FERestartDiagnostic::Run()
{
	FEModel* fem = GetFEModel();

	while (fem->Solve() == false)
	{
		// reset output mode (it was turned off in restart_test_cb)
		fem->UnBlockLog();

		if (m_bok)
		{
			feLogEx(fem, "Reading restart file...");
			if (m_bfile)
			{
				// reopen the dump file for readin
				DumpFile ar(*fem);
				if (ar.Open(m_szdmp) == false)
				{
					feLogErrorEx(fem, "FAILED OPENING RESTART DUMP FILE.\n");
					return false;
				}

				fem->Serialize(ar);
			}
			else
			{
				m_dmp.Open(false, false);
				fem->Serialize(m_dmp);
			}
			m_bok = false;
			feLogEx(fem, "done!\n");

			// reopen the log file for appending
/*			const char* szlog = fem.GetLogfileName();
			if (felog.append(szlog) == false)
			{
				printf("WARNING: Could not reopen log file. A new log file is created\n");
				felog.open(szlog);
			}
*/
		}
		else return false;
	}

	return true;
}
