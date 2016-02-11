#include "FERestartDiagnostics.h"
#include <FEBioLib/FEBioModel.h>
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
			felog.printf("FAILED CREATING RESTART DUMP FILE.\n");
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
	felog.SetMode(Logfile::NEVER);

	return false;
}

//-----------------------------------------------------------------------------
// initialize the diagnostic
bool FERestartDiagnostic::Init(const char* sz)
{
	FEModel& fem = *GetFEModel();

	// copy the file name (if any)
	if (sz && (sz[0] != 0)) strcpy(m_szdmp, sz);

	// Make sure that restart flags are off for all steps
	// This is because we are hijacking restart and we don't
	// want regular restart to interfere.
	int NS = fem.Steps();
	for (int i=0; i<NS; ++i)
	{
		FEAnalysis* ps = fem.GetStep(i);
		ps->m_ndump = FE_DUMP_NEVER;
	}

	// Add the restart callback
	fem.AddCallback(restart_test_cb, CB_MAJOR_ITERS, this);

	// do the FE initialization
	return fem.Init();
}

//-----------------------------------------------------------------------------
// run the diagnostic
bool FERestartDiagnostic::Run()
{
	FEBioModel& fem = dynamic_cast<FEBioModel&>(*GetFEModel());

	while (fem.Solve() == false)
	{
		if (m_bok)
		{
			if (m_bfile)
			{
				// reopen the dump file for readin
				DumpFile ar(fem);
				if (ar.Open(m_szdmp) == false)
				{
					felog.printf("FAILED OPENING RESTART DUMP FILE.\n");
					return false;
				}

				fem.Serialize(ar);
			}
			else
			{
				m_dmp.Open(false, false);
				fem.Serialize(m_dmp);
			}
			m_bok = false;

			// reopen the log file for appending
			const char* szlog = fem.GetLogfileName();
			if (felog.append(szlog) == false)
			{
				printf("WARNING: Could not reopen log file. A new log file is created\n");
				felog.open(szlog);
			}

			// reset output mode
			felog.SetMode(Logfile::FILE_AND_SCREEN);
		}
		else return false;
	}

	return true;
}
