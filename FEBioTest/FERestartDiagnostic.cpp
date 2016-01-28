#include "FERestartDiagnostics.h"
#include <FECore/FEModel.h>
#include <FECore/DumpFile.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
FERestartDiagnostic::FERestartDiagnostic(FEModel*pfem) : FECoreTask(pfem)
{
	m_bok = false;
	strcpy(m_szdmp, "out.dmp");
}

//-----------------------------------------------------------------------------
bool restart_test_cb(FEModel* pfem, unsigned int nwen, void* pd)
{
	FERestartDiagnostic* ptask = (FERestartDiagnostic*) pd;

	// try to create archive
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
	FEModel& fem = *GetFEModel();

	while (fem.Solve() == false)
	{
		if (m_bok)
		{
			// reopen the dump file for readin
			DumpFile ar(fem);
			if (ar.Open(m_szdmp) == false)
			{
				felog.printf("FAILED OPENING RESTART DUMP FILE.\n");
				return false;
			}

			// reset output mode
			felog.SetMode(Logfile::FILE_AND_SCREEN);
	
			// read the model data
			fem.Serialize(ar);
			m_bok = false;
		}
		else return false;
	}

	return true;
}
