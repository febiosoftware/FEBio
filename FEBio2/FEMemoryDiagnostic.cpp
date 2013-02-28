#include "stdafx.h"
#include "FEMemoryDiagnostic.h"
#include "FECore/log.h"
#include "console.h"

FEMemoryDiagnostic::FEMemoryDiagnostic(FEModel& fem) : FEDiagnostic(fem)
{
	m_szfile[0] = 0;
	m_iters = 1;
}

FEMemoryDiagnostic::~FEMemoryDiagnostic(void)
{

}

bool FEMemoryDiagnostic::Init()
{
	FEModel& fem = m_fem;

	// try to open the file
	FEFEBioImport im;
	if (im.Load(fem, m_szfile) == false)
	{
		return false;
	}

	// make sure the iters is a positive number
	if (m_iters <= 0) return false;

	// turn off all output
	m_fem.GetCurrentStep()->SetPlotLevel(FE_PLOT_NEVER);
	m_fem.GetCurrentStep()->SetPrintLevel(FE_PRINT_NEVER);

	return true;
}

bool FEMemoryDiagnostic::Run()
{
	// the logfile is a shared resource between the master FEM and the RVE
	// in order not to corrupt the logfile we don't print anything for
	// the RVE problem.
	Logfile::MODE nmode = clog.GetMode();
	clog.SetMode(Logfile::NEVER);

	// turn the console off
	Console::GetHandle()->Deactivate();

	// run the problem
	for (int i=0; i<m_iters; ++i)
	{
		fprintf(stderr, "%d/%d: ...", i+1, m_iters);
		m_fem.Reset();
		bool b = m_fem.Solve();
//		system("ps -C febio.test -o rss,vsize h");
		fprintf(stderr, "%s\n", (b?"NT" : "ET"));
	}

	return true;
}

bool FEMemoryDiagnostic::ParseSection(XMLTag &tag)
{
	if (tag == "file") tag.value(m_szfile);
	else if (tag == "iters") tag.value(m_iters);
	else return false;
	return true;
}
