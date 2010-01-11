#include "stdafx.h"
#include "FEMemoryDiagnostic.h"
#include "log.h"
#include "console.h"

FEMemoryDiagnostic::FEMemoryDiagnostic(FEM& fem) : FEDiagnostic(fem)
{
	m_szfile[0] = 0;
	m_iters = 1;
}

FEMemoryDiagnostic::~FEMemoryDiagnostic(void)
{

}

bool FEMemoryDiagnostic::Init()
{
	// try to open the file
	FEFEBioImport im;
	if (im.Load(m_fem, m_szfile) == false)
	{
		return false;
	}

	// make sure the iters is a positive number
	if (m_iters <= 0) return false;

	// turn off all output
	m_fem.m_pStep->SetPlotLevel(FE_PLOT_NEVER);
	m_fem.m_pStep->SetPrintLevel(FE_PRINT_NEVER);

	return true;
}

bool FEMemoryDiagnostic::Run()
{
	// the logfile is a shared resource between the master FEM and the RVE
	// in order not to corrupt the logfile we don't print anything for
	// the RVE problem.
	Logfile& log = GetLogfile();
	Logfile::MODE nmode = log.GetMode();
	GetLogfile().SetMode(Logfile::NEVER);

	// turn the console off
	Console::GetHandle()->Deactivate();

	// run the problem
	for (int i=0; i<m_iters; ++i)
	{
		fprintf(stderr, "%d/%d: ...", i+1, m_iters);
		m_fem.Reset();
		bool b = m_fem.Solve();
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
