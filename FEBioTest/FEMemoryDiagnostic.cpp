#include "FEMemoryDiagnostic.h"
#include "FECore/log.h"

FEMemoryDiagnostic::FEMemoryDiagnostic(FEModel& fem) : FEDiagnostic(fem)
{
	m_szfile[0] = 0;
	m_iters = 1;

	FEAnalysis* pstep = new FEAnalysis(&fem);
	fem.AddStep(pstep);
	fem.SetCurrentStep(pstep);
}

FEMemoryDiagnostic::~FEMemoryDiagnostic(void)
{

}

bool FEMemoryDiagnostic::Init()
{
	FEModel& fem = GetFEModel();

	// try to open the file
	FEBioImport im;
	if (im.Load(fem, m_szfile) == false)
	{
		return false;
	}

	// make sure the iters is a positive number
	if (m_iters <= 0) return false;

	// turn off all output
	fem.GetCurrentStep()->SetPlotLevel(FE_PLOT_NEVER);
	fem.GetCurrentStep()->SetPrintLevel(FE_PRINT_NEVER);

	return true;
}

bool FEMemoryDiagnostic::Run()
{
	// the logfile is a shared resource between the master FEM and the RVE
	// in order not to corrupt the logfile we don't print anything for
	// the RVE problem.
	Logfile::MODE nmode = felog.GetMode();
	felog.SetMode(Logfile::LOG_NEVER);

	// run the problem
	FEModel& fem = GetFEModel();
	for (int i=0; i<m_iters; ++i)
	{
		fprintf(stderr, "%d/%d: ...", i+1, m_iters);
		fem.Reset();
		bool b = fem.Solve();
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
