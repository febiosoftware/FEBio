#include "stdafx.h"
#include "FEBioModel.h"
#include "FEBioPlot/FEBioPlotFile.h"
#include "FEBioXML/FEBioImport.h"
#include "FEBioXML/FERestartImport.h"
#include "FECore/NodeDataRecord.h"
#include "FECore/ElementDataRecord.h"
#include "FECore/ObjectDataRecord.h"
#include "FECore/NLConstraintDataRecord.h"
#include "FECore/FERigidBody.h"
#include "FECore/log.h"
#include "FECore/FECoreKernel.h"
#include "FECore/DumpFile.h"
#include "FECore/DOFS.h"
#include "version.h"

//-----------------------------------------------------------------------------
// echo the input data to the log file
extern void echo_input(FEBioModel& fem);

//-----------------------------------------------------------------------------
// Callback that guides FEBio output
bool output_cb(FEModel* pfem, unsigned int nwhen, void* pd)
{
	FEBioModel* pfebio = (FEBioModel*) pd;
	pfebio->Write(nwhen);
	return true;
}

//-----------------------------------------------------------------------------
// Constructor of FEBioModel class.
FEBioModel::FEBioModel()
{
	// --- I/O-Data ---
	m_szfile_title = 0;
	m_szfile[0] = 0;
	m_szplot[0] = 0;
	m_szlog[0] = 0;
	m_szdump[0] = 0;
	m_debug = false;
	m_becho = true;
	m_plot = 0;

	// Add the output callback
	AddCallback(output_cb, CB_ALWAYS, this);
}

//-----------------------------------------------------------------------------
FEBioModel::~FEBioModel()
{
	// close the plot file
	if (m_plot) { delete m_plot; m_plot = 0; }
}

//-----------------------------------------------------------------------------
Timer& FEBioModel::GetSolveTimer()
{
	return m_SolveTime;
}

//=============================================================================
//
//		FEBioModel: I-O Functions
//
//=============================================================================

//-----------------------------------------------------------------------------
//! Return the data store
DataStore& FEBioModel::GetDataStore()
{
	return m_Data;
}

//-----------------------------------------------------------------------------
//! Add a data record to the data store
void FEBioModel::AddDataRecord(DataRecord* pd)
{
	m_Data.AddRecord(pd); 
}

//-----------------------------------------------------------------------------
//! Get the plot file
PlotFile* FEBioModel::GetPlotFile()
{
	return m_plot;
}

//-----------------------------------------------------------------------------
//! Sets the name of the FEBio input file
void FEBioModel::SetInputFilename(const char* szfile)
{ 
	strcpy(m_szfile, szfile); 
	m_szfile_title = strrchr(m_szfile, '/');
	if (m_szfile_title == 0) 
	{
		m_szfile_title = strrchr(m_szfile, '\\'); 
		if (m_szfile_title == 0) m_szfile_title = m_szfile; else ++m_szfile_title;
	}
	else ++m_szfile_title;
}

//-----------------------------------------------------------------------------
//! Set the name of the log file
void FEBioModel::SetLogFilename(const char* szfile) 
{ 
	strcpy(m_szlog , szfile); 
}

//-----------------------------------------------------------------------------
//! Set the name of the plot file
void FEBioModel::SetPlotFilename(const char* szfile) 
{ 
	strcpy(m_szplot, szfile); 
}

//-----------------------------------------------------------------------------
//! Set the name of the restart archive (i.e. the dump file)
void FEBioModel::SetDumpFilename (const char* szfile) 
{ 
	strcpy(m_szdump, szfile); 
}

//-----------------------------------------------------------------------------
//! Return the name of the input file
const char* FEBioModel::GetInputFileName()
{ 
	return m_szfile; 
}

//-----------------------------------------------------------------------------
//! Return the name of the log file
const char* FEBioModel::GetLogfileName()
{ 
	return m_szlog;  
}

//-----------------------------------------------------------------------------
//! Return the name of the plot file
const char* FEBioModel::GetPlotFileName()
{ 
	return m_szplot; 
}

//-----------------------------------------------------------------------------
//! Return the dump file name.
const char* FEBioModel::GetDumpFileName()
{
	return	m_szdump;
}

//-----------------------------------------------------------------------------
//! get the file title (i.e. name of input file without the path)
//! \todo Do I actually need to store this?
const char* FEBioModel::GetFileTitle()
{ 
	return m_szfile_title; 
}

//=============================================================================
//    I N P U T
//=============================================================================

//-----------------------------------------------------------------------------
//! This routine reads in an input file and performs some initialization stuff.
//! The rest of the initialization is done in Init

bool FEBioModel::Input(const char* szfile)
{
	// start the timer
	TimerTracker t(m_InputTime);

	// create file reader
	FEBioImport fim;

	// Load the file
	if (fim.Load(*this, szfile) == false)
	{
		char szerr[256];
		fim.GetErrorMessage(szerr);
		fprintf(stderr, szerr);

		return false;
	}

	// set the input file name
	SetInputFilename(szfile);

	// see if user redefined output filenames
	if (fim.m_szdmp[0]) SetDumpFilename(fim.m_szdmp);
	if (fim.m_szlog[0]) SetLogFilename (fim.m_szlog);
	if (fim.m_szplt[0]) SetPlotFilename(fim.m_szplt);

	// set the plot file
	FEBioPlotFile* pplt = new FEBioPlotFile(*this);
	m_plot = pplt;

	// set compression
	pplt->SetCompression(fim.m_nplot_compression);

	// define the plot file variables
	FEMesh& mesh = GetMesh();
	int NP = (int) fim.m_plot.size();
	for (int i=0; i<NP; ++i)
	{
		FEBioImport::PlotVariable& var = fim.m_plot[i];

		vector<int> item = var.m_item;
		if (item.empty() == false)
		{
			// TODO: currently, this is only supported for domain variables, where
			//       the list is a list of materials
			vector<int> lmat = var.m_item;

			// convert the material list to a domain list
			mesh.DomainListFromMaterial(lmat, item);
		}

		// add the plot output variable
		if (pplt->AddVariable(var.m_szvar, item, var.m_szdom) == false)
		{
			felog.printf("FATAL ERROR: Output variable \"%s\" is not defined\n", var.m_szvar);
			return false;
		}
	}

	// add the data records
	int ND = fim.m_data.size();
	for (int i=0; i<ND; ++i) AddDataRecord(fim.m_data[i]);

	// we're done reading
	return true;
}

//=============================================================================
//    O U T P U T
//=============================================================================

//-----------------------------------------------------------------------------
//! Export state to plot file.
void FEBioModel::Write(unsigned int nwhen)
{
	TimerTracker t(m_IOTimer);

	// get the current step
	FEAnalysis* pstep = GetCurrentStep();

	if (m_plot)
	{
		// get the plot level
		int nplt = pstep->GetPlotLevel();

		// if we don't want to plot anything we return
		if (nplt != FE_PLOT_NEVER)
		{
			// try to open the plot file
			if (nwhen == CB_STEP_ACTIVE)
			{
				if (m_plot->IsValid() == false)
				{
					if (m_plot->Open(*this, m_szplot) == false)
					{
						felog.printf("ERROR : Failed creating PLOT database\n");
						delete m_plot;
						m_plot = 0;
					}

					// Since it is assumed that for the first timestep
					// there are no loads or initial displacements, the case n=0 is skipped.
					// Therefor we can output those results here.
					// TODO: Offcourse we should actually check if this is indeed
					//       the case, otherwise we should also solve for t=0
					if (m_plot) m_plot->Write(*this);
				}
			}
			else
			{
				// assume we won't be writing anything
				bool bout = false;

				// see if we need to output something
				bool bdebug = GetDebugFlag();

				// when debugging we always output
				// (this coule mean we may end up writing the same state multiple times)
				if (bdebug) bout = true;
				else
				{
					switch (nwhen)
					{
					case CB_MINOR_ITERS: if (nplt == FE_PLOT_MINOR_ITRS   ) bout = true; break;
					case CB_MAJOR_ITERS  : 
						if ((nplt == FE_PLOT_MAJOR_ITRS ) && (pstep->m_ntimesteps % pstep->m_nplot_stride == 0)) bout = true; 
						if ((nplt == FE_PLOT_MUST_POINTS) && (pstep->m_nmust >= 0)) bout = true;
						break;
					case CB_AUGMENT: if (nplt == FE_PLOT_AUGMENTATIONS) bout = true; break;
					case CB_SOLVED : if (nplt == FE_PLOT_FINAL) bout = true;
					}
				}

				// output the state if requested
				if (bout) 
				{
					if (m_plot) m_plot->Write(*this);
				}
			}
		}
	}

	// Dump converged state to the archive
	int ndump = pstep->GetDumpLevel();
	if (ndump != FE_DUMP_NEVER)
	{
		bool bdump = false;
		if ((nwhen == CB_SOLVED     ) && (ndump == FE_DUMP_STEP      )) bdump = true;
		if ((nwhen == CB_MAJOR_ITERS) && (ndump == FE_DUMP_MAJOR_ITRS)) bdump = true;
		if (bdump) DumpData();
	}

	// write the output data
	int nout = pstep->GetOutputLevel();
	if (nout != FE_OUTPUT_NEVER)
	{
		bool bout = false;
		switch (nwhen)
		{
		case CB_MINOR_ITERS: if (nout == FE_OUTPUT_MINOR_ITRS) bout = true; break;
		case CB_MAJOR_ITERS:
			if (nout == FE_OUTPUT_MAJOR_ITRS) bout = true;
			if ((nout == FE_OUTPUT_MUST_POINTS) && (pstep->m_nmust >= 0)) bout = true;
			break;
		case CB_SOLVED:
			if (nout == FE_OUTPUT_FINAL) bout = true;
			break;
		}

		if (bout) WriteData();
	}
}

//-----------------------------------------------------------------------------
//! Write user data to the logfile
void FEBioModel::WriteData()
{
	m_Data.Write();
}

//-----------------------------------------------------------------------------
//! Dump state to archive for restarts
void FEBioModel::DumpData()
{
	DumpFile ar(*this);
	if (ar.Create(m_szdump) == false)
	{
		felog.printf("WARNING: Failed creating restart file (%s).\n", m_szdump);
	}
	else 
	{
		Serialize(ar);
		felog.printf("\nRestart point created. Archive name is %s\n", m_szdump);
	}
}

//=============================================================================
//    R E S T A R T
//=============================================================================

//-----------------------------------------------------------------------------
//!  Reads or writes the current state to/from a binary file
//!  This is used to restart the solution from a saved position
//!  or to create a restart point.
//!  A version number is written to file to make sure the same
//!  format is used for reading and writing.
//! \param[in] ar the archive to which the data is serialized
//! \sa DumpFile
bool FEBioModel::Serialize(DumpStream& ar)
{
	// don't need to do anything for running restarts
	if (ar.IsShallow())
	{
		// serialize model data
		FEModel::Serialize(ar);
	}
	else
	{
		if (ar.IsSaving())
		{
			// --- version number ---
			ar << (int) RSTRTVERSION;
		}
		else
		{
			// --- version ---
			int nversion;
			ar >> nversion;

			// make sure it is the right version
			if (nversion != RSTRTVERSION) return false;
		}

		// serialize model data
		FEModel::Serialize(ar);

		// serialize data store
		SerializeDataStore(ar);

		// --- Save IO Data
		SerializeIOData(ar);
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Serialization of FEBioModel data
void FEBioModel::SerializeIOData(DumpStream &ar)
{
	if (ar.IsSaving())
	{
		// file names
		ar << m_szfile << m_szplot << m_szlog << m_szdump;
		ar << m_sztitle;

		// plot file
		int npltfmt = 2;
		ar << npltfmt;

		// data records
		SerializeDataStore(ar);
	}
	else
	{
		// file names
		ar >> m_szfile >> m_szplot >> m_szlog >> m_szdump;
		ar >> m_sztitle;

		// don't forget to call store the input file name so
		// that m_szfile_title gets initialized
		SetInputFilename(m_szfile);

		// get the plot file format (should be 2)
		int npltfmt = 0;
		ar >> npltfmt;
		assert(npltfmt == 2);

		// remove the plot file (if any)
		if (m_plot) { delete m_plot; m_plot = 0; }

		// create the plot file and open it for appending
		m_plot = new FEBioPlotFile(*this);
		if (m_plot->Append(*this, m_szplot) == false)
		{
			printf("FATAL ERROR: Failed reopening plot database %s\n", m_szplot);
			throw "FATAL ERROR";
		}

		// data records
		SerializeDataStore(ar);
	}
}

//-----------------------------------------------------------------------------
void FEBioModel::SerializeDataStore(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		int N = m_Data.Size();
		ar << N;
		for (int i=0; i<N; ++i)
		{
			DataRecord* pd = m_Data.GetDataRecord(i);

			int ntype = -1;
			if (dynamic_cast<NodeDataRecord*       >(pd)) ntype = FE_DATA_NODE;
			if (dynamic_cast<ElementDataRecord*    >(pd)) ntype = FE_DATA_ELEM;
			if (dynamic_cast<ObjectDataRecord*     >(pd)) ntype = FE_DATA_RB;
			if (dynamic_cast<NLConstraintDataRecord*>(pd)) ntype = FE_DATA_NLC;
			assert(ntype != -1);
			ar << ntype;
			pd->Serialize(ar);
		}
	}
	else
	{
		int N;
		m_Data.Clear();
		ar >> N;
		for (int i=0; i<N; ++i)
		{
			int ntype;
			ar >> ntype;

			DataRecord* pd = 0;
			switch(ntype)
			{
			case FE_DATA_NODE: pd = new NodeDataRecord        (this, 0); break;
			case FE_DATA_ELEM: pd = new ElementDataRecord     (this, 0); break;
			case FE_DATA_RB  : pd = new ObjectDataRecord      (this, 0); break;
			case FE_DATA_NLC : pd = new NLConstraintDataRecord(this, 0); break;
			}
			assert(pd);
			pd->Serialize(ar);
			m_Data.AddRecord(pd);
		}
	}
}

//=============================================================================
//    I N I T I A L I Z A T I O N
//=============================================================================

//-----------------------------------------------------------------------------
// Forward declarations
int Hello();

//-----------------------------------------------------------------------------
//! This function performs one-time-initialization stuff. All the different 
//! modules are initialized here as well. This routine also performs some
//! data checks

bool FEBioModel::Init()
{
	TimerTracker t(m_InitTime);

	// Open the logfile
	if (InitLogFile() == false) return false;

	// open plot database file
	if (m_pStep->GetPlotLevel() != FE_PLOT_NEVER)
	{
		if (m_plot == 0) 
		{
			m_plot = new FEBioPlotFile(*this);
		}

		// see if a valid plot file name is defined.
		const char* szplt = GetPlotFileName();
		if (szplt[0] == 0)
		{
			// if not, we take the input file name and set the extension to .xplt
			char sz[1024] = {0};
			strcpy(sz, GetInputFileName());
			char *ch = strrchr(sz, '.');
			if (ch) *ch = 0;
			strcat(sz, ".xplt");
			SetPlotFilename(sz);
		}
	}

	// initialize model data
	if (FEModel::Init() == false) 
	{
		felog.printf("\nFATAL ERROR: Model initialization failed\n");
		const char* szerr = fecore_get_error_string();
		if (szerr) felog.printf("REASON: %s\n\n", szerr);
		return false;
	}

	// see if a valid dump file name is defined.
	const char* szdmp = this->GetDumpFileName();
	if (szdmp[0] == 0)
	{
		// if not, we take the input file name and set the extension to .dmp
		char sz[1024] = {0};
		strcpy(sz, GetInputFileName());
		char *ch = strrchr(sz, '.');
		if (ch) *ch = 0;
		strcat(sz, ".dmp");
		SetDumpFilename(sz);
	}

	// echo fem data to the logfile
	// we do this here (and not e.g. directly after input)
	// since the data can be changed after input, which is the case,
	// for instance, in the parameter optimization module
	if (m_becho) echo_input(*this);

	// Alright, all initialization is done, so let's get busy !
	return true;
}

//-----------------------------------------------------------------------------
// Opens the log file
bool FEBioModel::InitLogFile()
{
	// Only do this if the log file is not valid
	if (!felog.is_valid()) 
	{
		// see if a valid log file name is defined.
		const char* szlog = GetLogfileName();
		if (szlog[0] == 0)
		{
			// if not, we take the input file name and set the extension to .log
			char sz[1024] = {0};
			strcpy(sz, GetInputFileName());
			char *ch = strrchr(sz, '.');
			if (ch) *ch = 0;
			strcat(sz, ".log");
			SetLogFilename(sz);
		}
		
		if (felog.open(m_szlog) == false)
		{
			felog.printbox("FATAL ERROR", "Failed creating log file");
			return false;
		}

		// make sure we have a step
		if (m_pStep == 0)
		{
			felog.printf("FATAL ERROR: No step defined\n\n");
			return false;
		}

		// if we don't want to output anything we only output to the logfile
		if (m_pStep->GetPrintLevel() == FE_PRINT_NEVER) felog.SetMode(Logfile::FILE_ONLY);

		// print welcome message to file
		Logfile::MODE m = felog.SetMode(Logfile::FILE_ONLY);
		Hello();
		felog.SetMode(m);
	}

	return true;
}

//-----------------------------------------------------------------------------
//! This function resets the FEM data so that a new run can be done.
//! This routine is called from the optimization routine.

bool FEBioModel::Reset()
{
	// Reset model data
	FEModel::Reset();

	// re-initialize the log file
	if (InitLogFile() == false) return false;

	// open plot database file
	if (m_pStep->GetPlotLevel() != FE_PLOT_NEVER)
	{
		if (m_plot == 0) 
		{
			m_plot = new FEBioPlotFile(*this);
		}

		if (m_plot->Open(*this, m_szplot) == false)
		{
			felog.printf("ERROR : Failed creating PLOT database\n");
			return false;
		}
	}

	// do the callback
	DoCallback(CB_INIT);

	// All data is reset successfully
	return true;
}

//=============================================================================
//                               S O L V E
//=============================================================================

//-----------------------------------------------------------------------------
//! This is the main solve method. This function loops over all analysis steps
//! and solves each one in turn. 
//! \sa FEAnalysis

bool FEBioModel::Solve()
{
	// start the total time tracker
	m_SolveTime.start();

	// solve the FE model
	bool bconv = FEModel::Solve();

	// stop total time tracker
	m_SolveTime.stop();

	// get and print elapsed time
	if (felog.GetMode() != Logfile::NEVER)
	{
		char sztime[64];
		Logfile::MODE old_mode = felog.SetMode(Logfile::SCREEN_ONLY);
		m_SolveTime.time_str(sztime);
		felog.printf("\n Elapsed time : %s\n\n", sztime);

		// print more detailed timing info to the log file
		felog.SetMode(Logfile::FILE_ONLY);

		// sum up all the times spend in the linear solvers
		double total_time = 0.0;
		double input_time   = m_InputTime.GetTime(); total_time += input_time;
		double init_time    = m_InitTime.GetTime (); total_time += init_time;
		double solve_time   = m_SolveTime.GetTime(); total_time += solve_time;
		double io_time      = m_IOTimer.GetTime  ();
		double total_linsol = 0.0;
		double total_reform = 0.0;
		double total_stiff  = 0.0;
		double total_rhs    = 0.0;
		double total_update = 0.0;
		int NS = Steps();
		for (int i = 0; i<NS; ++i)
		{
			FEAnalysis* pstep = GetStep(i);
			FESolver* psolve = pstep->GetFESolver();
			if (psolve) 
			{
				total_linsol += psolve->m_SolverTime.GetTime();
				total_reform += psolve->m_ReformTime.GetTime();
				total_stiff  += psolve->m_StiffnessTime.GetTime();
				total_rhs    += psolve->m_RHSTime.GetTime();
				total_update += psolve->m_UpdateTime.GetTime();
			}
		}


		felog.printf(" T I M I N G   I N F O R M A T I O N\n\n");
		Timer::time_str(input_time  , sztime); felog.printf("\tInput time ...................... : %s (%lg sec)\n\n", sztime, input_time  );
		Timer::time_str(init_time   , sztime); felog.printf("\tInitialization time ............. : %s (%lg sec)\n\n", sztime, init_time   );
		Timer::time_str(solve_time  , sztime); felog.printf("\tSolve time ...................... : %s (%lg sec)\n\n", sztime, solve_time  );
		Timer::time_str(io_time     , sztime); felog.printf("\t   IO-time (plot, dmp, data) .... : %s (%lg sec)\n\n", sztime, io_time     );
		Timer::time_str(total_reform, sztime); felog.printf("\t   reforming stiffness .......... : %s (%lg sec)\n\n", sztime, total_reform);
		Timer::time_str(total_stiff , sztime); felog.printf("\t   evaluating stiffness ......... : %s (%lg sec)\n\n", sztime, total_stiff );
		Timer::time_str(total_rhs   , sztime); felog.printf("\t   evaluating residual .......... : %s (%lg sec)\n\n", sztime, total_rhs   );
		Timer::time_str(total_update, sztime); felog.printf("\t   model update ................. : %s (%lg sec)\n\n", sztime, total_update);
		Timer::time_str(total_linsol, sztime); felog.printf("\t   time in linear solver ........ : %s (%lg sec)\n\n", sztime, total_linsol);
		Timer::time_str(total_time  , sztime); felog.printf("\tTotal elapsed time .............. : %s (%lg sec)\n\n", sztime, total_time  );


		felog.SetMode(old_mode);

		if (bconv)
		{
			felog.printf("\n N O R M A L   T E R M I N A T I O N\n\n");
		}
		else
		{
			felog.printf("\n E R R O R   T E R M I N A T I O N\n\n");
		}

		// flush the log file
		felog.flush();
	}

	// close the log file
	felog.close();

	// close the plot file
	if (m_plot) m_plot->Close();

	// We're done !
	return bconv;
}
